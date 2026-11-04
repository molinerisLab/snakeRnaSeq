"""Restore target taxa removed by second-stage feature filtering.

Snakemake input:  original, filtered, targets, resolution
Snakemake output: restored, report
Optional params: n_metadata_cols (default: 3)

- `original` is the selected-sample count matrix BEFORE feature_filter.
- `filtered` is the selected-sample count matrix AFTER feature_filter.
- `targets` is qc/target_taxon_matches.tsv from filter_taxa.py.
- `resolution` is qc/target_taxon_resolution.tsv from filter_taxa.py.

Only targets with rescue_enabled=True in the resolution report are protected.
The output contains the union of feature-filter survivors and protected taxa,
ordered as in the original, with UNMODIFIED pre-filter counts.
An organism absent from `original` cannot be recovered by this script.

Call --help for optional standalone command-line execution.
"""

import argparse
import csv
import math
from collections import defaultdict
from pathlib import Path


def fail(message):
    raise ValueError(message)


def snakemake_param(smk, name, default=None):
    """Access a named Snakemake parameter (Namedlist has no dict.get)."""
    try:
        return smk.params[name]
    except (AttributeError, IndexError, KeyError, TypeError):
        return default


def as_taxid(value, where):
    v = str(value).strip()
    if not v or not v.isdecimal():
        fail(f"Invalid taxonomy_id {value!r} in {where}; expected a nonnegative integer")
    number = int(v)
    if number < 0:
        fail(f"Invalid taxonomy_id {value!r} in {where}; expected >=0")
    return number


def load_dicts(path, columns):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            fail(f"Missing TSV header in {path}")
        if len(set(reader.fieldnames)) != len(reader.fieldnames):
            fail(f"Duplicate column names in {path}")
        missing = set(columns) - set(reader.fieldnames)
        if missing:
            fail(f"Missing columns {sorted(missing)} in {path}")
        result = []
        for line_no, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                fail(f"Malformed TSV row in {path}:{line_no}")
            result.append(row)
        return result


def boolean(value, where):
    value = str(value).strip().lower()
    if value in ("true", "1", "yes"):
        return True
    if value in ("false", "0", "no"):
        return False
    fail(f"Unrecognized Boolean value {value!r} in {where}")


def open_matrix_header(path, n_metadata_cols):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        header = next(csv.reader(handle, delimiter="\t"), None)
    if not header or len(header) <= n_metadata_cols:
        fail(f"Matrix {path} needs >=1 sample after {n_metadata_cols} annotation columns")
    if any(not name for name in header) or len(set(header)) != len(header):
        fail(f"Duplicate/blank column names in {path}")
    if "taxonomy_id" not in header[:n_metadata_cols]:
        fail(f"Matrix {path} requires 'taxonomy_id' among first {n_metadata_cols} columns")
    return header


def matrix_rows(path, expected_cols):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader)
        for lineno, row in enumerate(reader, 2):
            if not row or not any(row):
                continue
            if len(row) != expected_cols:
                fail(f"{path}:{lineno}: expected {expected_cols} columns, found {len(row)}")
            yield row


def count_stats(row, sample_indexes, path):
    total = 0.0
    prevalence = 0
    for i in sample_indexes:
        try:
            x = float(row[i])
        except ValueError:
            fail(f"Non-numeric count {row[i]!r} in {path}")
        if not math.isfinite(x) or x < 0:
            fail(f"Invalid count {row[i]!r} in {path}")
        total += x
        prevalence += x > 0
    return prevalence, total


def run(*, original, filtered, targets, resolution, restored, report,
        n_metadata_cols=3):
    n_metadata_cols = int(n_metadata_cols)
    if n_metadata_cols < 1:
        fail("n_metadata_cols must be >= 1")

    original_header = open_matrix_header(original, n_metadata_cols)
    filtered_header = open_matrix_header(filtered, n_metadata_cols)
    annotation_header = original_header[:n_metadata_cols]
    original_samples = original_header[n_metadata_cols:]
    filtered_samples = filtered_header[n_metadata_cols:]
    if set(original_samples) != set(filtered_samples):
        fail("Filtered count matrix does not have the same sample IDs as original. "
             "Verify apply_qc_exclusions and feature_filter inputs.")
    if set(annotation_header) != set(filtered_header[:n_metadata_cols]):
        fail("Filtered and original matrices have different annotation columns")

    original_taxid_idx = original_header.index("taxonomy_id")
    filtered_taxid_idx = filtered_header.index("taxonomy_id")
    filtered_sample_idx = [filtered_header.index(s) for s in original_samples]
    original_sample_idx = list(range(n_metadata_cols, len(original_header)))

    resolution_rows = load_dicts(resolution, ("target", "category", "rescue_enabled"))
    enabled = {}
    for row in resolution_rows:
        key = row["target"]
        if not key or key in enabled:
            fail(f"Empty/duplicate target label in {resolution}: {key!r}")
        enabled[key] = (boolean(row["rescue_enabled"], f"target {key}"), row["category"])

    # The match table can include taxa with preserve=False. Only rescue those
    # whose label is explicitly rescue_enabled in the resolution report.
    matched_rows = load_dicts(targets, ("target", "taxonomy_id", "matched_taxon"))
    protected = defaultdict(lambda: {"targets": set(), "categories": set(), "names": set()})
    for row in matched_rows:
        label = row["target"]
        if label not in enabled:
            fail(f"Match target {label!r} absent from taxonomy resolution report")
        if not enabled[label][0]:
            continue
        taxid = as_taxid(row["taxonomy_id"], targets)
        protected[taxid]["targets"].add(label)
        protected[taxid]["categories"].add(enabled[label][1])
        if row["matched_taxon"]:
            protected[taxid]["names"].add(row["matched_taxon"])

    # First pass: record feature-filter survivors. We never use their numeric
    # values in the final matrix; original counts are the authoritative source.
    filtered_taxids = set()
    filtered_metadata = {}
    for row in matrix_rows(filtered, len(filtered_header)):
        taxid = as_taxid(row[filtered_taxid_idx], filtered)
        if taxid in filtered_taxids:
            fail(f"Duplicate taxid {taxid} in feature-filtered matrix")
        filtered_taxids.add(taxid)
        filtered_metadata[taxid] = {
            col: row[filtered_header.index(col)] for col in annotation_header
        }
        # Validate counts; also catch unexpected parsing failures early.
        count_stats(row, filtered_sample_idx, filtered)

    seen_original = set()
    selected_filtered = set()
    selected_targets = set()
    audit_rows = []
    retained = 0
    rescued = 0
    Path(restored).parent.mkdir(parents=True, exist_ok=True)
    Path(report).parent.mkdir(parents=True, exist_ok=True)
    with open(restored, "w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow(original_header)
        for row in matrix_rows(original, len(original_header)):
            taxid = as_taxid(row[original_taxid_idx], original)
            if taxid in seen_original:
                fail(f"Duplicate taxid {taxid} in original matrix; cannot rescue unambiguously")
            seen_original.add(taxid)
            in_filtered = taxid in filtered_taxids
            is_target = taxid in protected
            if not (in_filtered or is_target):
                continue
            if in_filtered:
                selected_filtered.add(taxid)
                for field in annotation_header:
                    if str(row[original_header.index(field)]) != str(filtered_metadata[taxid][field]):
                        fail(f"Taxid {taxid} has inconsistent {field!r} between matrices")
            if is_target:
                selected_targets.add(taxid)
            prevalence, total = count_stats(row, original_sample_idx, original)
            writer.writerow(row)
            retained += 1
            if is_target and not in_filtered:
                rescued += 1
            if is_target:
                details = protected[taxid]
                audit_rows.append([
                    taxid,
                    row[original_header.index("name")] if "name" in annotation_header else "",
                    ";".join(sorted(details["targets"])),
                    ";".join(sorted(details["categories"])),
                    str(in_filtered).lower(),
                    str(not in_filtered).lower(),
                    prevalence,
                    format(total, ".15g"),
                    str(total == 0).lower(),
                ])

    missing_filtered = filtered_taxids - selected_filtered
    if missing_filtered:
        fail("Filtered matrix includes taxids missing from original: "
             + ", ".join(str(x) for x in sorted(missing_filtered)[:15]))
    missing_targets = set(protected) - selected_targets
    if missing_targets:
        fail("Taxonomy target matches missing from selected-sample count matrix: "
             + ", ".join(str(x) for x in sorted(missing_targets)[:15]))

    with open(report, "w", newline="", encoding="utf-8") as output:
        writer = csv.writer(output, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "taxonomy_id", "name", "targets", "categories",
            "passed_feature_filter", "rescued_by_target",
            "detected_samples", "total_estimated_reads", "zero_after_exclusions"
        ])
        writer.writerows(audit_rows)

    print(f"rescue_target_taxa: original taxa={len(seen_original)}, "
          f"feature-filter survivors={len(filtered_taxids)}, "
          f"protected observed taxa={len(selected_targets)}, "
          f"newly rescued={rescued}, final taxa={retained}")
    if not selected_targets:
        print("WARNING: no preservation-enabled target taxa matched the selected count matrix")
    if any(row[-1] == "true" for row in audit_rows):
        print("WARNING: some retained target taxa have zero counts in all selected samples; "
              "see zero_after_exclusions in the rescue report")


def run_snakemake(smk):
    run(
        original=str(smk.input["original"]),
        filtered=str(smk.input["filtered"]),
        targets=str(smk.input["targets"]),
        resolution=str(smk.input["resolution"]),
        restored=str(smk.output["restored"]),
        report=str(smk.output["report"]),
        n_metadata_cols=snakemake_param(smk, "n_metadata_cols", 3),
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for name in ("original", "filtered", "targets", "resolution", "restored", "report"):
        ap.add_argument("--" + name, required=True)
    ap.add_argument("--n-metadata-cols", type=int, default=3)
    run(**vars(ap.parse_args()))


if __name__ == "__main__":
    if "snakemake" in globals():
        run_snakemake(snakemake)
    else:
        main()
