"""Apply a reviewed sample-exclusion manifest to a Bracken count matrix and metadata.

Snakemake input:  abundances, metadata, exclusions
Snakemake output: abundances, metadata, audit
Optional params: sample_column, condition_column, g1, g2, min_per_group,
                 n_metadata_cols

The exclusion file is TSV with columns: sample, reason. A header-only file
means no exclusions. No sample/group is excluded automatically.

The output count matrix has the same taxon rows as the input, in the same
order. Only sample columns change. If headers have a Bracken count suffix
(e.g. BD01.bracken_num), these are normalized to metadata IDs (BD01).

Run through Snakemake's script: directive or call --help for standalone usage.
"""

import argparse
import csv
import math
from pathlib import Path


SAMPLE_COLUMNS = ("sample", "sample_id", "sampleid", "sample_name", "id")
BRACKEN_SUFFIXES = (".bracken_num", "_num")


def fail(message):
    raise ValueError(message)


def snakemake_param(smk, name, default=None):
    """Access a named Snakemake parameter (Namedlist has no dict.get)."""
    try:
        return smk.params[name]
    except (AttributeError, IndexError, KeyError, TypeError):
        return default


def read_table(path, delimiter=None):
    """Read small metadata / manifest tables; accept TSV, CSV, or whitespace."""
    with open(path, encoding="utf-8-sig", newline="") as handle:
        first = handle.readline()
        if not first.strip():
            fail(f"Empty table: {path}")
        sep = delimiter or ("\t" if "\t" in first else "," if "," in first else None)
        handle.seek(0)
        rows = list(csv.reader(handle, delimiter=sep)) if sep else [
            line.split() for line in handle if line.strip()
        ]
    if not rows or not rows[0]:
        fail(f"Missing header: {path}")
    header = [h.strip() for h in rows[0]]
    if not all(header) or len(set(header)) != len(header):
        fail(f"Blank or duplicate header columns in {path}")
    data = []
    for line_no, row in enumerate(rows[1:], 2):
        if not row or not any(str(cell).strip() for cell in row):
            continue
        if len(row) != len(header):
            fail(f"{path}:{line_no}: expected {len(header)} columns, found {len(row)}")
        data.append(dict(zip(header, (value.strip() for value in row))))
    return header, data


def find_column(header, requested, aliases, required=True):
    names = {name.lower(): name for name in header}
    if requested:
        match = names.get(str(requested).lower())
        if not match:
            fail(f"Requested column {requested!r} not found. Available: {header}")
        return match
    for name in aliases:
        if name.lower() in names:
            return names[name.lower()]
    if required:
        fail(f"Could not locate column among {aliases}; available columns: {header}")
    return None


def canonical_sample(column, metadata_ids):
    """Only strip a known Bracken suffix if it produces a known metadata ID."""
    if column in metadata_ids:
        return column
    for suffix in BRACKEN_SUFFIXES:
        if column.endswith(suffix) and column[:-len(suffix)] in metadata_ids:
            return column[:-len(suffix)]
    fail(f"Matrix sample column {column!r} has no matching metadata sample ID")


def make_parent(path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def write_tsv(path, header, rows):
    make_parent(path)
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def run(*, abundances, metadata, exclusions, out_abundances, out_metadata,
        audit, sample_column=None, condition_column="condition", g1="Resistant",
        g2="Control", min_per_group=2, n_metadata_cols=3):
    n_metadata_cols = int(n_metadata_cols)
    min_per_group = int(min_per_group)
    if n_metadata_cols < 1 or min_per_group < 1:
        fail("n_metadata_cols and min_per_group must be positive")
    if str(g1) == str(g2):
        fail("g1 and g2 must identify distinct comparison conditions")

    meta_header, meta_rows = read_table(metadata)
    sample_col = find_column(meta_header, sample_column, SAMPLE_COLUMNS, required=False)
    if sample_col is None:
        # Select the sample column by coverage after reading the abundance header.
        candidates = [c for c in meta_header if len(set(row[c] for row in meta_rows)) == len(meta_rows)]
    else:
        candidates = [sample_col]
    if not meta_rows:
        fail("Metadata has no samples")

    excl_header, excl_rows = read_table(exclusions, delimiter="\t")
    exc_sample_col = find_column(excl_header, None, ("sample",))
    exc_reason_col = find_column(excl_header, None, ("reason",))
    excluded = {}
    for row in excl_rows:
        sid, reason = row[exc_sample_col], row[exc_reason_col]
        if not sid or not reason:
            fail("Every exclusion must have a nonempty sample and reason")
        if sid in excluded:
            fail(f"Duplicate exclusion for sample {sid!r}")
        excluded[sid] = reason

    with open(abundances, encoding="utf-8-sig", newline="") as source:
        reader = csv.reader(source, delimiter="\t")
        header = next(reader, None)
    if not header or len(header) <= n_metadata_cols:
        fail(f"Count matrix needs {n_metadata_cols} taxon columns and >=1 sample: {abundances}")
    if any(not h for h in header) or len(set(header)) != len(header):
        fail("Count matrix has blank or duplicate column names")
    taxon_header = header[:n_metadata_cols]
    sample_headers = header[n_metadata_cols:]

    if sample_col is None:
        matches = []
        for candidate in candidates:
            values = set(row[candidate] for row in meta_rows)
            try:
                _ = [canonical_sample(h, values) for h in sample_headers]
                matches.append(candidate)
            except ValueError:
                pass
        if len(matches) != 1:
            fail("Cannot uniquely identify sample ID column in metadata. "
                 f"Candidate columns matching the matrix: {matches}. "
                 "Set params.sample_column explicitly.")
        sample_col = matches[0]

    meta_lookup = {}
    for row in meta_rows:
        sid = row[sample_col]
        if not sid or sid in meta_lookup:
            fail(f"Blank/duplicate metadata sample ID: {sid!r}")
        meta_lookup[sid] = row

    samples = [canonical_sample(h, set(meta_lookup)) for h in sample_headers]
    if len(set(samples)) != len(samples):
        fail("Two matrix columns map to the same sample ID after suffix handling")
    missing_exclusions = set(excluded) - set(samples)
    if missing_exclusions:
        fail("Exclusion manifest references samples absent from the count matrix: "
             + ", ".join(sorted(missing_exclusions)))

    condition_col = find_column(meta_header, condition_column, ("condition",))
    kept_indices = [i for i, sid in enumerate(samples) if sid not in excluded]
    if not kept_indices:
        fail("All samples are excluded")
    kept_ids = [samples[i] for i in kept_indices]
    group_counts = {
        group: sum(meta_lookup[s][condition_col] == group for s in kept_ids)
        for group in (str(g1), str(g2))
    }
    for group, count in group_counts.items():
        if count < min_per_group:
            fail(f"After exclusions, group {group!r} has {count} samples; "
                 f"minimum required is {min_per_group}")

    # Build both outputs in the original sample order, not metadata file order.
    write_tsv(out_metadata, meta_header,
              [[meta_lookup[s][field] for field in meta_header] for s in kept_ids])
    audit_rows = []
    for col, sid in zip(sample_headers, samples):
        audit_rows.append([
            sid, col, meta_lookup[sid][condition_col],
            "excluded" if sid in excluded else "retained", excluded.get(sid, "")
        ])
    write_tsv(audit, ["sample", "input_count_column", "condition", "status", "reason"], audit_rows)

    make_parent(out_abundances)
    taxa_count = 0
    with open(abundances, encoding="utf-8-sig", newline="") as source, \
            open(out_abundances, "w", newline="", encoding="utf-8") as dest:
        reader = csv.reader(source, delimiter="\t")
        writer = csv.writer(dest, delimiter="\t", lineterminator="\n")
        next(reader)
        writer.writerow(taxon_header + kept_ids)
        for line_no, row in enumerate(reader, 2):
            if not row or not any(row):
                continue
            if len(row) != len(header):
                fail(f"{abundances}:{line_no}: expected {len(header)} fields, got {len(row)}")
            for idx in kept_indices:
                value = row[n_metadata_cols + idx]
                try:
                    count = float(value)
                except ValueError:
                    fail(f"Nonnumeric count in row {line_no}, sample {samples[idx]}: {value!r}")
                if not math.isfinite(count) or count < 0:
                    fail(f"Invalid count in row {line_no}, sample {samples[idx]}: {value!r}")
            writer.writerow(row[:n_metadata_cols] + [row[n_metadata_cols + i] for i in kept_indices])
            taxa_count += 1

    print(f"apply_qc_exclusions: {taxa_count} taxa unchanged, "
          f"{len(kept_ids)}/{len(samples)} samples retained, {len(excluded)} excluded")
    print(f"Comparison sample counts: {group_counts}")


def run_snakemake(smk):
    run(
        abundances=str(smk.input["abundances"]),
        metadata=str(smk.input["metadata"]),
        exclusions=str(smk.input["exclusions"]),
        out_abundances=str(smk.output["abundances"]),
        out_metadata=str(smk.output["metadata"]),
        audit=str(smk.output["audit"]),
        sample_column=snakemake_param(smk, "sample_column"),
        condition_column=snakemake_param(smk, "condition_column", "condition"),
        g1=snakemake_param(smk, "g1", "Resistant"),
        g2=snakemake_param(smk, "g2", "Control"),
        min_per_group=snakemake_param(smk, "min_per_group", 2),
        n_metadata_cols=snakemake_param(smk, "n_metadata_cols", 3),
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for required in ("abundances", "metadata", "exclusions", "out-abundances", "out-metadata", "audit"):
        ap.add_argument("--" + required, required=True)
    ap.add_argument("--sample-column")
    ap.add_argument("--condition-column", default="condition")
    ap.add_argument("--g1", default="Resistant")
    ap.add_argument("--g2", default="Control")
    ap.add_argument("--min-per-group", type=int, default=2)
    ap.add_argument("--n-metadata-cols", type=int, default=3)
    args = vars(ap.parse_args())
    run(**{k.replace("-", "_"): v for k, v in args.items()})


if __name__ == "__main__":
    if "snakemake" in globals():
        run_snakemake(snakemake)
    else:
        main()
