"""Taxonomy-aware filtering for a combined Bracken abundance matrix.

Snakemake integration: expects input[matrix,nodes,names], output[filtered,summary,
resolution,matches], params[min_prevalence,min_total,targets,targets_enabled,
preserve_targets] and log[0].

Important: names.dmp and nodes.dmp must be readable NCBI-format files representing
THE SAME TAXONOMY as the Kraken2 database that generated the Bracken data.
"""
from __future__ import annotations

import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd


def normalise(value):
    return " ".join(str(value).split()).casefold()


def load_nodes(path):
    """Return {taxid: parent_taxid} and {taxid: rank} from NCBI nodes.dmp."""
    parents, ranks = {}, {}
    with open(path, encoding="utf-8") as handle:
        for line_num, line in enumerate(handle, 1):
            fields = [item.strip() for item in line.split("|")]
            if len(fields) < 3:
                raise ValueError(f"Malformed nodes.dmp line {line_num}: {line[:80]!r}")
            child, parent = int(fields[0]), int(fields[1])
            parents[child] = parent
            ranks[child] = fields[2]
    if not parents:
        raise ValueError(f"No taxonomy nodes in {path}")
    return parents, ranks


def load_requested_names(path, requested_names):
    """Only index names requested by the YAML targets, not all NCBI names.

    Scientific names take precedence over synonyms; ambiguous matches fail
    rather than selecting an arbitrary taxid.
    """
    synonyms = defaultdict(set)
    primary = defaultdict(set)
    scientific_name = {}
    requested = {normalise(n) for n in requested_names if n}
    with open(path, encoding="utf-8") as handle:
        for line_num, line in enumerate(handle, 1):
            fields = [item.strip() for item in line.split("|")]
            if len(fields) < 4:
                raise ValueError(f"Malformed names.dmp line {line_num}: {line[:80]!r}")
            taxid = int(fields[0])
            name = fields[1]
            kind = fields[3]
            key = normalise(name)
            if key in requested:
                synonyms[key].add(taxid)
                if kind == "scientific name":
                    primary[key].add(taxid)
            if kind == "scientific name" and key in requested:
                scientific_name[taxid] = name
    return primary, synonyms, scientific_name


def resolve_name(name, scientific, aliases):
    key = normalise(name)
    candidates = scientific.get(key, set()) or aliases.get(key, set())
    if not candidates:
        return None, "not found in names.dmp"
    if len(candidates) != 1:
        raise ValueError(f"Ambiguous taxonomy name {name!r}: candidate taxids {sorted(candidates)}")
    return next(iter(candidates)), "scientific name" if scientific.get(key) else "synonym"


def has_ancestor(taxid, ancestor, parents):
    """True for the taxid itself or any descendant; safe for malformed cycles."""
    current = taxid
    seen = set()
    while current in parents and current not in seen:
        if current == ancestor:
            return True
        seen.add(current)
        next_node = parents[current]
        if next_node == current:
            break
        current = next_node
    return False


def row_count_columns(df):
    """Support upstream Bracken combine implementations using either suffix."""
    dotted = [c for c in df.columns if c.endswith(".bracken_num")]
    underscored = [c for c in df.columns if c.endswith("_num")]
    if dotted and underscored:
        raise ValueError("Both '.bracken_num' and '_num' count columns detected; choose one format")
    num_cols = dotted or underscored
    if not num_cols:
        raise ValueError("No Bracken counts found (expected *_num or *.bracken_num columns)")
    return num_cols


def run(snakemake):
    # Snakemake passes output/log values as paths or string-like values.
    log_path = Path(str(snakemake.log[0]))
    log_path.parent.mkdir(parents=True, exist_ok=True)
    with log_path.open("w", encoding="utf-8") as logfile:
        old_out, old_err = sys.stdout, sys.stderr
        sys.stdout = sys.stderr = logfile
        try:
            _execute(snakemake)
        finally:
            sys.stdout, sys.stderr = old_out, old_err


def _execute(smk):
    matrix_path = str(smk.input["matrix"])
    node_path = str(smk.input["nodes"])
    name_path = str(smk.input["names"])
    targets = smk.params["targets"] or {}
    min_prevalence = int(smk.params["min_prevalence"])
    min_total = float(smk.params["min_total"])
    enabled = bool(smk.params["targets_enabled"])
    preserve_targets = bool(smk.params["preserve_targets"])
    if not isinstance(targets, dict):
        raise ValueError("TARGET_TAXA.targets must be a mapping keyed by target label")

    print(f"Matrix: {matrix_path}")
    df = pd.read_csv(matrix_path, sep="\t")
    if "taxonomy_id" not in df.columns:
        raise ValueError("Combined Bracken matrix requires 'taxonomy_id' column")
    count_cols = row_count_columns(df)
    print(f"Input taxa: {len(df)}; samples: {len(count_cols)}")
    counts = df[count_cols].apply(pd.to_numeric, errors="raise")
    if counts.isna().any().any() or (counts < 0).any().any():
        raise ValueError("Bracken abundance values must be nonmissing and nonnegative")
    taxids = pd.to_numeric(df["taxonomy_id"], errors="coerce").astype("Int64")
    prevalence = (counts > 0).sum(axis=1)
    totals = counts.sum(axis=1)
    ordinary = (prevalence >= min_prevalence) & (totals >= min_total)

    # Load just the relevant taxonomy information. When targets are entirely
    # disabled, the taxonomy input paths remain declared but do not need parse.
    parents, ranks, scientific, aliases = {}, {}, {}, {}
    if enabled:
        parents, ranks = load_nodes(node_path)
        names = [t.get("name") for t in targets.values()
                 if isinstance(t, dict) and t.get("enabled", True) and t.get("match", "taxonomy") == "taxonomy"]
        scientific, aliases, _ = load_requested_names(name_path, names)

    is_target = pd.Series(False, index=df.index, dtype=bool)
    eligible_for_rescue = pd.Series(False, index=df.index, dtype=bool)
    resolution = []
    matches = []
    for key, entry in targets.items():
        if not isinstance(entry, dict):
            raise ValueError(f"Target {key!r} should contain a YAML mapping")
        label = str(entry.get("name", key.replace("_", " ")))
        category = str(entry.get("category", "unspecified"))
        active = enabled and bool(entry.get("enabled", True))
        rescues = preserve_targets and bool(entry.get("preserve", True))
        include_children = bool(entry.get("include_children", False))
        mode = entry.get("match", "taxonomy")
        resolved_id = None
        status = "disabled" if not active else "unresolved"
        note = ""
        flag = pd.Series(False, index=df.index, dtype=bool)
        if active:
            if mode == "taxonomy":
                config_taxid = entry.get("taxid")
                requested_id = int(config_taxid) if config_taxid is not None else None
                name_id, name_message = resolve_name(label, scientific, aliases)
                if requested_id is not None and requested_id not in parents:
                    note = f"Configured taxid {requested_id} absent from supplied nodes.dmp"
                elif requested_id is not None and name_id is not None and name_id != requested_id:
                    raise ValueError(f"{key}: name {label!r} resolves to {name_id}, not configured {requested_id}")
                else:
                    resolved_id = requested_id if requested_id is not None else name_id
                    if resolved_id is not None and resolved_id not in parents:
                        note = f"Taxid {resolved_id} missing from supplied nodes.dmp"
                        resolved_id = None
                    if resolved_id is not None:
                        expected_rank = entry.get("rank")
                        actual_rank = ranks.get(resolved_id)
                        if expected_rank and normalise(expected_rank) != normalise(actual_rank):
                            raise ValueError(f"{key}: configured rank {expected_rank!r} differs from taxonomy rank {actual_rank!r}")
                        flag = taxids.map(lambda x: False if pd.isna(x) else (
                            has_ancestor(int(x), resolved_id, parents) if include_children else int(x) == resolved_id
                        )).astype(bool)
                        status = "resolved" if flag.any() else "resolved_no_hits"
                        note = f"Resolved from {name_message}; include_children={include_children}"
                    else:
                        # Exact matrix-name fallback is intentionally limited to the
                        # displayed name. This does NOT imply a resolved lineage.
                        if "name" in df.columns and label:
                            flag = df["name"].fillna("").map(normalise).eq(normalise(label))
                        status = "unresolved_name_fallback" if flag.any() else "unresolved_no_hits"
                        note = (note + "; " if note else "") + (
                            "Only exact matrix-name matching was available; no descendant matching"
                        )
            elif mode == "matrix_name_regex":
                if "name" not in df.columns:
                    raise ValueError(f"{key}: 'name' column needed for matrix_name_regex")
                pattern = entry.get("pattern")
                if not pattern or not isinstance(pattern, str):
                    raise ValueError(f"{key}: regex mode requires a non-empty 'pattern'")
                try:
                    compiled = re.compile(pattern)
                except re.error as exc:
                    raise ValueError(f"{key}: bad regex {pattern!r}: {exc}") from exc
                flag = df["name"].fillna("").astype(str).map(lambda s: compiled.search(s) is not None)
                status = "regex_hits" if flag.any() else "regex_no_hits"
                note = "Unverified matrix-name regex; does not resolve taxonomy or infer ancestry"
            else:
                raise ValueError(f"{key}: match must be 'taxonomy' or 'matrix_name_regex'")

        is_target |= flag
        if rescues:
            eligible_for_rescue |= flag
        match_count = int(flag.sum())
        resolution.append({
            "target": key, "configured_name": label, "category": category,
            "mode": mode, "resolved_taxid": "" if resolved_id is None else resolved_id,
            "status": status, "matching_rows": match_count,
            "total_estimated_reads": float(totals[flag].sum()),
            "rescue_enabled": bool(active and rescues), "note": note,
        })
        for idx in df.index[flag]:
            matches.append({
                "target": key, "category": category,
                "matched_taxon": str(df.at[idx, "name"]) if "name" in df else "",
                "taxonomy_id": df.at[idx, "taxonomy_id"],
                "prevalence": int(prevalence.at[idx]),
                "total_estimated_reads": float(totals.at[idx]),
                "passed_standard": bool(ordinary.at[idx]),
                "rescued_by_target": bool(active and rescues and not ordinary.at[idx]),
                "match_mode": mode,
            })
        print(f"{key}: {status}; taxid={resolved_id}; matched={match_count}; note={note}")

    # Critical distinction: BOTH numerical thresholds, or rescue of a target.
    rescued = enabled and preserve_targets
    rescued_mask = (eligible_for_rescue & ~ordinary) if rescued else pd.Series(False, index=df.index)
    keep = ordinary | (eligible_for_rescue if rescued else False)
    df.loc[keep].to_csv(str(smk.output["filtered"]), sep="\t", index=False)

    fields_resolution = ["target", "configured_name", "category", "mode", "resolved_taxid", "status",
                         "matching_rows", "total_estimated_reads", "rescue_enabled", "note"]
    fields_matches = ["target", "category", "matched_taxon", "taxonomy_id", "prevalence",
                      "total_estimated_reads", "passed_standard", "rescued_by_target", "match_mode"]
    for path, rows, fields in [(smk.output["resolution"], resolution, fields_resolution),
                               (smk.output["matches"], matches, fields_matches)]:
        Path(str(path)).parent.mkdir(parents=True, exist_ok=True)
        with open(str(path), "w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    unresolved = sum(x["status"].startswith("unresolved") for x in resolution)
    text = [
        f"Input taxa: {len(df)}",
        f"Prevalence >= {min_prevalence} AND total estimated reads >= {min_total}",
        f"Passed standard filtering: {int(ordinary.sum())}",
        f"Detected target taxa (unique matrix rows): {int(is_target.sum())}",
        f"Target-rescued taxa (unique matrix rows): {int(rescued_mask.sum())}",
        f"Retained taxa: {int(keep.sum())}",
        f"Removed taxa: {int((~keep).sum())}",
        f"Unresolved target labels: {unresolved}",
        "Unresolved names do not imply absence. Check database sequence coverage.",
        "Target preservation here does not override later feature_filter rules.",
    ]
    Path(str(smk.output["summary"])).write_text("\n".join(text) + "\n", encoding="utf-8")
    for line in text:
        print(line)


if __name__ == "__main__":
    # For use with Snakemake script: directive only.
    run(snakemake)
