#!/usr/bin/env python3
import argparse
import csv
import glob
import os
import re
import sys

def normalize_sample(s: str) -> str:
    s = os.path.basename(s.strip())
    s = re.sub(r"\.k2report$", "", s)
    s = re.sub(r"_unmapped_R[12]\.fastq(\.gz)?$", "", s)
    s = re.sub(r"_R[12]\.fastq(\.gz)?$", "", s)
    return s

def read_depth_csv(path: str, value_col: str):
    data = {}
    with open(path, "r", encoding="utf-8", newline="") as f:
        reader = csv.DictReader(f)
        if "name" not in (reader.fieldnames or []):
            raise ValueError(f"{path} must contain 'name' column")
        if value_col not in (reader.fieldnames or []):
            raise ValueError(f"{path} must contain '{value_col}' column")

        rows = list(reader)
        for row in rows:
            sample = normalize_sample(row["name"])
            data[sample] = row[value_col]
    return rows, data

def read_kreport_counts(kreport_dir: str):
    counts = {}
    for path in glob.glob(os.path.join(kreport_dir, "*.k2report")):
        classified = 0
        unclassified = 0
        human_9606 = 0

        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue

                cols = line.split("\t")
                if len(cols) < 6:
                    cols = line.split()

                if len(cols) >= 8:
                    # unfiltered k2report with minimizer columns
                    # pct, clade, taxon, minimizers, distinct_min, rank, taxid, name
                    rank = cols[5]
                    taxid = cols[6]
                elif len(cols) >= 6:
                    # filtered k2report
                    # pct, clade, taxon, rank, taxid, name
                    rank = cols[3]
                    taxid = cols[4]
                else:
                    continue

                clade_reads = int(float(cols[1]))

                if rank == "R" and taxid == "1":
                    classified = clade_reads
                elif rank == "U" and taxid == "0":
                    unclassified = clade_reads
                elif taxid == "9606":
                    human_9606 = clade_reads

        sample = normalize_sample(path)
        counts[sample] = {
            "kraken_input_reads": classified + unclassified,
            "kraken_unclassified_reads": unclassified,
            "kraken_classified_reads": classified,
            "kraken_human_9606_reads": human_9606,
        }
    return counts

def main():
    parser = argparse.ArgumentParser(
        description="Merge sequencing_depth + unmapped_depth and append kraken_classified_reads from k2report."
    )
    parser.add_argument("--sequencing-csv", required=True)
    parser.add_argument("--unmapped-csv", required=True)
    parser.add_argument("--kreport-dir", required=True)
    parser.add_argument("--output-csv", required=True)
    args = parser.parse_args()

    sequencing_rows, _ = read_depth_csv(args.sequencing_csv, "sequencing_depth")
    _, unmapped_map = read_depth_csv(args.unmapped_csv, "unmapped_depth")
    kraken_map = read_kreport_counts(args.kreport_dir)

    fieldnames = [
        "name",
        "sequencing_depth",
        "unmapped_depth",
        "kraken_input_reads",
        "kraken_unclassified_reads",
        "kraken_classified_reads",
        "kraken_human_9606_reads",
    ]
    out_rows = []

    missing_unmapped = []
    missing_kraken = []

    for row in sequencing_rows:
        sample = normalize_sample(row["name"])
        unmapped = unmapped_map.get(sample, "")
        kraken = kraken_map.get(sample)
        if kraken is None:
            missing_kraken.append(sample)
            kraken = {
                "kraken_input_reads": "",
                "kraken_unclassified_reads": "",
                "kraken_classified_reads": "",
                "kraken_human_9606_reads": "",
            }

        if unmapped == "":
            missing_unmapped.append(sample)
        if kraken == "":
            missing_kraken.append(sample)

        out_rows.append({
            "name": row["name"],
            "sequencing_depth": row["sequencing_depth"],
            "unmapped_depth": unmapped,
            "kraken_input_reads": kraken["kraken_input_reads"],
            "kraken_unclassified_reads": kraken["kraken_unclassified_reads"],
            "kraken_classified_reads": kraken["kraken_classified_reads"],
            "kraken_human_9606_reads": kraken["kraken_human_9606_reads"],
        })

    with open(args.output_csv, "w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(out_rows)

    print(f"Wrote: {args.output_csv}")
    print(f"Rows: {len(out_rows)}")
    print(f"Matched unmapped_depth: {len(out_rows) - len(missing_unmapped)}/{len(out_rows)}")
    print(f"Matched kraken_classified_reads: {len(out_rows) - len(missing_kraken)}/{len(out_rows)}")

    if missing_unmapped:
        print("Missing unmapped_depth for:")
        for s in sorted(set(missing_unmapped)):
            print(f"  {s}")

    if missing_kraken:
        print("Missing k2report for:")
        for s in sorted(set(missing_kraken)):
            print(f"  {s}")

if __name__ == "__main__":
    main()

