import os
import re
import csv
import gzip

log_dir = "star/human"
unmapped_dir = "fastq/unmapped/"
k2report_dir = "kreports"  
output_csv = "reads_summary.csv"
contaminants_file = "Reports.csv" # Path to your contaminants file
metadata_file = "metadata.txt"  # Path to your metadata file

# Read contaminant data and aggregate by sample
contaminants_by_sample = {}
if os.path.exists(contaminants_file):
    with open(contaminants_file, 'r') as f:
        reader = csv.reader(f, delimiter='\t')
        headers = next(reader)  # Skip header
        for row in reader:
            if len(row) >= 3:  
                sample = row[0]
                reads = row[2]
                try:
                    reads_count = int(reads)
                    if sample in contaminants_by_sample:
                        contaminants_by_sample[sample] += reads_count
                    else:
                        contaminants_by_sample[sample] = reads_count
                except ValueError:
                    # Skip if reads is not an integer
                    pass

# Read metadata and map sample to sample_type
sample_type_by_sample = {}
if os.path.exists(metadata_file):
    with open(metadata_file, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            sample = row.get('sample')
            sample_type = row.get('sample_type')
            if sample and sample_type:
                sample_type_by_sample[sample] = sample_type

rows = []

for filename in os.listdir(log_dir):
    if filename.endswith(".Log.final.out"):
        sample = filename.replace(".Log.final.out", "")
        unmapped_fastq = os.path.join(unmapped_dir, f"{sample}_unmapped_R1.fastq.gz")
        # Calculate unmapped reads by counting lines in gzipped FASTQ and dividing by 4
        if not os.path.exists(unmapped_fastq):
            unmapped_reads = ""
        else:
            with gzip.open(unmapped_fastq, "rt") as f_unmapped:
                line_count = sum(1 for _ in f_unmapped)
                unmapped_reads = str(line_count // 4)
        # Get classified kraken reads and human assigned reads from .k2report
        k2report_file = os.path.join(k2report_dir, f"{sample}.k2report")
        classified_kraken = ""
        kraken_human_assigned_reads = ""
        if os.path.exists(k2report_file):
            with open(k2report_file) as f_k2:
                lines = f_k2.readlines()
                # classified_kraken: second row, second column
                if len(lines) >= 2:
                    cols = lines[1].strip().split("\t")
                    if len(cols) >= 2:
                        classified_kraken = cols[1]
                # kraken_human_assigned_reads: find row with 5th col == "9606", get 2nd col
                for line in lines:
                    cols = line.strip().split("\t")
                    if len(cols) >= 6 and cols[4 ] == "9606":
                        kraken_human_assigned_reads = cols[1]
                        break
        with open(os.path.join(log_dir, filename)) as f:
            content = f.read()
            # STAR log tells us input reads, but not initial pre-trimming reads.
            # We will use STAR input as "trimmed reads" and leave "initial reads" blank or duplicate it.
            m = re.search(r"Number of input reads \|\t(\d+)", content)
            if m:
                initial_reads = m.group(1)  # Or leave blank "" if you only want trimmed
                trimmed_reads = m.group(1)
                try:
                    mapped_reads = str(int(trimmed_reads) - int(unmapped_reads))
                except ValueError:
                    mapped_reads = ""
                
                # Calculate percentage columns
                percent_trimmed = ""
                percent_mapped = ""
                percent_classified_kraken = ""
                
                try:
                    if initial_reads and trimmed_reads:
                        percent_trimmed = f"{(int(trimmed_reads) / int(initial_reads)) * 100:.2f}"
                except (ValueError, ZeroDivisionError):
                    pass
                
                try:
                    if mapped_reads and trimmed_reads:
                        percent_mapped = f"{(int(mapped_reads) / int(trimmed_reads)) * 100:.2f}"
                except (ValueError, ZeroDivisionError):
                    pass
                
                try:
                    if classified_kraken and unmapped_reads:
                        percent_classified_kraken = f"{(int(classified_kraken) / int(unmapped_reads)) * 100:.2f}"
                except (ValueError, ZeroDivisionError):
                    pass
                
                # Get the contaminant reads for this sample
                contaminant_reads = contaminants_by_sample.get(sample, 0)
                sample_type = sample_type_by_sample.get(sample, "")
                rows.append([
                    sample, sample_type,  # Insert sample_type as second column
                    initial_reads, trimmed_reads, unmapped_reads, mapped_reads,
                    classified_kraken, kraken_human_assigned_reads,
                    percent_trimmed, percent_mapped, percent_classified_kraken,
                    contaminant_reads
                ])

reports_dir = "Reports"
os.makedirs(reports_dir, exist_ok=True)

with open(os.path.join(reports_dir, output_csv), "w", newline="") as csvfile:
    writer = csv.writer(csvfile, delimiter="\t")
    writer.writerow([
        "sample", "sample_type",  # Add sample_type as second header
        "initial_reads", "trimmed_reads", "unmapped_reads",
        "mapped_reads", "classified_kraken", "Kraken_human_assigned_reads",
        "%trimmed", "%mapped", "%classified_kraken", "contaminant_reads"
    ])
    writer.writerows(rows)

print(f"CSV file '{os.path.join(reports_dir, output_csv)}' created with {len(rows)} samples.")

try:
    split_by = snakemake.config.get("split_by")
except NameError:
    split_by = None

if split_by:
    # Read metadata and get all unique values for split_by
    unique_values = set()
    sample_by_value = {}
    if os.path.exists(metadata_file):
        with open(metadata_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                value = row.get(split_by)
                if isinstance(value, list):
                    value = value[0]
                if value:
                    unique_values.add(value)
                    sample_by_value.setdefault(value, set()).add(row['sample'])

    # For each unique value, filter rows and write a CSV
    for value in unique_values:
        filtered_samples = sample_by_value[value]
        filtered_rows = [row for row in rows if row[0] in filtered_samples]
        filtered_csv = os.path.join(reports_dir, f"{value}_reads_summary.csv")
        with open(filtered_csv, "w", newline="") as csvfile:
            writer = csv.writer(csvfile, delimiter="\t")
            writer.writerow([
                "sample", "sample_type",
                "initial_reads", "trimmed_reads", "unmapped_reads",
                "mapped_reads", "classified_kraken", "Kraken_human_assigned_reads",
                "%trimmed", "%mapped", "%classified_kraken", "contaminant_reads"
            ])
            writer.writerows(filtered_rows)
        print(f"Filtered CSV '{filtered_csv}' created with {len(filtered_rows)} samples.")