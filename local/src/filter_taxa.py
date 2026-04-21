import pandas as pd

# Access Snakemake variables
input_matrix = snakemake.input["matrix"]
output_filtered = snakemake.output["filtered"]
output_summary = snakemake.output["summary"]

min_prevalence = snakemake.params["min_prevalence"]
min_total = snakemake.params["min_total"]

log_file = snakemake.log[0]

# Redirect prints to log file
import sys
sys.stdout = open(log_file, "w")
sys.stderr = sys.stdout

print("Loading matrix...")
df = pd.read_csv(input_matrix, sep="\t")

# Identify columns
meta_cols = df.columns[:3]
num_cols = [c for c in df.columns if c.endswith(".bracken_num")]

print(f"Detected {len(num_cols)} sample columns")

# Compute metrics
df["prevalence"] = (df[num_cols] > 0).sum(axis=1)
df["total_abundance"] = df[num_cols].sum(axis=1)

# Apply filtering
filtered = df[
    (df["prevalence"] >= min_prevalence) &
    (df["total_abundance"] >= min_total)
].copy()

# Save filtered matrix
filtered.drop(columns=["prevalence", "total_abundance"]) \
        .to_csv(output_filtered, sep="\t", index=False)

# Summary stats
total_taxa = len(df)
retained_taxa = len(filtered)
removed_taxa = total_taxa - retained_taxa

print("Writing summary...")
with open(output_summary, "w") as f:
    f.write(f"Total taxa: {total_taxa}\n")
    f.write(f"Retained taxa: {retained_taxa}\n")
    f.write(f"Removed taxa: {removed_taxa}\n\n")
    f.write("Filtering criteria:\n")
    f.write(f"- Prevalence >= {min_prevalence} samples\n")
    f.write(f"- Total abundance >= {min_total}\n")

print("Filtering complete.")