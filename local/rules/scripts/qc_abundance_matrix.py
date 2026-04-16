import pandas as pd
import numpy as np
import re

matrix_path = snakemake.input["matrix"]
sample_totals_path = snakemake.input["sample_totals"]

checked_num_matrix_path = snakemake.output["checked_num_matrix"]
checked_frac_matrix_path = snakemake.output["checked_frac_matrix"]
sample_name_map_path = snakemake.output["sample_name_map"]
missing_summary_path = snakemake.output["missing_summary"]
sample_summary_path = snakemake.output["sample_summary"]
taxa_summary_path = snakemake.output["taxa_summary"]
notes_path = snakemake.output["notes"]

min_sample_presence = int(snakemake.params["min_sample_presence"])
min_total_abundance = int(snakemake.params["min_total_abundance"])


# -----------------------------
# Read inputs
# -----------------------------
df = pd.read_csv(matrix_path, sep="\t", dtype=str)
sample_totals = pd.read_csv(sample_totals_path, sep=r"\s+", engine="python")

sample_totals.columns = ["Sample", "Total_Reads"]
sample_totals["Sample"] = sample_totals["Sample"].astype(str)
sample_totals["Total_Reads"] = pd.to_numeric(sample_totals["Total_Reads"], errors="coerce")

expected_samples = set(sample_totals["Sample"])

required_meta = ["name", "taxonomy_id", "taxonomy_lvl"]
missing_meta = [c for c in required_meta if c not in df.columns]
if missing_meta:
    raise ValueError(f"Missing required metadata columns: {missing_meta}")


# -----------------------------
# Identify sample columns
# -----------------------------
num_cols = [c for c in df.columns if c.endswith(".bracken_num")]
frac_cols = [c for c in df.columns if c.endswith(".bracken_frac")]

if len(num_cols) == 0:
    raise ValueError("No '.bracken_num' columns found in abundance matrix.")
if len(frac_cols) == 0:
    raise ValueError("No '.bracken_frac' columns found in abundance matrix.")

def strip_suffix(colname, suffix):
    return re.sub(re.escape(suffix) + r"$", "", colname)

num_samples = [strip_suffix(c, ".bracken_num") for c in num_cols]
frac_samples = [strip_suffix(c, ".bracken_frac") for c in frac_cols]

# map sample names found in matrix
sample_name_map = pd.DataFrame({
    "num_column": num_cols,
    "sample_from_num": num_samples
})

frac_map = pd.DataFrame({
    "frac_column": frac_cols,
    "sample_from_frac": frac_samples
})

sample_name_map = sample_name_map.merge(
    frac_map,
    left_on="sample_from_num",
    right_on="sample_from_frac",
    how="outer"
)

sample_name_map["has_num_column"] = sample_name_map["num_column"].notna()
sample_name_map["has_frac_column"] = sample_name_map["frac_column"].notna()

sample_name_map["sample_name"] = sample_name_map["sample_from_num"].combine_first(
    sample_name_map["sample_from_frac"]
)

sample_name_map["in_sample_totals"] = sample_name_map["sample_name"].isin(expected_samples)

samples_only_in_matrix = sorted(set(sample_name_map["sample_name"]) - expected_samples)
samples_only_in_totals = sorted(expected_samples - set(sample_name_map["sample_name"]))

duplicate_num_samples = pd.Series(num_samples).duplicated().sum()
duplicate_frac_samples = pd.Series(frac_samples).duplicated().sum()

# Keep only samples that have both num and frac columns for paired QC
paired_samples = sorted(set(num_samples).intersection(frac_samples))
paired_num_cols = [f"{s}.bracken_num" for s in paired_samples]
paired_frac_cols = [f"{s}.bracken_frac" for s in paired_samples]

if len(paired_samples) == 0:
    raise ValueError("No samples have both '.bracken_num' and '.bracken_frac' columns.")


# -----------------------------
# Build checked matrices
# -----------------------------
meta_df = df[required_meta].copy()

num_df_raw = df[paired_num_cols].copy()
frac_df_raw = df[paired_frac_cols].copy()

num_df = num_df_raw.apply(pd.to_numeric, errors="coerce")
frac_df = frac_df_raw.apply(pd.to_numeric, errors="coerce")

# Missingness before fill
num_missing_mask = num_df.isna()
frac_missing_mask = frac_df.isna()

num_missing_total = int(num_missing_mask.sum().sum())
frac_missing_total = int(frac_missing_mask.sum().sum())

# Fill missing with 0 for checked matrices
num_df_checked = num_df.fillna(0)
frac_df_checked = frac_df.fillna(0)

checked_num_matrix = pd.concat([meta_df, num_df_checked], axis=1)
checked_frac_matrix = pd.concat([meta_df, frac_df_checked], axis=1)

# Duplicate taxa checks
duplicate_taxon_names = pd.Series(df["name"].astype(str)).duplicated().sum()
duplicate_taxonomy_ids = pd.Series(df["taxonomy_id"].astype(str)).duplicated().sum()


# -----------------------------
# Sample-level QC
# -----------------------------
sample_qc_rows = []
for sample in paired_samples:
    num_col = f"{sample}.bracken_num"
    frac_col = f"{sample}.bracken_frac"

    matrix_total = num_df_checked[num_col].sum()
    detected_taxa = int((num_df_checked[num_col] > 0).sum())
    n_missing_num = int(num_missing_mask[num_col].sum())
    n_missing_frac = int(frac_missing_mask[frac_col].sum())
    frac_sum = frac_df_checked[frac_col].sum()

    bracken_total_row = sample_totals.loc[sample_totals["Sample"] == sample, "Total_Reads"]
    bracken_total = bracken_total_row.iloc[0] if len(bracken_total_row) > 0 else np.nan

    sample_qc_rows.append({
        "Sample": sample,
        "Matrix_Total_Num": matrix_total,
        "Bracken_Total": bracken_total,
        "Total_Difference": matrix_total - bracken_total if pd.notna(bracken_total) else np.nan,
        "Detected_Taxa": detected_taxa,
        "Missing_Num_Values": n_missing_num,
        "Missing_Frac_Values": n_missing_frac,
        "Frac_Sum": frac_sum,
        "Frac_Sum_Close_to_1": np.isclose(frac_sum, 1.0, atol=0.01),
        "Name_Match_in_sample_totals": sample in expected_samples
    })

sample_qc = pd.DataFrame(sample_qc_rows)

# Outlier flag from count totals
valid_totals = sample_qc["Matrix_Total_Num"].dropna()
if len(valid_totals) >= 4:
    q1 = valid_totals.quantile(0.25)
    q3 = valid_totals.quantile(0.75)
    iqr = q3 - q1
    lower = q1 - 1.5 * iqr
    upper = q3 + 1.5 * iqr
    sample_qc["Total_Outlier_Flag"] = (
        (sample_qc["Matrix_Total_Num"] < lower) |
        (sample_qc["Matrix_Total_Num"] > upper)
    )
else:
    sample_qc["Total_Outlier_Flag"] = False


# -----------------------------
# Taxon-level QC
# Based on *.bracken_num columns
# -----------------------------
taxa_prevalence = (num_df_checked > 0).sum(axis=1)
taxa_total_abundance = num_df_checked.sum(axis=1)
taxa_missing_num = num_missing_mask.sum(axis=1)
taxa_missing_frac = frac_missing_mask.sum(axis=1)

taxa_qc = meta_df.copy()
taxa_qc["Samples_Present"] = taxa_prevalence
taxa_qc["Total_Abundance"] = taxa_total_abundance
taxa_qc["Missing_Num_Values"] = taxa_missing_num
taxa_qc["Missing_Frac_Values"] = taxa_missing_frac
taxa_qc["All_Zero_Flag"] = taxa_qc["Total_Abundance"] == 0
taxa_qc["Low_Prevalence_Flag"] = taxa_qc["Samples_Present"] < min_sample_presence
taxa_qc["Low_Abundance_Flag"] = taxa_qc["Total_Abundance"] < min_total_abundance


# -----------------------------
# Missing-value summary
# -----------------------------
missing_summary = pd.DataFrame([
    ["n_rows_taxa", df.shape[0]],
    ["n_metadata_columns", len(required_meta)],
    ["n_num_columns", len(num_cols)],
    ["n_frac_columns", len(frac_cols)],
    ["n_paired_samples", len(paired_samples)],
    ["num_missing_total", num_missing_total],
    ["frac_missing_total", frac_missing_total],
    ["samples_with_num_missing", int((num_missing_mask.sum(axis=0) > 0).sum())],
    ["samples_with_frac_missing", int((frac_missing_mask.sum(axis=0) > 0).sum())],
    ["taxa_with_num_missing", int((num_missing_mask.sum(axis=1) > 0).sum())],
    ["taxa_with_frac_missing", int((frac_missing_mask.sum(axis=1) > 0).sum())],
    ["duplicate_num_sample_names", int(duplicate_num_samples)],
    ["duplicate_frac_sample_names", int(duplicate_frac_samples)],
    ["duplicate_taxon_name_rows", int(duplicate_taxon_names)],
    ["duplicate_taxonomy_id_rows", int(duplicate_taxonomy_ids)],
    ["samples_only_in_matrix", ",".join(samples_only_in_matrix) if samples_only_in_matrix else "None"],
    ["samples_only_in_sample_totals", ",".join(samples_only_in_totals) if samples_only_in_totals else "None"]
], columns=["metric", "value"])


# -----------------------------
# Write outputs
# -----------------------------
checked_num_matrix.to_csv(checked_num_matrix_path, sep="\t", index=False)
checked_frac_matrix.to_csv(checked_frac_matrix_path, sep="\t", index=False)
sample_name_map.to_csv(sample_name_map_path, sep="\t", index=False)
missing_summary.to_csv(missing_summary_path, sep="\t", index=False)
sample_qc.to_csv(sample_summary_path, sep="\t", index=False)
taxa_qc.to_csv(taxa_summary_path, sep="\t", index=False)


# -----------------------------
# Filtering notes
# -----------------------------
n_all_zero = int(taxa_qc["All_Zero_Flag"].sum())
n_low_prev = int(taxa_qc["Low_Prevalence_Flag"].sum())
n_low_abund = int(taxa_qc["Low_Abundance_Flag"].sum())
n_total_outliers = int(sample_qc["Total_Outlier_Flag"].sum())
n_frac_not_close = int((~sample_qc["Frac_Sum_Close_to_1"]).sum())

with open(notes_path, "w") as fh:
    fh.write("Abundance matrix QC notes\n")
    fh.write("=========================\n\n")

    fh.write("Matrix structure\n")
    fh.write("----------------\n")
    fh.write(f"Rows (taxa): {df.shape[0]}\n")
    fh.write(f"Metadata columns: {len(required_meta)}\n")
    fh.write(f"Sample '.bracken_num' columns: {len(num_cols)}\n")
    fh.write(f"Sample '.bracken_frac' columns: {len(frac_cols)}\n")
    fh.write(f"Paired samples with both num and frac columns: {len(paired_samples)}\n\n")

    fh.write("Sample name verification\n")
    fh.write("------------------------\n")
    fh.write(f"Samples present only in matrix: {', '.join(samples_only_in_matrix) if samples_only_in_matrix else 'None'}\n")
    fh.write(f"Samples present only in sample_totals.txt: {', '.join(samples_only_in_totals) if samples_only_in_totals else 'None'}\n")
    fh.write(f"Duplicate num sample names: {duplicate_num_samples}\n")
    fh.write(f"Duplicate frac sample names: {duplicate_frac_samples}\n\n")

    fh.write("Missing values\n")
    fh.write("--------------\n")
    fh.write(f"Total missing values in num block: {num_missing_total}\n")
    fh.write(f"Total missing values in frac block: {frac_missing_total}\n\n")

    fh.write("Initial filtering notes\n")
    fh.write("-----------------------\n")
    fh.write(f"- {n_all_zero} taxa have zero counts across all samples.\n")
    fh.write(f"- {n_low_prev} taxa are present in fewer than {min_sample_presence} samples.\n")
    fh.write(f"- {n_low_abund} taxa have total abundance < {min_total_abundance}.\n")
    fh.write(f"- {n_total_outliers} samples are flagged as total-abundance outliers.\n")
    fh.write(f"- {n_frac_not_close} samples have fraction-column sums not close to 1.0.\n\n")

    fh.write("Interpretation\n")
    fh.write("--------------\n")
    fh.write("- Use sample_qc_summary.tsv to compare matrix totals against source Bracken totals.\n")
    fh.write("- Consider removing taxa with All_Zero_Flag == True.\n")
    fh.write("- Consider filtering taxa with Low_Prevalence_Flag == True for downstream analyses.\n")
    fh.write("- Inspect samples with large Total_Difference or Total_Outlier_Flag before exclusion.\n")
    fh.write("- Inspect samples where Frac_Sum_Close_to_1 == False, as this may indicate matrix-construction issues.\n")