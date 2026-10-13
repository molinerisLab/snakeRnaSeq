import pandas as pd
import sys
from collections import defaultdict

# Access Snakemake variables
input_matrix = snakemake.input["matrix"]
taxonomy_nodes = snakemake.input["taxonomy_nodes"]

output_filtered = snakemake.output["filtered"]
output_summary = snakemake.output["summary"]

min_prevalence = int(snakemake.params["min_prevalence"])
min_total = float(snakemake.params["min_total"])

targets = snakemake.params["targets"]
targets_enabled = bool(snakemake.params["targets_enabled"])
preserve_targets = bool(snakemake.params["preserve_targets"])

log_file = snakemake.log[0]

# Redirect output to log
sys.stdout = open(log_file, "w")
sys.stderr = sys.stdout

# Helper: read NCBI taxonomy parent-child structure
def load_taxonomy_parents(nodes_file):
    parents = {}

    with open(nodes_file) as handle:
        for line in handle:
            fields = [x.strip() for x in line.split("|")]

            taxid = int(fields[0])
            parent_taxid = int(fields[1])

            parents[taxid] = parent_taxid

    return parents


def is_target_or_descendant(taxid, target_ids, parents):
    """Return True if taxid is a target or descends
    from a target configured with include_children=True.
    """
    if pd.isna(taxid):
        return False

    taxid = int(taxid)

    # Exact target match is handled separately.
    current = taxid
    visited = set()

    while current in parents and current not in visited:
        if current in target_ids:
            return True

        visited.add(current)
        parent = parents[current]

        # Root reached
        if parent == current:
            break

        current = parent

    return False

# Load merged abundance matrix
print("Loading merged Bracken matrix...")
df = pd.read_csv(input_matrix, sep="\t")

# Identify sample columns and taxonomic ID column
# meta_cols = df.columns[:3]
# num_cols = [c for c in df.columns if c.endswith(".bracken_num")]

# Adjust this name if your matrix uses a different taxid column
taxid_col = "taxonomy_id"

if taxid_col not in df.columns:
    raise ValueError(
        f"Taxonomic ID column '{taxid_col}' not found. "
        f"Available columns: {df.columns.tolist()}"
    )

num_cols = [col for col in df.columns if col.endswith(".bracken_num")]

if not num_cols:
    raise ValueError("No sample abundance columns ending in '.bracken_num' found.")

print(f"Detected {len(num_cols)} sample columns")
print(f"Loaded {len(df)} taxa")


# Compute filtering metrics
df[num_cols] = df[num_cols].apply(
    pd.to_numeric, errors="raise"
)

df["prevalence"] = (df[num_cols] > 0).sum(axis=1)
df["total_abundance"] = df[num_cols].sum(axis=1)

# Filtering conditions
passes_prevalence_filter = (df["prevalence"] >= min_prevalence)
passes_abundance_filter = (df["total_abundance"] >= min_total)
passes_standard_filter = (passes_prevalence_filter & passes_abundance_filter)

# Identify configured target taxa
taxid_values = pd.to_numeric(df[taxid_col], errors="coerce")
taxid_values = taxid_values.astype("Int64")

exact_target_ids = set()
descendant_target_ids = set()

if targets_enabled and preserve_targets:

    for target_name, target in targets.items():
        target_id = int(target["taxid"])

        # Always preserve the exact configured taxon.
        exact_target_ids.add(target_id)

        print(
            f"Configured target: {target_name} "
            f"(taxid={target_id}, "
            f"rank={target.get('rank', 'unspecified')}, "
            f"include_children="
            f"{target.get('include_children', False)})"
        )

        # Descendants are included only when requested.
        if target.get("include_children", False):
            descendant_target_ids.add(target_id)

    # Exact target matches
    taxon_is_target = taxid_values.isin(exact_target_ids)

    # Include descendants of selected targets
    if descendant_target_ids:
        parents = load_taxonomy_parents(taxonomy_nodes)

        descendant_mask = taxid_values.map(
            lambda x: is_target_or_descendant(
                x, descendant_target_ids, parents
            )
        )

        taxon_is_target = taxon_is_target | descendant_mask

else:
    taxon_is_target = pd.Series(
        False, index=df.index
    )

# Combine standard filtering and preservation
keep_taxon = (passes_prevalence_filter | taxon_is_target)
filtered = df.loc[keep_taxon].copy()

# Remove temporary columns
filtered = filtered.drop(columns=["prevalence", "total_abundance"])
filtered.to_csv(output_filtered, sep="\t", index=False)


# Summary stats
total_taxa = len(df)
retained_taxa = int(keep_taxon.sum())
removed_taxa = total_taxa - retained_taxa

retained_by_target_rescue = int(
    (taxon_is_target & ~passes_standard_filter).sum()
)
target_taxa_present = int(taxon_is_target.sum())

print("Writing summary...")

with open(output_summary, "w") as handle:
    handle.write(f"Total taxa: {total_taxa}\n")
    handle.write(f"Retained taxa: {retained_taxa}\n")
    handle.write(f"Removed taxa: {removed_taxa}\n\n")

    handle.write("Standard filtering criteria:\n")
    handle.write(f"- Prevalence >= {min_prevalence} samples AND\n")
    handle.write(f"- Total abundance >= {min_total}\n\n")

    handle.write("Target-taxon preservation:\n")
    handle.write(f"- Target system enabled: {targets_enabled}\n")
    handle.write(f"- Preserve targets: {preserve_targets}\n")
    handle.write(
        f"- Target taxa/descendants present: "
        f"{target_taxa_present}\n"
)
    handle.write(
        f"- Taxa rescued from standard filtering: "
        f"{retained_by_target_rescue}\n"
    )

print("Filtering complete.")