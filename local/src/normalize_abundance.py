import pandas as pd
import numpy as np

df = pd.read_csv(snakemake.input["matrix"], sep="\t")

meta_cols = df.columns[:3]
num_cols = [c for c in df.columns if c.endswith(".bracken_num")]

counts = df[num_cols]

# ----------------------
# 1. Relative abundance (per sample)
# ----------------------
rel = counts.div(counts.sum(axis=0), axis=1)

rel_df = pd.concat([df[meta_cols], rel], axis=1)
rel_df.to_csv(snakemake.output["relative"], sep="\t", index=False)

# ----------------------
# 2. CLR transformation (per sample)
# ----------------------
pseudocount = 1
counts_pc = counts + pseudocount

geom_mean = np.exp(np.log(counts_pc).mean(axis=0))
clr = np.log(counts_pc.div(geom_mean, axis=1))

clr_df = pd.concat([df[meta_cols], clr], axis=1)
clr_df.to_csv(snakemake.output["clr"], sep="\t", index=False)