"""Split merged Bracken matrix (post-target filtering) into count/fraction TSVs.

Works with either standard Bracken combine *_num/*_frac columns or the
*.bracken_num/*.bracken_frac variant. Run via Snakemake script: directive.
"""
import pandas as pd

input_tsv = str(snakemake.input[0])
result_num = str(snakemake.output['num'])
result_frac = str(snakemake.output['frac'])
df = pd.read_csv(input_tsv, sep='\t')
metadata = ['name', 'taxonomy_id', 'taxonomy_lvl']
if any(col not in df.columns for col in metadata):
    raise ValueError(f'Expected columns {metadata} in {input_tsv}; observed {df.columns.tolist()}')

patterns = [('.bracken_num', '.bracken_frac'), ('_num', '_frac')]
# Match the more specific dotted suffix first (it also ends in "_num").
if any(c.endswith('.bracken_num') for c in df.columns):
    num_suffix, frac_suffix = patterns[0]
elif any(c.endswith('_num') for c in df.columns):
    num_suffix, frac_suffix = patterns[1]
else:
    raise ValueError('No recognized Bracken sample count columns')
num_cols = [c for c in df.columns if c.endswith(num_suffix)]
frac_cols = [c for c in df.columns if c.endswith(frac_suffix)]
num_samples = [c[:-len(num_suffix)] for c in num_cols]
frac_samples = [c[:-len(frac_suffix)] for c in frac_cols]
if len(num_samples) == 0 or set(num_samples) != set(frac_samples):
    raise ValueError('Inconsistent count/fraction columns after Bracken merging')

# Same three metadata columns as the original grep_columns step.
df.loc[:, metadata + num_cols].rename(columns={f'{s}{num_suffix}': s for s in num_samples}).to_csv(
    result_num, sep='\t', index=False)
# Reorder fraction columns to match count columns.
frac_cols = [f'{s}{frac_suffix}' for s in num_samples]
df.loc[:, metadata + frac_cols].rename(columns={f'{s}{frac_suffix}': s for s in num_samples}).to_csv(
    result_frac, sep='\t', index=False)
