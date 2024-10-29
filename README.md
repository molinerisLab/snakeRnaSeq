# snakeRnaSeq
(not) just another RNA-seq pipeline


# Configuration

## Reference data

Clone the GENCODE repo https://bitbucket.org/irccit/gencode

## Raw data 

Fastq data must be linked in `<work_dir>/fastq/` and must be named like `<sample_name>_R[1|2].fastq.gz`

You can use the `rename` command to normalize the file names. 
<!-- TODO link for rename -->

## Metadata

Create the tab-separated file `metadata.txt` with sample specifications.
The first column must be named `sample` and contain the `<sample_name>` of the fastq files.

## Environment 

Most rules work in a common conda environment defined in `local/env/env.yaml`.
You should create the corresponding env with 
`conda create --prefix=local/env/conda --file=local/env/env.yaml`

The environment is directly activated running `direnv allow` when prompted.

You can **avoid** the `--use-conda` option when launching snakemake.


# Usage

To get all the alignments, quality controls and gene counts run:
```
snakemake -p -j 1 all
```


# Pipeline steps

## Alignment with STAR

## Clean

## DGE analysis
For any DGE analysis, please set the dge_tool in the config.yaml. Then you should run snakemake with the --use-conda flag and ask for the file to be made in the DGE folder
    e.g. snakemake -j 1 --use-conda DGE/edger.toptable_clean.ALL_contrast.gz
