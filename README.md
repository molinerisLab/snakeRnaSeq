# snakeRnaSeq
(not) just another RNA-seq pipeline

# DGE analysis
For any DGE analysis, please set the dge_tool in the config.yaml. Then you should run snakemake with the --use-conda flag and ask for the file to be made in the DGE folder
    e.g. snakemake -j 1 --use-conda DGE/edger.toptable_clean.ALL_contrast.gz