# snakeRnaSeq
(not) just another RNA-seq pipeline

This GitHub repository for this pipeline is: https://github.com/molinerisLab/snakeRnaSeq

# goal of the pipeline
Currently, this pipeline does not fully support single end sequencing and expects a forward and a reverse file. 
The final output of the pipeline will be a series of GEP.count.{...} files containing the gene expression counts. 
In particular, a file named GEP.count.exp_filter.ltmm.gz will contain the filtered and normalized reads. 

# file location
The Snakefile expects to find the fastq files in /path/to/folder_name/dataset/v1/fastq.  
However, these files are extremely heavy and so it is useful to store them in a section of the server that is not backed up, and that the folder contains only symlinks to the actual files. 

# file names
The rules expect the fastq files to end in {...}R1.fastq.gz for the forward read, and {...}R2.fastq.gz for the reverse read. 
The rule concerning the name of the files is "star_pe_multi" which can be found in the star_alignment.smk file in /path/to/folder_name/local/rules
This folder also contains all the rules referenced in the Snakefile. 

# environment
Currently, the environment needed to run the pipeline must be activated manually using the following path: 
conda activate /path/to/folder_name/local/env/conda

# how to run the pipeline
- The first step in the pipeline is to run the rule "ALL_infer_experiment" from the Snakefile. 
- This will yield a file with the same name. This file has four columns, the first of which lists the file that was analysed, and the other three show the quantity of reads that were either forward, reverse or unstranded, in this order: 
file name |forward | reverse | unstranded 
- The strandness can be determined based on which of the three is the highest. This is needed to ensure the correct setting is on for the final step
- Once the strandness has been determined, in the config.yaml file (in the same directory as the Snakefile) the "STRANDNESS" setting should be changed accordingly. 
NOTE: the order of the numbers corresponding to the different settings does NOT match the one in the ALL_infer_experiment file. 
- It is now possible to run the "all" rule with Snakemake.
