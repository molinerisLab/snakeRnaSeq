# snakeRnaSeq
(not) just another RNA-seq pipeline

This GitHub repository for this pipeline is: https://github.com/molinerisLab/snakeRnaSeq

The working directory for this pipeline is /home/molinerislab/RepeatCancer/PignochinoSarcomi/dataset/v1

The Snakefile expects to find the Fastq files in /home/molinerislab/RepeatCancer/PignochinoSarcomi/dataset/v1/fastq
However, the fastq files take up a lot of space, and are therefore stored in a section of the system that is not backed up and should only be pointed to with symlinks. 
The folder where the files should be uploaded to is: /home/nobackup/molinerislab/RepeatCancer/PignochinoSarcomi/dataset/v1/fastq

The rules also expect the fastq files to end in {...}R1.fastq.gz for the forward read, and {...}R2.fastq.gz for the reverse read. 
the rule concerning the name of the files is "star_pe_multi" which can be found in the star_alignment.smk file in /home/molinerislab/RepeatCancer/PignochinoSarcomi/local/rules
This folder also contains all the rules referenced in the Snakefile. 

Currently, the environment needed to run the pipeline must be activated manually using the following path: 
conda activate /home/molinerislab/RepeatCancer/PignochinoSarcomi/local/env/conda

- The first step in the pipeline is to run the rule "ALL_infer_experiment" from the Snakefile. 
- This will yield a file with the same name. This file has four columns, the first of which lists the file that was analysed, and the other three show the quantity of reads that were either forward, reverse or unstranded, in this order: 
file name |forward | reverse | unstranded 
- the strandness can be determined based on which of the three is the highest. This is needed to ensure the correct setting is on for the final step
- once the strandness has been determined, in the config.yaml file (in the same directory as the snakefile) the "STRANDNESS" setting should be changed accordingly. 
NOTE: the order of the numbers corrisponding to the different settings does NOT match the one in the ALL_infer_experiment file. 
- it is now possible to run the "all" rule with Snakemake. 


