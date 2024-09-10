# snakeRnaSeq
(not) just another RNA-seq pipeline

fastq are located in /home/nobackup/molinerislab/RepeatCancer/PignochinoSarcomi/dataset/v1/fastq
we used a symlink in order to save backup space 

pipeline was cloned from [git@github.com:molinerisLab/snakeRnaSeq.git](https://github.com/molinerisLab/snakeRnaSeq)

the conda environment must be activated manually. 
it is located in: /home/molinerislab/RepeatCancer/PignochinoSarcomi/local/env/conda

rule ribo_read_counts.smk expects fastq files to end in: {sample}.fastq.gz

pipeline: 
- run the snakemake rule ALL_infer_experiment
- file ALL_infer_experiment has information on the strandness: columns are as it follows: 
   forward | reverse | unstranded 
- based on which one is the most abundant, strandness is determined. 
- change  the strandness setting in the config.yaml
- you are NOW able to run the all rule with snakemake