# 🐍 snakeRnaSeq
(not) just another RNA-seq pipeline


## 🛠️ Prepare your working directory:

a) git clone the following repository:
```bash
git clone -b develop --single-branch [https://github.com/molinerisLab/snakeRNASeq.git](https://github.com/molinerisLab/snakeRNASeq.git)
```
b) 📜 Clone the SnakeReferences repository: 
```bash
git clone -b develop https://github.com/molinerisLab/SnakeReferences.git
```
---
## ⚙️ Singularity and Snakemake configuration
Before running the pipeline, you must configure which host directories will be mounted inside the Singularity containers. Update the following profile configuration files accordingly:

### 📦 Singularity (default profile)

Edit `profiles/default/config.yaml` and update the `singularity-args` field by listing all directories that need to be bound inside the Singularity container, for example:

```yaml
singularity-args: "-B /path/to/reference/genome -B /path/to/fastq"
```
---
## 🧬 Raw data

Fastq data must be linked in `<work_dir>/fastq/` and must be named like `<sample_name>_R[1|2].fastq.gz`

You can use the `rename` command to normalize the file names.

Create the tab-separated file `metadata.txt` with sample specifications.
The first column must be named `sample` and contain the `<sample_name>` of the fastq files.

---
## 🌐 Environment 
Inside the Singularity container is present a conda environment in which all the necessary tools are present, but since the pipeline is designed to run locally without the container, users can create the corresponding env with 
`conda env create --prefix=local/env/conda --file=local/env/env.yaml`

The environment is directly activated running `direnv allow` when prompted.

You can **avoid** the `--use-conda` option when launching snakemake.

---

## 🚀 Running Snakemake with Singularity
The pipeline can be executed locally, by running snakemake, while the individual rules are executed inside a Singularity container. 

***How it Works*** 
- ***Snakemake*** needs to be installed locally
- **Singularity container:** Defined in the Snakefile, it encapsulates the software dependencies required for running individual workflow rules.

***Usage*** 
To run Snakemake rules, simply run: 

```bash
snakemake examplerule
```
---

## 📊 DGE analysis
For any DGE analysis, please set the dge_tool in the config.yaml. Then you should run snakemake with the --use-conda flag and ask for the file to be made in the DGE folder
    e.g. snakemake -j 1 --use-conda DGE/edger.toptable_clean.ALL_contrast.gz
