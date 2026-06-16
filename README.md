# 🐍 snakeRnaSeq
*(not) just another RNA-seq pipeline*

A comprehensive, fully-automated workflow powered by Snakemake to process RNA-seq data from raw reads to Differential Gene Expression (DGE), Isoform analysis, and Metagenomic profiling.

### ✨ Features
* **Quality Control**: Automated read filtering and adapter trimming (`fastp`).
* **Alignment & Quantification**: High-performance host alignment (`STAR`) and transcript quantification.
* **Isoform Discovery & Analysis**: Comprehensive profiling, quantification, and differential usage of transcript isoforms.
* **Differential Gene Expression (DGE)**: Automated contrast analysis utilizing standard R packages.
* **Metagenomic Profiling**: Repurposes unmapped reads to discover and quantify the sample's microbiome (bacteria, viruses, fungi) using `Kraken2` and `Bracken`.

### 🖥️ Hardware Prerequisites
* **Environment**: Linux server or HPC cluster with **Snakemake** and **Conda** installed.
* **CPU**: Highly parallelizable; **16–32+ compute cores** recommended for acceptable runtimes.
* **RAM**: 
  * **~32GB+** for standard RNA-seq analysis (STAR host alignment).
  * **~700GB+** if executing the optional *Metagenomics module* with the GTDB database.

---

## 🛠️ Prepare your working directory

a) git clone the following repository:
```bash
git clone -b develop --single-branch [https://github.com/molinerisLab/snakeRNASeq.git](https://github.com/molinerisLab/snakeRNASeq.git)
```
b) 📜 Clone the SnakeReferences repository: 
```bash
git clone -b develop https://github.com/molinerisLab/SnakeReferences.git
```

---

## 🧬 Raw data

Fastq data must be linked in `<work_dir>/fastq/` and must be named like `<sample_name>_R[1|2].fastq.gz`

You can use the `rename` command to normalize the file names.

Create the tab-separated file `metadata.txt` with sample specifications.
The first column must be named `sample` and contain the `<sample_name>` of the fastq files.

---

## 🌐 Environment 
The pipeline is designed to run locally using conda environments. You can create the base environment with:
```bash
conda env create --prefix=local/env/conda --file=local/env/env.yaml
```

The environment is directly activated by running `direnv allow` when prompted.

---

## 🚀 Running Snakemake
To run Snakemake rules, simply execute:

```bash
snakemake examplerule
```

---

## 📊 DGE analysis
For any DGE analysis, please set the `dge_tool` in the `config.yaml`. Then you should run snakemake  and ask for the file to be made in the DGE folder
    e.g. `snakemake -j 1  DGE/edger.toptable_clean.ALL_contrast.gz`

---

---
## 🧬 Isoform Discovery & Analysis Pipeline
The workflow includes an advanced pipeline for novel isoform assembly (`PsiCLASS`) and Differential Transcript Usage (`DRIMSeq`) analysis.
[👉 View the Isoform Pipeline Overview](README_Isoform.md)

---

## 🦠 Metagenomics Pipeline
In addition to standard RNA-seq analysis, this workflow includes a metagenomics module to discover and quantify microbial populations from unmapped reads. 

[👉 View the Metagenomics Pipeline Overview](README_Metagenome.md)

---

<details>
<summary><b>📦 Legacy: Running with Singularity</b></summary>
<br>

*Note: Singularity containers are not actively updated for all new components (such as the Metagenomics module). Conda is the recommended environment manager.*

### ⚙️ Singularity and Snakemake configuration
Before running the pipeline with Singularity, you must configure which host directories will be mounted inside the Singularity containers. Update the following profile configuration files accordingly:

#### 📦 Singularity (default profile)
Edit `profiles/default/config.yaml` and update the `singularity-args` field by listing all directories that need to be bound inside the Singularity container, for example:

```yaml
singularity-args: "-B /path/to/reference/genome -B /path/to/fastq"
```

### 🚀 Running Snakemake with Singularity
The pipeline can be executed locally, by running snakemake, while the individual rules are executed inside a Singularity container. 

***How it Works*** 
- ***Snakemake*** needs to be installed locally
- **Singularity container:** Defined in the Snakefile, it encapsulates the software dependencies required for running individual workflow rules.

You can avoid the `--use-conda` option when launching snakemake with singularity profiles.
</details>
