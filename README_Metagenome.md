# 🦠 Metagenomics Pipeline Overview

The primary objective of this pipeline is to discover and quantify microbial populations (such as bacteria, viruses, and fungi) present in biological samples, leveraging data that is already generated during standard sequencing protocols.

While sequencing experiments (like RNA-seq) are often designed to study a host organism, such as a human patient or a mouse model, these samples simultaneously capture genetic material from the *microbiome*---the community of microorganisms residing in or on the host. This pipeline takes advantage of the "leftover" sequencing reads that fail to align to the host genome, repurposing them to characterize the microbial composition and estimate their relative abundances.

### 🐍 How It Works
The pipeline is based on *Snakemake*, a workflow management system designed for reproducible and scalable data analyses. In Snakemake, the entire analysis is broken down into modular components called _rules_. Each rule defines a specific processing step along with its required input files and expected output files. 

By automatically tracking the dependencies between these files, Snakemake connects the rules into a continuous chain. This ensures that every tool is executed in the exact correct order, seamlessly linking them together into a streamlined, fully automated workflow.

The overall process follows a structured, step-by-step path to filter, classify, and quantify microbial taxa from the raw sequencing data.
---
## 🚀 Pipeline at a Glance
```mermaid
flowchart TD
linkStyle default stroke:#333,stroke-width:2px,color:#000;
    %% Core Pipeline
    A[📦 Raw FASTQ Reads] --> B(🛠️ <b>fastp</b>: QC & Adapter Trimming)
    B --> C(🧬 <b>STAR</b>: Host Read Alignment & Separation)
    C --> D[🔬 Unmapped non-host reads]
    D --> E(🗂️ <b>Kraken2</b>: Initial taxonomic classification)
    E --> F(🧹 <b>KrakenTools</b>: Residual host removal)
    F --> G(🗂️ <b>Kraken2</b>: Second pass classification <br> NCBI PlusPF + GTDB)
    G --> H(📊 <b>Bracken</b>: Abundance re-estimation & noise filtering)
    H --> I[(📈 Cohort-wide Abundance Matrix)]
    
    %% Forking paths for downstream analysis
    I --> J(🎨 <b>Krona</b> + <b>R</b>: Interactive plots)
    I --> K(📉 <b>R Stats</b>: Feature filtering & Diff. Abundance)
    K --> L((🎯 Significant Microbial Biomarkers))

    %% Structural Anchor: Tethering to the subgraph container itself
    L ~~~ Validation

    %% Optional Validation Module enclosed in a subgraph
    subgraph Validation ["Targeted Validation Module (Optional)"]
        direction LR
        M(✂️ <b>KrakenTools</b>: Taxon-specific read extraction) --> N(🗺️ <b>Minimap2</b>: Alignment against reference)
        N --> O(👀 <b>IGV</b>: Genomic coverage visual inspection)
        O --> P((✅ Verified Microbial Presence))
    end
    
    %% Styling - Explicitly locking text to black (color:#000) for Dark Mode compatibility
    classDef default fill:#f9f9f9,stroke:#333,stroke-width:1px,color:#000;
    classDef matrix fill:#e1f5fe,stroke:#0288d1,stroke-width:2px,color:#000;
    classDef endpoint fill:#e8f5e9,stroke:#388e3c,stroke-width:2px,color:#000;
    classDef toolNode fill:#fff8e1,stroke:#ffb300,stroke-width:2px,color:#000;
    
    %% Applying Classes
    class I matrix;
    class L,P endpoint;
    class B,C,E,F,G,H,J,K,M,N,O toolNode;
```



<details>
<summary><b>📖 Click to read detailed pipeline step descriptions</b></summary>

<br>

**1. Quality Control (`fastp`)**
* Before biological analysis, the raw sequencing data undergoes rigorous quality control. *fastp* ensures the removal of low-quality bases and technical artifacts (such as adapter sequences) to provide a clean dataset for downstream alignment.

**2. Host Read Separation (`STAR`)**
* **Host Alignment:** Because the initial sequencing targeted the host organism (e.g., Human or Mouse), the vast majority of the reads belong to the host. *STAR* aligns all the sequences against the respective host reference genome. Any read that successfully maps to the host is separated and set aside for standard transcriptomic analysis.
* **Microbiome Isolation:** What remains is a collection of *unmapped* reads. This leftover data serves as the core input for the metagenomic analysis, as it contains the genetic signatures of the microbiome, alongside potential uncharacterized biological elements.

**3. Initial Taxonomic Classification (`Kraken2`)**
* **The Process:** Once the host data is depleted, the remaining unmapped sequences are queried by *Kraken2* against comprehensive databases of known microbial genomes to achieve a robust taxonomic classification.
* **Dual-Database Strategy:** To ensure the highest accuracy and mitigate the risk of false positives, the pipeline cross-references the data against two distinct databases:
  * **NCBI PlusPF Database:** A broad, standard database covering bacteria, archaea, viruses, and fungi. *(Limitation: Due to the over-representation of clinical isolates in public repositories, it is prone to erroneously assigning novel environmental sequences to well-studied taxa—database bias).*
  * **GTDB (Genome Taxonomy Database):** A highly standardized database focused specifically on high-quality bacterial and archaeal genomes. *(Limitation: It relies exclusively on prokaryotic phylogenies, making it inherently blind to viruses and fungi).*

**4. Secondary Host Depletion (`KrakenTools`)**
* **Flagging Residuals:** Despite the initial alignment step, some residual host reads can still escape detection. During the first Kraken2 classification against the NCBI PlusPF database, any sequence that is classified as _Homo sapiens_ (TaxID: 9606) is explicitly flagged.
* **Extraction:** These residual host reads are computationally extracted and removed from the unmapped dataset.

**5. Second-Pass Classification (`Kraken2`)**
* Kraken2 is executed a second time on this strictly depleted dataset. This rigorous two-pass approach ensures that the final microbial profile is not skewed by human sequences.

**6. Abundance Estimation & Noise Filtering (`Bracken`)**
* **Resolving Ambiguity:** After the initial taxonomic classification, the pipeline estimates the precise *relative abundance* of the identified microbial populations. Because Kraken2 is a conservative classifier, reads that share sequence similarities with multiple species are assigned to a higher taxonomic rank (such as the genus or family level) rather than guessing a specific species. 
* **Probabilistic Reallocation:** To resolve this, *Bracken* employs a Bayesian probabilistic model to evaluate these ambiguously classified reads. It mathematically redistributes them down to the species level based on the relative abundance of reads that were unambiguously assigned. 
* **Noise Reduction:** Concurrently, *Bracken* applies essential noise filtering by discarding taxa that appear below a defined abundance threshold. This combined approach of probabilistic reallocation and noise reduction ensures that the final output reports only confident, reliable microbial identifications rather than background computational artifacts.

**7. Data Aggregation & Differential Abundance (Custom Scripts + R)**
* **Matrix Generation:** Once the individual microbial profiles are refined, the pipeline aggregates the data across the entire cohort. It merges the sample-level abundances into a unified, cohort-wide matrix.
* **Statistical Analysis:** To identify biologically meaningful differences, the pipeline applies stringent prevalence and expression filters to remove sparse taxa. Finally, it integrates experimental metadata and performs statistical testing (e.g., Wilcoxon rank-sum test) to pinpoint which microbial populations are significantly enriched or depleted between different biological conditions.

**8. Interactive Visualizations (`Krona` & `R/ggplot2`)**
* *Krona* renders dynamic, multi-layered pie charts that allow researchers to visually explore the microbiome. 
* Customized *R/ggplot2* scripts generate publication-ready graphics to summarize the cohort, including:
  * **Top 10 Taxa:** Summary barplots highlighting the most abundant microbes across all samples.
  * **Abundance Heatmaps:** Visualizations to identify microbial clustering patterns across the cohort.
  * **Read Tracking:** Stacked barplots detailing the proportion of initial raw reads, unmapped reads, and successfully classified reads.

**9. Targeted Genomic Validation (Optional)**
* **Read Extraction:** To confirm the biological presence of a specific microbe of interest (e.g., a suspected pathogen), the pipeline includes an orthogonal validation module. First, reads that were computationally assigned to the taxon of interest are explicitly extracted and converted back into FASTQ format using *KrakenTools*.
* **Re-alignment & Inspection:** These targeted reads are then rigorously aligned against the specific reference genome of that species using *Minimap2*. Finally, the resulting alignments are merged and prepared for visualization in *IGV*. This allows researchers to visually inspect the genomic coverage and physically confirm whether the reads map evenly across the microbe's genome, ruling out localized sequence artifacts or false positives.

</details>

---
### 🧩 Alternative Profiling Strategies (Optional)
While the core pipeline relies on k-mer based classification with *Kraken2*, it also natively integrates alternative classifiers to provide rigorous cross-validation and orthogonal profiling:

* **`MetaPhlAn4`**: Employs a curated database of unique clade-specific marker genes, providing highly precise species-level identification for bacteria and archaea.
* **`Kaiju`**: Translates unmapped sequencing reads into amino acids and classifies them at the protein level. This approach is highly effective for identifying divergent or novel organisms that may evade standard nucleotide-based detection.



## ⚙️ Technical Requirements & Database Sizes

This pipeline relies on several memory-intensive steps, particularly during host alignment and taxonomic classification. Below are the estimated sizes for the required databases and the recommended hardware configurations.

### Core Database Sizes & RAM Requirements

| Tool | Component / Database | Disk Space (Index/DB) | Expected RAM Usage |
| :--- | :--- | :--- | :--- |
| **`STAR`** | Host Genome (e.g., GRCh38/mm10) | ~27 GB | ~30 GB |
| **`Kraken2`** | NCBI PlusPF Database | ~80–104 GB | ~128 GB |
| **`Kraken2`** | GTDB v226 Database | ~650 GB | ~700 GB |
| **`Kraken2`** | NCBI core_nt Database (Optional) | ~316 GB | ~320 GB |
| **`MetaPhlAn4`** | `mpa_vJun23_CHOCOPhlAnSGB` (Optional)| ~19 GB | ~24 GB |
| **`Kaiju`** | `nr_euk` Database (Optional) | ~40–50 GB | ~60 GB |
| **`Minimap2`** | Microbial Reference (per species) | < 1 GB | ~4 GB |

### General Tool Requirements
* **`fastp`, `KrakenTools`, `Bracken`, `Krona`, `R Stats`**: These utilities are relatively lightweight and will run comfortably on standard compute nodes with **8–16 GB RAM**.
* **Storage**: In addition to the database sizes listed above, ensure you have sufficient temporary disk space to handle intermediate fastq files (unmapped reads) and `.sam`/`.bam` alignments.

## 🔬 Visualizing Results with IGV (Remote Server Setup)

Once the targeted validation module successfully maps extracted reads against a reference genome, you can visually inspect the `.bam` alignments. If your pipeline runs on a headless remote server, you can securely stream the files to your local browser or desktop app without downloading massive BAM files.

### Step 1: Open a Secure SSH Tunnel
Bind a local port on your machine to the remote server to route the data securely. Run this on your **local machine's terminal**:

```bash
ssh -L 8080:localhost:8080 user@server_address
```

### Step 2: Start the Smart Python Server
IGV relies on "Byte-Range" requests to stream only the specific genomic coordinates you are viewing. The default Python HTTP server will crash when attempting this. Instead, use `RangeHTTPServer`.

Navigate to the root directory of your pipeline on the **remote server** and run:

```bash
# Install if not already present: pip install RangeHTTPServer
python3 -m RangeHTTPServer 8080
```
*(Tip: You can run this inside a `tmux` session to keep the server alive in the background permanently).*

### Step 3: Stream to IGV
Open the **[IGV Web App](https://igv.org/app/)** (or the IGV Desktop App) on your local machine.

**1. Load the Reference Genome (and Annotations):**
* Click **Tracks** > **URL...** (or *File > Load from URL* in Desktop).
* **Fasta URL:** `http://localhost:8080/Resources/genomes/<TaxID>/genome.fna`
* **Index URL:** `http://localhost:8080/Resources/genomes/<TaxID>/genome.fna.fai`
* *(Optional)* **Annotation GFF:** `http://localhost:8080/Resources/genomes/<TaxID>/annotation.gff`

**2. Load the Read Alignments:**
* Click **Tracks** > **URL...**
* **Track URL (BAM):** `http://localhost:8080/alignments_merged/<TaxID>/merged_all_samples.bam`
* **Index URL (BAI):** `http://localhost:8080/alignments_merged/<TaxID>/merged_all_samples.bam.bai`

### Step 4: Locating High-Coverage Contigs
If your reference genome is highly fragmented (e.g., hundreds of contigs), do not guess where your reads landed. Use `samtools` on your server to instantly rank the contigs by read depth:

```bash
samtools idxstats alignments_merged/<TaxID>/merged_all_samples.bam | sort -k3,3nr | head -n 10
```
