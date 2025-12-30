# FungiPhylogen 🧬🍄‍🟫🍄
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

**FungiPhyloGen (FPG)** is a comprehensive, production-grade bioinformatics pipeline designed for genomic epidemiology and phylogenetic analysis of fungal pathogens. It is optimized for *Candida auris*, *Candida parapsilosis*, and *Cryptococcus* species.

Built on **Nextflow**, FPG leverages containerized **Conda** environments to ensure reproducibility and portability across different computing platforms (Local Laptop vs. HPC/Slurm).

## 🚀 Key Features

* **Scatter-Gather Architecture:** Uses a scalable BCFtools-based variant calling workflow capable of handling 100+ samples efficiently.
* **Robust Quality Control:** Automated adapter trimming (TrimGalore) and quality assessment (FastQC, MultiQC).
* **Scientific Filtering:** Implements a **"Filter Profile"** strategy:
    * **Hard Filter:** Removes low-depth noise (`DP < 5`) per sample to optimize merging.
    * **Soft Filter:** Applies complex, organism-specific quality thresholds (e.g., `MQ`, `QUAL`, `FS`) via config profiles.
* **Phylogenetics:** Generates Maximum Likelihood trees (**IQ-TREE**) and Neighbor-Joining trees (**RapidNJ**) with bootstrap support.
* **Resistance Profiling:** Automated annotation of resistance mutations in key genes (*ERG11*, *FKS1*, *FUR1*) using **SnpEff**.

## 🛠️ Installation

### 1. Prerequisites
Ensure you have the following installed on your system:
* [Conda](https://docs.conda.io/en/latest/) (or Mamba/Micromamba for speed)
* [Git](https://git-scm.com/)
* [Nextflow](https://www.nextflow.io/)

### 2. Clone Repository
```bash
mkdir -p $HOME/github
cd $HOME/github

git clone https://github.com/stanikae/FungiPhyloGen.git
cd FungiPhyloGen
```

### 3. Setup Environments
FPG uses modular environments. You should create these in a central location (especially for HPC usage) and reference that path in nextflow.config.

```bash
conda env create --file lib/fpgtrimReads.yml --solver=libmamba -y
conda env create --file lib/fpgtrimReads.yml --solver=libmamba -y
conda env create --file lib/fpgDenovo.yml --solver=libmamba -y
conda env create --file lib/vcftools.yml --solver=libmamba -y
conda env create --file lib/align.yml --solver=libmamba -y
conda env create --file lib/fpgVcf2FastaEnv.yml --solver=libmamba -y
conda env create --file lib/vcfkit.yml --solver=libmamba -y
conda env create --file lib/callVar.yml --solver=libmamba -y
conda env create --file lib/phylo.yml --solver=libmamba -y
```

## ⚙️ Configuration (Crucial)

To ensure the pipeline runs smoothly on different infrastructures, you must configure the **Conda Paths** and **Filter Profiles** in `nextflow.config`.

### 1. Tool Paths (Conda & Picard)
The pipeline requires absolute paths to the Conda environments and the Picard JAR file.
* **Local Runs:** The default is set to `${HOME}/anaconda3/envs`.
* **HPC Runs:** You **must** override the `condaCacheDir` in the `slurm` profile within nextflow.config i.e. `Fungiphylogen.config`.

```groovy
// In nextflow.config -> profiles -> slurm
slurm {
    // ...
    // UPDATE THIS PATH to your central HPC conda environment location e.g.
    params.CacheDir = "/scratch/package/anaconda/anaconda3_2023-02-27/envs" 
    
    // This ensures the pipeline finds the correct Picard JAR
    params.PICARD = "${params.CacheDir}/fpgAlign/share/picard-2.27.4-0/picard.jar"
}
```

### 2. Filter Profiles
FungiPhyloGen uses pre-defined profiles to switch between strict and lenient filtering logic. Use the `--filter_profile` flag at runtime.

| Profile | Target Organism | Filter Logic & Use Case |
| :--- | :--- | :--- |
| `general` | **All Fungi** | **Recommended Default.** Balanced filters with purity check (AF > 0.75). Use for *Histoplasma*, *Emergomyces*, etc. |
| `cauris_small` | *C. auris* | **Standard Purity.** Identical to `general`. |
| `cauris_batch` | *C. auris* | **Core Consensus.** Allows heterozygous calls but removes sites with >20% missing data. |
| `wanomalus_small` | *W. anomalus* | **High Stringency.** Very strict Quality/Depth ratio (QD > 0.1). |
| `wanomalus_batch` | *W. anomalus* | **Ultra-High Purity.** Requires AF > 0.85. Ideal for defining reference SNPs. |
| `cparapsilosis` | *C. parapsilosis* | **Custom Metric.** Uses specialized `QUAL/MAX(AD)` ratio. |
| `cneoformans` | *C. neoformans* | **Standard Purity.** Identical to `general`. |


## 🏃 Usage

### 1. Prepare Reference Genome
You need a FASTA file and a GenBank format annotation file for your target organism.
* Example: `ref/Clade_I_B8441.fna` and `ref/Clade_I_B8441.gbff`.

### 2. Prepare Samplesheet
Create a CSV file (e.g., `samplesheet.csv`) containing your sample metadata and read paths. 

**Crucial:** The headers must be exactly `sampleID`, `read1`, and `read2`.

```csv
sampleID,read1,read2
Sample_A,/full/path/to/Sample_A_R1.fastq.gz,/full/path/to/Sample_A_R2.fastq.gz
Sample_B,/full/path/to/Sample_B_R1.fastq.gz,/full/path/to/Sample_B_R2.fastq.gz
