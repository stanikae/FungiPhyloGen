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
