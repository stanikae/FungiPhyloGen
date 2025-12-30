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
git clone [https://github.com/stanikae/FungiPhyloGen.git](https://github.com/stanikae/FungiPhyloGen.git)
cd FungiPhyloGen


### 3. Setup Environments
FPG uses modular environments. You should create these in a central location (especially for HPC usage) and reference that path in nextflow.config.

```
# Example: Creating environments in a central directory
# Adjust path /spaces/stanford/anaconda3/envs to your preference
conda env create -f lib/fpgtrimReads.yml --prefix /spaces/stanford/anaconda3/envs/fpgtrimReads
conda env create -f lib/align.yml --prefix /spaces/stanford/anaconda3/envs/fpgAlign
# ... repeat for other yml files in lib/
```
