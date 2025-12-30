# FungiPhylogen 🧬🍄‍🟫🍄
[![DOI](https://zenodo.org/badge/547808488.svg)](https://doi.org/10.5281/zenodo.18097753)
[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A521.10.3-23aa62.svg?labelColor=000000)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

**FungiPhyloGen (FPG)** is a comprehensive, production-grade bioinformatics pipeline designed for genomic epidemiology and phylogenetic analysis of fungal pathogens. It is optimized for *Candida auris*, *Candida parapsilosis*, and *Cryptococcus* species.

Built on **Nextflow**, FPG leverages containerized **Conda** environments to ensure reproducibility and portability across different computing platforms (Local Laptop vs. HPC/Slurm).

## 🚀 Key Features

* **Scalable Parallelization (Scatter-Gather) Architecture:** Uses a scalable BCFtools-based variant calling workflow capable of handling 100+ samples efficiently.
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
FPG uses modular environments. You should create these in a central location (especially for HPC usage) and reference that path in `FungiPhyloGen.config`.

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

To ensure the pipeline runs smoothly on different infrastructures, you must configure the **Conda Paths** and **Filter Profiles** in `FungiPhyloGen.config`.

### 1. Tool Paths (Conda & Picard)
The pipeline requires absolute paths to the Conda environments and the Picard JAR file.
* **Local Runs:** The default is set to `${HOME}/anaconda3/envs`.
* **HPC Runs:** You **must** override the `condaCacheDir` in the `slurm` profile within nextflow.config i.e. `FungiPhyloGen.config`.

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
```

> **Tip:** Use absolute paths (e.g., `/home/user/data/...`) for read files to avoid "File not found" errors, especially when running on an HPC.


### 3. Run the Pipeline

**Option A: HPC Cluster (Slurm)**
* **Profile:** `-profile slurm`
* **Resources:** Uses the high-memory configuration (up to 128GB) and submits jobs to the batch queue.
* **Conda:** Uses the central Conda path defined in `nextflow.config`.

```bash
nextflow run main_fpg.nf \
    -c FungiPhyloGen.config \
    -profile slurm \
    -resume \
    --samplesheet ./samplesheet.csv \
    --filter_profile cauris_small \
    --prjName "FPG_Run01"
```

**Option B: Local Machine (Testing)**
* **Profile:** `-profile standard`
* **Resources:** Caps usage at 4 CPUs and 8GB RAM to prevent crashing your laptop.

```bash
nextflow run main_fpg.nf \
    -c FungiPhyloGen.config \
    -profile standard \
    --samplesheet ./samplesheet.csv \
    --filter_profile cauris_small
```
## 📂 Outputs
Results are saved in the `results/` folder (or directory specified by `--resultsDir`).

| Directory | Key File | Description |
| :--- | :--- | :--- |
| `trimmed/` | `*.fq.gz` | Cleaned FASTQ files. |
| `variants/` | `final_vcf/final.pass.vcf.gz` | **The Final Output.** High-quality SNPs only. |
| `snpeff_annotation/` | `snpeff_ann.vcf` | VCF annotated with gene changes (Synonymous/Missense). |
| `iqtree_phylogeny/` | `*.treefile` | Maximum Likelihood tree (Newick format). |
| `snp_distances/` | `snp_distance_matrix.tsv` | Pairwise SNP count matrix. |
| `multiqc/` | `multiqc_report.html` | Aggregate quality report. |


## 🔬 Downstream Analysis: Targeted Mutation Report (*C. auris*)

FungiPhyloGen includes helper scripts to perform targeted analysis of drug resistance genes (e.g., *ERG11*, *FKS1*). This workflow extracts missense variants from specific genes and generates a publication-ready CSV summary.

### 📋 Phase 1: Reference Data (Gene Map)

The analysis relies on a specific set of *Candida auris* Locus Tags. Create a file named `genelist.txt` containing the Locus Tags for the genes you wish to screen.

**Reference Gene Map (Clade I - B8441 Reference):**

| Gene Name | Locus Tag | Function / Description |
| :--- | :--- | :--- |
| **ERG11** | `B9J08_03698` | Azole resistance target (Lanosterol 14-alpha-demethylase) |
| **FKS1** | `B9J08_02922` | Echinocandin resistance target (1,3-beta-glucan synthase) |
| **FUR1** | `B9J08_01933` | Flucytosine resistance (Uracil phosphoribosyltransferase) |
| **TAC1b** | `B9J08_04780` | Transcription factor (Azole resistance) |
| **MRR1** | `B9J08_01918` | Transcription factor (Multi-drug resistance) |
| **CDR1** | `B9J08_02123` | ABC Transporter (Azole efflux pump) |
| **ERG3** | `B9J08_01595` | C-5 sterol desaturase |
| **MEC3** | `B9J08_00960` | DNA damage checkpoint protein |
| **PEA2** | `B9J08_03606` | Cell polarization protein |
| **FLO8** | `B9J08_02359` | Transcription factor (Filamentation) |
| **CIS2** | `B9J08_01093` | Gamma-glutamylcysteine synthetase |
| **rpsU** | `B9J08_05449` | Ribosomal protein S21 (Control/Housekeeping) |

### 🛠️ Phase 2: Execution

This process uses two scripts located in the `scripts/` directory:
1. `fpg_targeted_gene_mutations.sh`: Filters the Annotated VCF for missense variants in your target list.
2. `clean_report_v2.py`: Formats the output into a readable CSV, mapping Locus Tags to Gene Names.

#### Step 1: Create your Gene List
Create a plain text file listing the Locus Tags you want to check (one per line).

```bash
# Example: Create genelist.txt for ERG11, FKS1, rpsU etc - you can use the locus tags from the table above
nano genelist.txt
```

#### Step 2: Run the Extraction Script
This step requires the `fpgCallVariants` environment (or any environment with SnpSift and BCFtools).

```bash
# 1. Start an interactive session (if on HPC)
srun -J SNPann --pty bash
```

```bash
# 2. Activate the environment
conda activate fpgCallVariants
```

```bash
# 3. Run the extraction script
# Usage: bash fpg_targeted_gene_mutations.sh <Input_Ann_VCF> <Gene_List> <Output_VCF> <Output_Raw_CSV>
bash ~/github/FungiPhyloGen/scripts/fpg_targeted_gene_mutations.sh \
    results/snpeff_annotation/snpeff_ann.vcf \
    genelist.txt \
    results/targeted_output.vcf \
    results/targeted_raw.csv
```

#### Step 3: Generate Final Clean Report
Run the Python script to map the tags to gene names and produce the final summary.

```bash

conda activate fpgCallVariants

# Usage: python clean_report_v2.py <Input_Raw_CSV> > <Final_Report.csv>
python ~/github/FungiPhyloGen/scripts/clean_report_v2.py \
    results/targeted_raw.csv > results/resistance_report_final.csv
```


## 🆘 Troubleshooting

**1. "CRITICAL ERROR: Invalid Filter Profile"**
* **Cause:** You passed a profile name to `--filter_profile` that doesn't exist in the `filters` block of `nextflow.config`.
* **Fix:** Check `nextflow.config` for the list of valid keys (e.g., `cauris_small`).

**2. "Argument of file function cannot be null"**
* **Cause:** The samplesheet often has an empty trailing line or incorrect headers.
* **Fix:** Ensure headers are exactly `sampleID,read1,read2` and remove empty rows.

**3. "Picard JAR not found"**
* **Cause:** The pipeline is looking for `picard.jar` in the wrong location (e.g., local home dir) while running on Slurm.
* **Fix:** Check the `slurm` profile in `nextflow.config`. Ensure `params.CacheDir` is explicitly overridden to point to your HPC environment path.

**4. Pipeline ignores config / "Params are null"**
* **Cause:** Nextflow is not reading your configuration file because it is named `FungiPhyloGen.config` instead of the default `nextflow.config`.
* **Fix:** Rename the file: `mv FungiPhyloGen.config nextflow.config`, or explicitly pass it with `-c FungiPhyloGen.config`.


## 👥 Credits & Citations

<img src="https://github.com/Terra-Informatix-Pty-Ltd/Intro-To-Linux/blob/main/Wits-Mycology-Logo_Screenshot%202025-08-10%20115556.png?raw=true" align="right" width="150">
<img src="https://github.com/Terra-Informatix-Pty-Ltd/Intro-To-Linux/blob/main/BioInfoX-log.png?raw=true" align="right" width="150">

**FungiPhyloGen** was developed by **Stanford Kwenda**

### How to Cite
If you use this pipeline in your research, please cite this repository:
> Kwenda, S., Mwamba, T. M., Naicker, S., Maphanga, T. G., Nzimande, S. P., Jallow, S., & Govender, N. P. (2025). FungiPhyloGen: A Nextflow pipeline for fungal genomic epidemiology (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.18097754

### Software References
This pipeline leverages the following excellent open-source tools. Please cite them in your methods section:
* **Nextflow:** Di Tommaso, P., et al. (2017). *Nature Biotechnology*.
* **BCFtools:** Danecek, P., et al. (2021). *GigaScience*.
* **IQ-TREE:** Minh, B. Q., et al. (2020). *Molecular Biology and Evolution*.
* **SnpEff:** Cingolani, P., et al. (2012). *Fly*.
* **MultiQC:** Ewels, P., et al. (2016). *Bioinformatics*.
* **RapidNJ:** Simonsen, M., et al. (2011). *WABI*.

## 💰 Funding

This work was supported by **[]** under grant number **[]**.

* **Principal Investigator:** [Professor Nelesh Govender](https://www.wits.ac.za/people/academic-a-z-listing/g/neleshgovenderwitsacza/)
* **Grant Title:** []
