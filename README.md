# FungiPhyloGen 🧬🍄
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
* **Auditable SNP alignments:** Converts the final multi-sample PASS VCF directly into a canonical SNP alignment, retained-site manifest, alignment QC report, pairwise SNP-distance matrix, and pairwise comparable-site matrix. The converter does not apply additional biological filtering.
* **Resistance Profiling:** Automated annotation of resistance mutations in key genes (*ERG11*, *FKS1*, *FUR1*) using **SnpEff**.

## 🛠️ Installation

### 1. Prerequisites
Ensure you have the following installed on your system:
* [Conda](https://docs.conda.io/en/latest/) (or Mamba/Micromamba for speed)
* [Git](https://git-scm.com/)
* [Nextflow](https://www.nextflow.io/) and a Java runtime compatible with the installed Nextflow version

### 2. Clone the repository

For a reproducible installation, clone the release tag. The example uses `v1.2.1` and keeps the checkout directory named `FungiPhyloGen`.

```bash
mkdir -p "$HOME/github" &&
cd "$HOME/github" &&
git clone --branch v1.2.1 https://github.com/stanikae/FungiPhyloGen.git FungiPhyloGen &&
cd FungiPhyloGen
```

If `FungiPhyloGen` already exists, use a different parent directory or verify the existing checkout before reusing it. Do not overwrite an existing installation.

### 3. Setup Environments
FPG uses modular environments. To streamline installation, automated setup scripts are provided for both Linux and macOS. These scripts will automatically build all required Conda environments.

(Note: Environments should ideally be created in a central location for HPC usage, and that path must be referenced in `FungiPhyloGen.config`).

For Linux Users:
```bash
# Run the standard automated setup script
./setup_envs.sh
```

For macOS Users (Apple Silicon & Intel):
Because specialized bioinformatics tools require Intel-compiled packages, macOS users must use the OSX-specific script, which automatically handles architecture emulation.

```bash
# Make the script executable
chmod +x setup_envs_osx.sh

# Run the OSX automated setup script
./setup_envs_osx.sh
```

## ⚙️ Configuration

After cloning, open **`FungiPhyloGen.config`** and update the settings for the system that will execute the workflow. You can edit this file or make a site-specific copy. Set `CONFIG` in the launch template to the file you actually edited. Keep the config filename as supplied; pass it explicitly with `-c "$CONFIG"`.

### 1. Conda paths, tool paths, and Slurm queue

Run `conda env list` to locate the installed `fpg*` environments. Use the **absolute path to the parent directory** containing those environments, not the path to one individual environment.

- **Local (`standard`):** Set `condaCacheDir` in the top-level `params` block. The `standard` profile uses that path through `params.cacheDir`.
- **Cluster (`slurm`):** In `profiles { slurm { ... } }`, change **both** `process.queue` to your cluster's Slurm partition and `params.condaCacheDir` to the environment directory visible to compute nodes. The Slurm Conda path overrides the top-level setting. Keep `params.cacheDir = params.condaCacheDir`.

Edit these lines inside the **existing** `slurm` profile; do not replace the rest of the block:

```groovy
process.queue        = 'your_slurm_partition'
params.condaCacheDir = '/absolute/path/to/conda/envs'
params.cacheDir      = params.condaCacheDir
```

Verify that `PICARD` points to the installed Picard JAR and that `FASTQC_CONF`, `FASTQC_CONTAMINANTS`, and `FASTQC_ADAPTERS` resolve to the installed FastQC configuration directory and files. If the derived paths differ from the installation, update them in the settings used by the selected execution profile. For Slurm runs, the checkout, config, FASTQs, reference files, environments, and work and results directories must be accessible to the compute nodes. Start the Nextflow driver from a location permitted by your cluster.

### 2. Haploid population-SNP filter profiles

FungiPhyloGen applies species-appropriate population-SNP filters to the normalized multi-sample BCF. Select a profile at runtime with:

```bash
--filter_profile <profile>
```

| Profile | Population AF threshold | Description |
|---|---:|---|
| `general` | ≥ 0.25 | Default population-level haploid SNP profile |
| `cauris` | ≥ 0.25 | Canonical *Candida auris* population-SNP profile |
| `cauris_small` | ≥ 0.25 | Backward-compatible alias for `cauris` |
| `cauris_batch` | ≥ 0.25 | Backward-compatible alias for `cauris` |
| `wanomalus_small` | ≥ 0.75 | Conservative *Wickerhamomyces anomalus* population-SNP profile |
| `wanomalus_batch` | ≥ 0.85 | Ultra-high-purity *W. anomalus* population-SNP profile |
| `crypto` | AC ≥ 2 | Recurrent ALT alleles for conservative genomic epidemiology |
| `crypto_all` | None | Quality-filtered SNPs, including private variants, for exploratory phylogenetics |
| `crypto_common` | AF ≥ 0.25 | Analyses intentionally restricted to common ALT alleles |
| `cparapsilosis` | Species-specific | *Candida parapsilosis* profile using its dedicated depth-normalized expression |

For new *C. auris* analyses, use:

```bash
--filter_profile cauris
```

`cauris_small` and `cauris_batch` are backward-compatible aliases for `cauris`. All three profile names apply the same filtering expression.

#### Population-SNP criteria

The `general` and `cauris` profiles exclude a record when any of the following is true:

```text
TYPE!="snp"
MQ<40
QUAL<30
FS>60
MQSBZ<-2.0 or MQSBZ>2.0
INFO/AF<0.25
```

The *W. anomalus* profiles use the same core criteria with higher population allele-frequency thresholds:

```text
wanomalus_small: INFO/AF>=0.75
wanomalus_batch: INFO/AF>=0.85
```

These thresholds construct conservative, population-representative SNP sets. They are not intended for analyses that must retain rare, private, or low-frequency variants.

#### Population allele annotations

After normalization, FungiPhyloGen uses `bcftools +fill-tags` to materialize `INFO/AC`, `INFO/AN`, and `INFO/AF`. For haploid datasets, `AF = AC / AN`.

#### QUAL-by-depth change

FungiPhyloGen does not use merged-callset `QUAL/INFO/DP` as a filtering criterion. Samples are called independently and their BCF files are subsequently merged. The merged `QUAL` value is not recomputed as a joint-cohort confidence score, while `INFO/DP` increases with cohort depth; consequently, this ratio is cohort-depth dependent and is not equivalent to GATK QD.

FungiPhyloGen remains strictly BCFtools-based and uses explicit BCFtools annotations for site quality, mapping quality, strand bias, and population allele representation.

#### Missing-genotype filtering

The population-SNP profiles do not apply `F_MISSING`. The workflow uses `bcftools merge -0`, which can represent absent sample records as reference genotypes and therefore does not provide a reliable basis for cohort missingness filtering.

#### Upstream per-sample filtering

The preliminary per-sample filter is:

```text
TYPE="snp" && QUAL>=20 && INFO/DP>=5
```

This permissive filter is applied before merging. Final population-SNP selection is performed on the normalized multi-sample callset with the selected profile.

### Cryptococcus analysis profiles

Choose the profile by analysis purpose, not by cohort size. All three profiles apply the same SNP quality criteria: `MQ >= 40`, `QUAL >= 30`, `FS <= 60`, and `MQSBZ` between `-2` and `2` when present. `crypto` requires the same ALT allele in at least two called haploid isolates (`INFO/AC >= 2`). `crypto_all` does not apply an AC or AF threshold. `crypto_common` requires `INFO/AF >= 0.25` in the sampled cohort.

Select the mapping reference separately from the filter profile. Use a FASTA and matching GenBank annotation from the same assembly. The VNB and VNBII audits used H99 (VNI; not a VNB-matched reference). The VGIV audit used IND107 (*Cryptococcus tetragattii*, `GCF_000835755.1`). IND107 is not a universal reference for the *C. gattii* species complex.

For molecular-type investigations, compare query isolates with a labelled panel on a suitable shared reference coordinate system. Use phylogenetic placement and SNP distances together, and assess callable overlap before assigning a type. With sparse calls and `bcftools merge -0`, an absent record can be represented as REF; zero `N` bases does not establish callability. The current audits validate filter behaviour for the tested cohorts, not automated molecular-type assignment or transmission cutoffs. SNP distances from different filters or reference assemblies are not directly interchangeable. Use `--filter_profile crypto` for the recurrent rule previously tested under the interim `cneoformans` and `cgattii` keys.

The deprecated `cneoformans` key remains available for compatibility with the v1.2.0 `AF >= 0.25` behaviour. For new runs, use `crypto_common` to select that rule explicitly.

## 🏃 Usage

### 1. Check the workflow and runtime versions

Run from the environment where Nextflow is installed (for example, Conda `base`). Nextflow activates the task-specific environments during the run.

```bash
REPO="$HOME/github/FungiPhyloGen"
RELEASE="v1.2.1"
(
  set -e
  printf 'FungiPhyloGen version: %s\n' "$(cat "$REPO/VERSION")"
  CHECKOUT=$(git -C "$REPO" rev-parse HEAD)
  TAGGED=$(git -C "$REPO" rev-parse --verify "refs/tags/${RELEASE}^{commit}")
  printf 'Checkout commit: %s\nRelease commit:  %s\n' "$CHECKOUT" "$TAGGED"
  if [ "$(cat "$REPO/VERSION")" != "${RELEASE#v}" ] || [ "$CHECKOUT" != "$TAGGED" ]; then
    printf 'STOP: checkout does not match %s.\n' "$RELEASE" >&2
    exit 1
  fi
  git -C "$REPO" status --short
  nextflow -version
  java -version
  conda --version
)
```

`VERSION` reports the **FungiPhyloGen release**; `nextflow -version` reports the **Nextflow runtime**. The Git comparison verifies the committed revision, not local edits. Review the status output before launch. Local path edits in `FungiPhyloGen.config` are expected if you edited it directly; do not overlook changes to filter rules or workflow code.

### 2. Prepare the inputs

Create a run directory and a `samplesheet.csv` with one row per isolate, unique sample IDs, and absolute paths to both FASTQs. Use this exact header:

```csv
sampleID,read1,read2
Sample_A,/absolute/path/to/Sample_A_R1.fastq.gz,/absolute/path/to/Sample_A_R2.fastq.gz
Sample_B,/absolute/path/to/Sample_B_R1.fastq.gz,/absolute/path/to/Sample_B_R2.fastq.gz
```

Keep the reference FASTA and GenBank annotation from the **same assembly**. Set their paths explicitly in the launch command. The filter profile does not select a reference.

| Example analysis | Reference pair | `FILTER_PROFILE` |
|---|---|---|
| Known *C. auris* Clade I cohort | B8441 FASTA and its matching GenBank annotation | `cauris` |
| *Cryptococcus* VGIV cohort | IND107 FASTA and its matching GenBank annotation | `crypto` |

For other cohorts, choose the reference and filter independently. See [Cryptococcus analysis profiles](#cryptococcus-analysis-profiles) for reference scope and the `crypto_all` and `crypto_common` alternatives.

### 3. Launch on a local machine or Slurm

Edit the variables at the top, then paste the **entire block**, including the parentheses, into Bash or zsh. Use a separate run directory per dataset. For Slurm, complete the queue and Conda setup above before launching.

```bash
(
  set -e
  set -o pipefail

  REPO="$HOME/github/FungiPhyloGen"
  RUN="/absolute/path/to/run"
  CONFIG="$REPO/FungiPhyloGen.config"
  WORK="$RUN/work"
  OUT="$RUN/results"
  SHEET="$RUN/samplesheet.csv"
  REF="$RUN/ref/reference.fna"
  GBK="$RUN/ref/genomic.gbff"

  EXEC_PROFILE="slurm"
  FILTER_PROFILE="cauris"
  PROJECT_NAME="FPG_Run01"
  RUN_MODE="fresh"
  RESUME_SESSION=""

  command -v nextflow >/dev/null 2>&1 || {
    printf 'STOP: Nextflow is not available in this shell.\n' >&2
    exit 1
  }

  for FILE in "$REPO/main_fpg.nf" "$REPO/VERSION" "$CONFIG" "$SHEET" "$REF" "$GBK"; do
    if [ ! -s "$FILE" ]; then
      printf 'STOP: missing or empty file: %s\n' "$FILE" >&2
      exit 1
    fi
  done

  case "$EXEC_PROFILE" in
    standard|slurm) ;;
    *) printf 'STOP: choose standard or slurm.\n' >&2; exit 1 ;;
  esac

  RESUME_ARGS=()
  case "$RUN_MODE" in
    fresh) ;;
    resume)
      if [ ! -d "$RUN/.nextflow/cache" ] || [ ! -d "$WORK" ]; then
        printf 'STOP: resume requires the original Nextflow cache and work directory.\n' >&2
        exit 1
      fi
      RESUME_ARGS=(-resume)
      if [ -n "$RESUME_SESSION" ]; then
        RESUME_ARGS+=("$RESUME_SESSION")
      fi
      ;;
    *) printf 'STOP: RUN_MODE must be fresh or resume.\n' >&2; exit 1 ;;
  esac

  if [ -e "$OUT" ] || [ -L "$OUT" ]; then
    printf 'STOP: choose a new results directory: %s\n' "$OUT" >&2
    exit 1
  fi

  cd "$RUN"
  CHECK=$(mktemp -d "$RUN/run-logs.XXXXXX")
  printf 'FungiPhyloGen: %s\nLogs: %s\nResults: %s\n' \
    "$(cat "$REPO/VERSION")" "$CHECK" "$OUT"

  NXF_SYNTAX_PARSER=v1 nextflow run "$REPO/main_fpg.nf" \
    -c "$CONFIG" \
    -profile "$EXEC_PROFILE" \
    "${RESUME_ARGS[@]}" \
    -work-dir "$WORK" \
    --inputDir "$RUN" \
    --resultsDir "$OUT" \
    --samplesheet "$SHEET" \
    --refseq "$REF" \
    --gbk "$GBK" \
    --filter_profile "$FILTER_PROFILE" \
    --prjName "$PROJECT_NAME" \
    -with-trace "$CHECK/trace.txt" \
    -with-report "$CHECK/report.html" \
    2>&1 | tee "$CHECK/run.log"
)
```

| Setting | What to use |
|---|---|
| `EXEC_PROFILE` | `standard` for local execution or `slurm` for a configured cluster. This determines where tasks run. |
| `FILTER_PROFILE` | Selects the SNP-filtering rule, independently of the execution profile and reference. |
| `RUN_MODE` | `fresh` for a new run. For `resume`, reuse the original launch directory, work directory, and input paths. Keep the original settings unless intentionally changing them. |
| `RESUME_SESSION` | Leave empty to resume the last run. If the launch directory contains multiple runs, use `nextflow log` from that directory to identify the session to resume and set its ID here. |
| `OUT` | A new published-results directory for this attempt. For a resume, change it to a new name, such as `$RUN/results-resumed`. |
| `CHECK` | Created automatically. It stores `run.log`, `trace.txt`, and `report.html` separately from published results and the work cache. |

### 4. Monitor and verify the run

At startup, confirm that the banner shows the expected **FASTA, GenBank annotation, filter profile, and exclusion expression**. In a separate terminal, use the directory printed under `Logs:`:

```bash
tail -f "/path/printed/under/Logs/run.log"
```

`Ctrl+C` in the monitoring terminal stops `tail`, not Nextflow. If the workflow fails, check `run.log` and the driver log (`.nextflow.log` in the launch directory). `tee` copies output to the terminal and log; `pipefail` preserves a failing Nextflow exit status. When pasting into an interactive shell, the surrounding parentheses keep `set -e` scoped to the launch block. A GNU Screen window started directly with a script can close when that script exits, so open an interactive Screen shell first if you want to retain the prompt after the run.

Check the trace for completed tasks and, on resumed runs, for tasks marked `CACHED`. Resume does not guarantee that all tasks will be reused. After completion, review the final PASS VCF and the alignment QC and distance outputs described below. Keep the work directory while you still need to resume or access published symlinks.

## 📂 Outputs
Results are saved in the `results/` folder (or directory specified by `--resultsDir`).

| Directory | Key File | Description |
| :--- | :--- | :--- |
| `trimmed/` | `*.fq.gz` | Cleaned FASTQ files. |
| `variants/` | `final_vcf/final.pass.vcf.gz` | **The Final Output.** High-quality SNPs only. |
| `snpeff_annotation/` | `snpeff_ann.vcf` | VCF annotated with gene changes (Synonymous/Missense). |
| `variants/msa/` | `fpg_snp_aln.fa` | Canonical SNP-only alignment used by IQ-TREE, RapidNJ, and SNP-distance analysis. |
| `variants/msa/` | `fpg_snp_sites.tsv` | Alignment-column provenance, genomic coordinates, reconstructed alternate alleles, source-record counts, and ambiguity counts. |
| `variants/msa/` | `fpg_snp_alignment_qc.tsv` | Alignment dimensions, source-record totals, reconstructed split positions, and per-taxon ambiguous-base counts. |
| `variants/msa/` | `fpg_snp_distances.tsv` | Raw pairwise SNP-difference matrix. |
| `variants/msa/` | `fpg_snp_comparable_sites.tsv` | Pairwise number of comparable A/C/G/T alignment sites used to interpret each SNP distance. |
| `iqtree_phylogeny/` | `*.treefile` | Maximum Likelihood tree (Newick format). |
| `snp_distances/snpdists/` | `FungiPhylogen.snpdist.csv` | Pairwise SNP count matrix. |
| `multiqc/` | `multiqc_report.html` | Aggregate quality report. |

## SNP alignment and distance construction

FungiPhyloGen constructs its canonical downstream SNP alignment directly from `final.pass.vcf.gz`. All biological variant-quality filtering is completed upstream by the `bcftools` processes and the selected filter profile. The alignment converter does not reapply `QUAL`, `MQ`, `FS`, `MQSBZ`, depth, allele-count, allele-frequency, or missingness thresholds.

The converter reconstructs normalized split multiallelic records into one biological genomic position, preserves explicit missing haploid genotypes as `N`, rejects non-haploid or contradictory genotype states, and records every retained column in `fpg_snp_sites.tsv`. Pairwise distances exclude comparisons involving `N`; `fpg_snp_comparable_sites.tsv` reports the denominator for every sample pair.

The same canonical `fpg_snp_aln.fa` is supplied to IQ-TREE and RapidNJ. The generated `fpg_snp_distances.tsv` matrix is validated against `snp-dists` semantics.



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

**1. Invalid filter profile**

Use a key from `params.filters` in `FungiPhyloGen.config`, such as `cauris` or `crypto`. The filter profile is separate from the `standard` or `slurm` execution profile.

**2. A file argument is null or an input cannot be found**

Check the CSV header (`sampleID,read1,read2`), remove empty rows, and confirm the FASTQs and reference pair are readable. For Slurm, verify the paths are visible to compute nodes, not just the login node.

**3. Picard JAR or FastQC support files cannot be found**

Check `PICARD` and the `FASTQC_*` settings in the selected execution profile. For Slurm, confirm that `params.condaCacheDir` points to the correct environment directory and that the derived file paths match the installation.

**4. Configuration did not load or Slurm jobs fail to submit**

Pass the edited `FungiPhyloGen.config` (or its site-specific copy) with `-c "$CONFIG"`. Verify the config path in the startup banner. For Slurm, check that `process.queue` names a partition you can submit to and review the driver log for the scheduler error.

## 👥 Credits & Citations

<img src="https://github.com/Terra-Informatix-Pty-Ltd/Intro-To-Linux/blob/main/Wits-Mycology-Logo_Screenshot%202025-08-10%20115556.png?raw=true" align="right" width="150">
<img src="https://github.com/Terra-Informatix-Pty-Ltd/Intro-To-Linux/blob/main/BioInfoX-log.png?raw=true" align="right" width="150">

**FungiPhyloGen** was developed by [**Stanford Kwenda**](https://www.linkedin.com/in/stanford-kwenda-7bb87530/)

### How to Cite
If you use this pipeline in your research, please cite this repository:
> Kwenda, S., Mwamba, T. M., Naicker, S., Maphanga, T. G., Nzimande, S. P., Jallow, S., & Govender, N. P. (2026). FungiPhyloGen: A Nextflow pipeline for fungal genomic epidemiology (v1.1.0). Zenodo. https://doi.org/10.5281/zenodo.18121445

### Software References
This pipeline leverages the following excellent open-source tools. If you use FungiPhyloGen, please ensure you also cite the underlying tools in your methods section:

* **Nextflow:** Di Tommaso, P., et al. (2017). Nextflow enables reproducible computational workflows. *Nature Biotechnology*.
* **FastQC:** Andrews, S. (2010). FastQC: A Quality Control Tool for High Throughput Sequence Data.
* **Trim Galore:** Krueger, F. (2015). Trim Galore.
* **BWA:** Li, H., & Durbin, R. (2009). Fast and accurate short read alignment with Burrows-Wheeler transform. *Bioinformatics*.
* **SAMtools & BCFtools:** Danecek, P., et al. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*.
* **SnpEff:** Cingolani, P., et al. (2012). A program for annotating and predicting the effects of single nucleotide polymorphisms, SnpEff. *Fly*.
* **snp-dists:** Seemann, T. (2020). snp-dists: Pairwise SNP distance matrix from a FASTA sequence alignment.
* **IQ-TREE:** Minh, B. Q., et al. (2020). IQ-TREE 2: New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era. *Molecular Biology and Evolution*.
* **RapidNJ:** Simonsen, M., et al. (2011). Rapid neighbour-joining. *WABI*.
* **MultiQC:** Ewels, P., et al. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*.

## 💰 Funding

This work was supported by **[]** under grant number **[]**.

* **Principal Investigator:** [Professor Nelesh Govender](https://www.wits.ac.za/people/academic-a-z-listing/g/neleshgovenderwitsacza/)