# SunClock: an epigenetic age estimator trained in young individuals

This repository provides the official model coefficients and an optimized R pipeline to calculate DNA methylation age (DNAm age) and Age Acceleration using the **SunClock** and **SunClock-CP** epigenetic clocks. 

- **SunClock**: Developed using Whole-Genome Bisulfite Sequencing (WGBS) data.
- **SunClock-CP**: Tailored for traditional methylation microarrays (Illumina 450K/EPICv1/EPICv2/MSA), fully compatible with standard `cg` probe nomenclature.

---

## Repository Structure

- `SunClock_coef.csv`: Coefficient table for the WGBS model.
- `SunClock-CP_coef.csv`: Coefficient table for the Array model.
- `SunClock_script.R`: Main source script containing the core pipeline function.

---

## Input Data Format

### 1. Methylation Beta Matrix (`beta_matrix`)
A numeric matrix or data frame of DNA methylation beta values (ranging from 0 to 1).
- **Rows**: Sample IDs (set as `rownames`).
- **Columns**: CpG site identifiers. 
  - For `model_type = "WGBS"`, column names must match genomic coordinates (e.g., `chr12:21771411:21771412`).
  - For `model_type = "array"`, column names must match Illumina probe IDs (e.g., `cg00034076`).

> 🧬 **Genome Assembly Note (WGBS)**
> In `SunClock_coef.csv`, the default genomic coordinates in the `CpG_pos` column are based on the **hg38** reference genome. If your dataset is aligned to **hg19**, please use the provided `genome_version = "hg19"` option. *(Note: A minimal number of clock CpGs failed the liftOver conversion to hg19 and are marked as `NA`.)*

### 2. Phenotype Metadata (`pheno_df`)
*(Optional, only required if computing Age Acceleration)*
A data frame where sample names are provided as a standard column instead of row names. It **must** contain the following two columns:
- `Sample_ID`: Character string matching the `rownames` of your `beta_matrix`.
- `Age`: Numeric values representing the chronological age.
- *Note*: Any additional clinical or biological columns (e.g., `Sex`, `Group`, `Cell_Types`) will be automatically preserved and aligned in the final output.

> ⚠️ **Important Note on Data Preprocessing**
> This pipeline handles missing or `NA` features by assigning them a value of `0` (meaning they contribute nothing to the cumulative score). To prevent severe underestimation of epigenetic age due to high data missingness, it is **strongly recommended** that users perform proper baseline imputation (e.g., KNN or median imputation) on their beta matrix prior to running this script.

---

## Quick Start & Usage

Ensure that `SunClock_script.R` and the coefficient CSV files are placed in your working directory, then execute the following implementation in R:

```R
# Source the pipeline function
source("SunClock_script.R")

# Required Inputs Format:
# 1. my_betas : Matrix/DataFrame where RowNames = Sample IDs, ColNames = CpG IDs (cg... or chr...)
# 2. my_pheno : DataFrame containing columns "Sample_ID" and "Age" (plus any optional clinical columns)
# 3. CSV files: Keep "SunClock-CP_coef.csv" or "SunClock_coef.csv" in your working directory, or specific coef_dir
# --- Example 1: Basic DNAm Age Prediction (Default Array model, SunClock-CP) ---
final_output <- predict_epigenetic_age(beta_matrix = my_betas, model_type = "array")

# --- Example 2: WGBS Prediction using hg19 Coordinates (SunClock hg19) ---
final_output <- predict_epigenetic_age(
  beta_matrix    = my_betas,
  model_type     = "WGBS",
  genome_version = "hg19"
)

# --- Example 3: Full Pipeline (Calculate WGBS hg38 Age + Age Acceleration + Keep metadata) ---
final_output <- predict_epigenetic_age(
  beta_matrix    = my_betas,
  model_type     = "WGBS",
  genome_version = "hg38",
  coef_dir       = ".",
  compute_accel  = TRUE,
  pheno_df       = my_pheno
)

head(final_output)
