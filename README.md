# SunClock: an epigenetic age estimator trained in young individuals

This repository provides the official model coefficients and an optimized R pipeline to calculate DNA methylation age (DNAm age) and Age Acceleration using the **SunClock** and **SunClock-CP** epigenetic clocks. 

- **SunClock**: Developed using Whole-Genome Bisulfite Sequencing (WGBS) data, featuring single-nucleotide resolution.
- **SunClock-CP**: Tailored for traditional methylation microarrays (Illumina 450K/EPICv1/EPICv2/MSA), fully compatible with standard `cg` probe nomenclature.

---

## Repository Structure

- `SunClock_coef.csv`: Coefficient table for the WGBS model.
- `SunClock-CP_coef.csv`: Coefficient table for the Array model.
- `predict_epigenetic_age.R`: Main source script containing the core pipeline function.

---

## Input Data Format

### 1. Methylation Beta Matrix (`beta_matrix`)
A numeric matrix or data frame of DNA methylation beta values (ranging from 0 to 1).
- **Rows**: Sample IDs (set as `rownames`).
- **Columns**: CpG site identifiers. 
  - For `model_type = "WGBS"`, column names must match genomic coordinates (e.g., `chr12:21771411:21771412`).
  - For `model_type = "array"`, column names must match Illumina probe IDs (e.g., `cg00034076`).

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

Ensure that `predict_epigenetic_age.R` and the coefficient CSV files are placed in your working directory, then execute the following implementation in R:

```R
# Source the pipeline function
source("predict_epigenetic_age.R")

# Load your cohort data (replace with your actual data objects)
# my_betas <- as.matrix(your_imputed_beta_data)
# my_pheno <- your_phenotype_dataframe

# --- Scenario A: Basic DNAm Age Prediction (Default Array Model) ---
output_basic <- predict_epigenetic_age(
  beta_matrix = my_betas, 
  model_type  = "array"
)

# --- Scenario B: Full Pipeline (DNAm Age + Age Acceleration + Metadata Preservation) ---
output_full <- predict_epigenetic_age(
  beta_matrix   = my_betas, 
  model_type    = "array",      # Use "array" for SunClock-CP or "WGBS" for SunClock
  coef_dir      = ".",          # Directory path where coefficient CSVs are located
  compute_accel = TRUE,         # Enable regression residual-based age acceleration
  pheno_df      = my_pheno      # Phenotype metadata containing 'Sample_ID' and 'Age'
)

# Preview the standardized results table
head(output_full)
