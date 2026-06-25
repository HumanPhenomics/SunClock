# **SunClock: an epigenetic age estimator trained in young individuals**

This repository provides the official model coefficients and an optimized R pipeline to calculate DNA methylation age (DNAm age) and Age Acceleration using the **SunClock** and **SunClock-CP** epigenetic clocks.

- **SunClock**: Developed using Whole-Genome Bisulfite Sequencing (WGBS) data.
- **SunClock-CP**: Tailored for traditional methylation microarrays (Illumina 450K/EPICv1/EPICv2/MSA).

---

## **System Requirements**

### **Software Dependencies**

- **R** (version ≥ 4.0.0)

### **Operating Systems**

The software has been tested and verified on:

- **Linux** (Ubuntu 20.04 LTS, 22.04 LTS)
- **macOS** (Monterey 12, Ventura 13)
- **Windows** (10, 11) via RStudio or R GU

---

## **Installation Guide**

### **Step 1: Clone or Download the Repository**

```bash
git clone https://github.com/HumanPhenomics/SunClock.git
cd SunClock
```

Alternatively, download the ZIP file from the GitHub repository page and extract it to your working directory.

### **Step 2: Verify Installation**

```R
# Source the main script to verify everything is working
source("SunClock_script.R")
# If no error messages appear, installation is successful
```

---

## **Instructions for Use**

### **1. Prepare Your Input Data**

**Methylation Beta Matrix (`beta_matrix`)**

- A numeric matrix or data frame of DNA methylation beta values (ranging from 0 to 1).
- **Rows**: Sample IDs (set as `rownames`).
- **Columns**: CpG site identifiers:
    - For `model_type = "array"`: Column names must match Illumina probe IDs (e.g., `cg00034076`). Compatible with Illumina 450K, EPICv1, EPICv2, and MSA platforms.
    - For `model_type = "WGBS"`: Column names must match genomic coordinates (e.g., `chr12:21771411:21771412`).

> ⚠️ **Important Note on Data Preprocessing**
> 
> 
> This pipeline handles missing or `NA` features by assigning them a value of `0` (meaning they contribute nothing to the cumulative score). To prevent severe underestimation of epigenetic age due to high data missingness, it is **strongly recommended** that users perform proper baseline imputation (e.g., KNN or median imputation) on their beta matrix prior to running this script.
> 

**Phenotype Metadata (`pheno_df`)**

*(Optional, only required if computing Age Acceleration)*

- A data frame **must** contain the following two columns:
    - `Sample_ID`: Character string matching the `rownames` of your `beta_matrix`.
    - `Age`: Numeric values representing the chronological age.
- *Note*: Any additional clinical or biological columns (e.g., `Sex`, `Group`, `Cell_Types`) will be automatically preserved and aligned in the final output.

### **2. Run the Prediction**

The core function `predict_epigenetic_age()` accepts the following parameters:

**Function Parameters**

| **Parameter** | **Type** | **Required/Optional** | **Description** |
| --- | --- | --- | --- |
| `beta_matrix` | Numeric matrix/data.frame | **Required** | DNA methylation beta values (0-1). **Rows** = samples, **Columns** = CpG probes. Row names must be sample IDs. |
| `model_type` | Character | Optional (default: `"array"`) | Clock model to use: •  `"array"` •  `"WGBS"` |
| `genome_version` | Character | Optional (default: `"hg38"`) | Genome build for WGBS model only. Options: `"hg38"` or `"hg19"`. Ignored for array model. |
| `coef_dir` | Character | Optional (default: `"."`) | Directory path containing coefficient CSV files (`SunClock_coef.csv` or `SunClock-CP_coef.csv`). |
| `compute_accel` | Logical | Optional (default: `FALSE`) | If `TRUE`, calculates age acceleration (residual from regressing DNAm age on chronological age). |
| `pheno_df` | Data.frame | **Required if `compute_accel = TRUE`** | Sample metadata. Must contain columns: • `Sample_ID` - matches `rownames(beta_matrix)` • `Age` - chronological age in years*Additional columns are preserved in output.* |

**Example 1: Basic DNAm Age Prediction (Array model, SunClock-CP)**

```R
source("SunClock_script.R")
final_output <- predict_epigenetic_age(
  beta_matrix = my_betas,
  model_type  = "array"
)
```

**Example 2: WGBS Prediction with hg19 Coordinates**

```R
final_output <- predict_epigenetic_age(
  beta_matrix    = my_betas,
  model_type     = "WGBS",
  genome_version = "hg19"
)
```

**Example 3: Full Pipeline with Age Acceleration and Metadata**

```R
final_output <- predict_epigenetic_age(
  beta_matrix    = my_betas,
  model_type     = "WGBS",
  genome_version = "hg38",
  coef_dir       = ".",
  compute_accel  = TRUE,
  pheno_df       = my_pheno
)
```

### **3. Interpreting Results**

The `final_output` data frame contains:

| **Column** | **Description** |
| --- | --- |
| `Sample_ID` | Sample identifier |
| `DNAmAge` | Predicted epigenetic age in years |
| `DNAm_age_accel` | Age acceleration (residual from regressing DNAmAge on chronological age) |
| `Missing_CpG_Count` | Missing clock CpGs for each sample |
| *(other columns)* | All additional columns from your `pheno_df` |

---

## **Demo**

### **Demo Data**

The repository includes a simulated methylation beta matrix (`demo_methylation.csv`) with simulating Illumina array data of 4 samples.

### **Run the Demo**

```R
# Source the pipeline function
source("SunClock_script.R")

# Load the demo methylation data
demo_methylation_data <- read.csv("demo_methylation.csv", row.names = 1)
head(demo_methylation_data)

# Prepare data: transpose so rows = samples, columns = CpG probes
my_betas <- t(demo_methylation_data)

# Prepare phenotype data (requires Sample_ID and Age columns)
my_pheno <- data.frame(
  Sample_ID = c("Sample1", "Sample2", "Sample3"),
  Age = c(35, 36, 37)
)

# Run SunClock-CP (array version) on demo data with age acceleration
final_output <- predict_epigenetic_age(
  beta_matrix    = my_betas,
  model_type     = "array",
  coef_dir       = ".",
  compute_accel  = TRUE,
  pheno_df       = my_pheno
)

# View results
print(final_output)
```

### **Expected Output**

```
  Sample_ID Age DNAm_age DNAm_age_accel Missing_CpG_Count
1   Sample1  35 34.26206      0.2993821                89
2   Sample2  36 36.23559     -0.5987642                89
3   Sample3  37 40.00541      0.2993821                89
```

### **Expected Run Time for Demo**

**Less than 10 seconds** on a standard desktop computer (8GB RAM, Intel Core i5 or equivalent).

---

## **Citation**

---

## **Contact**

For questions or issues, please open a GitHub issue or contact:

- Yuxin Sun: sunyuxin2022@sinh.ac.cn

---

## **License**

This project is licensed under the MIT License - see the LICENSE file for details.



