################################################################################
# Supplementary Script: Epigenetic Age & Age Acceleration Prediction Pipeline
# 
# Description: This optimized script provides a user-friendly pipeline to calculate 
#              DNA methylation age (DNAm_age) and Age Acceleration using SunClock 
#              (WGBS) or SunClock-CP (Array) models. 
#
# RECOMMENDATION: 
#              This pipeline skips missing CpGs/NAs (treating them as 0 contributions).
#              To maximize accuracy, users are highly encouraged to perform baseline 
#              imputation (e.g., KNN) on their beta matrix prior to analysis.
################################################################################

#' Predict Epigenetic Age and Age Acceleration
#'
#' @param beta_matrix A numeric matrix or data frame of DNA methylation beta values 
#'                    (Format: rows = samples, columns = CpGs).
#' @param model_type A character string specifying the clock model type: "array" (Default, for SunClock-CP) or "WGBS" (for SunClock).
#' @param genome_version A character string specifying the genome build for WGBS model: "hg38" (Default) or "hg19". (Ignored for array model).
#' @param coef_dir A character string pointing to the folder directory where coefficient csv files are stored. Defaults to ".".
#' @param compute_accel Logical. If TRUE, calculates regression residual-based age acceleration.
#' @param pheno_df A data frame containing sample metadata. Required only if compute_accel = TRUE.
#'                 Must contain columns: "Sample_ID" (matching rownames of beta_matrix) and "Age".
#'
#' @return A data frame containing sample IDs, predicted DNAm_age, missing CpG counts, 
#'         and optionally, chronological age, calculated age acceleration, and all custom metadata columns.
predict_epigenetic_age <- function(beta_matrix, model_type = c("array", "WGBS"), 
                                   genome_version = c("hg38", "hg19"), coef_dir = ".",
                                   compute_accel = FALSE, pheno_df = NULL) {
  model_type <- match.arg(model_type)
  genome_version <- match.arg(genome_version)
  
  # 1. Parameter validation for Age Acceleration calculation
  if (compute_accel) {
    if (is.null(pheno_df)) stop("Error: 'pheno_df' must be provided when 'compute_accel = TRUE'.")
    if (!"Sample_ID" %in% colnames(pheno_df)) stop("Error: 'pheno_df' must contain a 'Sample_ID' column.")
    if (!"Age" %in% colnames(pheno_df)) stop("Error: 'pheno_df' must contain an 'Age' column.")
    if (is.null(rownames(beta_matrix))) stop("Error: 'beta_matrix' must have valid row names (Sample IDs).")
  }
  
  # 2. Automatically locate, verify, and read the matching coefficient file
  coef_filename <- if (model_type == "WGBS") "SunClock_coef.csv" else "SunClock-CP_coef.csv"
  coef_full_path <- file.path(coef_dir, coef_filename)
  
  if (!file.exists(coef_full_path)) {
    stop(paste0("Error: Coefficient file not found at [ ", coef_full_path, " ]. ",
                "Please verify the file location or adjust 'coef_dir'."))
  }
  coefs <- read.csv(coef_full_path, stringsAsFactors = FALSE)
  
  # 3. Dynamic filtration and extraction of Intercept and Coefficients
  if (model_type == "WGBS" && genome_version == "hg19") {
    # Verify that the hg19 coordinate column exists in the provided table
    if (!"CpG_pos_hg19" %in% colnames(coefs)) {
      stop("Error: 'CpG_pos_hg19' column not found in the SunClock coefficient table.")
    }
    # Filter out rows where the hg19 conversion failed (is NA or empty string)
    coefs <- coefs[!is.na(coefs$CpG_pos_hg19) & coefs$CpG_pos_hg19 != "", ]
  }
  
  # Extract intercept safely from the cleaned/filtered table
  intercept_row <- which(coefs$CpG_pos == "(Intercept)")
  if (length(intercept_row) == 0) stop("Error: Intercept row not found in the coefficient table.")
  
  intercept <- as.numeric(coefs[intercept_row, "Coefficient"])
  model_cpgs <- coefs[-intercept_row, ]
  
  # Select target ID column based on model configuration
  if (model_type == "array") {
    model_cpg_ids <- model_cpgs$CpG_name
  } else {
    model_cpg_ids <- if (genome_version == "hg38") model_cpgs$CpG_pos else model_cpgs$CpG_pos_hg19
  }
  
  # 4. Track missing features in the user's dataset
  data_cpg_ids <- colnames(beta_matrix)
  missing_cpgs_in_data <- setdiff(model_cpg_ids, data_cpg_ids)
  
  existing_model_cpgs <- intersect(model_cpg_ids, data_cpg_ids)
  na_counts <- if (length(existing_model_cpgs) > 0) {
    rowSums(is.na(beta_matrix[, existing_model_cpgs, drop = FALSE]))
  } else {
    rep(0, nrow(beta_matrix))
  }
  
  # 5. Initialize results data frame
  results <- data.frame(
    Sample_ID = rownames(beta_matrix),
    DNAm_age = numeric(nrow(beta_matrix)),
    Missing_CpG_Count = numeric(nrow(beta_matrix)),
    stringsAsFactors = FALSE
  )
  
  # 6. Perform Epigenetic Age calculation using Matrix Multiplication
  if (length(existing_model_cpgs) > 0) {
    beta_matched <- beta_matrix[, existing_model_cpgs, drop = FALSE]
    coef_indices <- match(existing_model_cpgs, model_cpg_ids)
    coefs_matched <- as.numeric(model_cpgs$Coefficient[coef_indices])
    
    # Set missing or NA values to 0 so they contribute nothing to the final score.
    beta_matched[is.na(beta_matched)] <- 0
    results$DNAm_age <- as.numeric(as.matrix(beta_matched) %*% coefs_matched + intercept)
  } else {
    warning("Warning: No overlapping CpG sites found between data and model coefficients.")
    results$DNAm_age <- intercept
  }
  results$Missing_CpG_Count <- length(missing_cpgs_in_data) + na_counts
  
  # 7. Optional: Compute Age Acceleration & Merge All Phenotype Columns via Sample_ID Column
  if (compute_accel) {
    match_idx <- match(results$Sample_ID, pheno_df$Sample_ID)
    if (any(is.na(match_idx))) stop("Error: Some Sample IDs in 'beta_matrix' could not be found in 'pheno_df$Sample_ID'.")
    
    pheno_matched <- pheno_df[match_idx, , drop = FALSE]
    chronological_age <- as.numeric(pheno_matched$Age)
    
    fit <- lm(DNAm_age ~ chronological_age, data = results)
    results$DNAm_age_accel <- as.numeric(residuals(fit))
    
    results_final <- data.frame(
      Sample_ID = results$Sample_ID,
      Age = chronological_age,
      DNAm_age = results$DNAm_age,
      DNAm_age_accel = results$DNAm_age_accel,
      Missing_CpG_Count = results$Missing_CpG_Count,
      stringsAsFactors = FALSE
    )
    
    other_metadata_cols <- setdiff(colnames(pheno_df), c("Age", "Sample_ID"))
    if (length(other_metadata_cols) > 0) {
      results_final <- cbind(results_final, pheno_matched[, other_metadata_cols, drop = FALSE])
    }
    results <- results_final
  }
  
  return(results)
}


################################################################################
# Quick Start & Usage Guide
################################################################################

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