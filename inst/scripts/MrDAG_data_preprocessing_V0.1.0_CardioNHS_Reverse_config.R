# =============================================================================
# MrDAG_data_preprocessing_V0.1.0_CardioNHS_Reverse_config.R
#
# All user-provided parameters for the CardioNHS Reverse MrDAG data preprocessing run:
# paths, analysis settings, trait URLs/sample sizes, Exposures/Outcomes.
# Source this, then MrDAG_data_preprocessing.R, to run the analysis
# (see MrDAG_data_preprocessing_V0.1.0_CardioNHS_Reverse.R).
# =============================================================================


# --- Working directory ---
work_dir <- "/Users/leonardobottolo/Documents/Articles/C_Zuber/Article_7/Data/Real/R_Scripts"

# --- Directories ---
gwas_dir   <- "gwas_data/CardioNHS"         # Where raw GWAS files are stored
ref_dir    <- "ld_reference"                # Where the 1KG LD reference panel is stored
output_dir <- "Output/CardioNHS/Reverse"   # Where all output files are written

ldsc_munged_dir <- file.path(gwas_dir, "ldsc_munged")   # Default; user can override
ldsc_ref_dir    <- file.path(ref_dir, "ldsc_reference") # LD scores + HapMap3 SNP list (§1); default; user can override

# --- Analysis parameters ---
clump_r2  <- 0.01    # LD clumping r2 threshold
p_thresh  <- 5e-8    # Genome-wide significance threshold
clump_kb  <- 10000   # Clumping window (kb)
clump_p   <- 0.99    # Clumping p threshold

# --- Output files ---
r2_tag               <- as.character(clump_r2)
output_iv_rdata      <- file.path(output_dir, "IV.Rdata")
output_clumped_rdata <- file.path(output_dir, paste0("IV_clump_r2_", r2_tag, ".Rdata"))
output_matrix_rdata  <- file.path(output_dir, "MR_input.RData")
output_iv_summary_rdata <- file.path(output_dir, paste0("IV_summary_clump_r2_", r2_tag, ".RData"))
output_ldsc_gencov_rdata <- file.path(output_dir, "LDSC_output.RData")
output_mrdag_input_rdata <- file.path(output_dir, "MrDAG_input.RData")   # data + MrDAGcheck + gencov, ready for MrDAG()

# --- Global download timeout ---
options(timeout = 3600)

# --- Broad UKB annotation file (maps variant -> rsid, chr, pos, ref, alt) ---
annot_url  <- "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/annotations/variants.tsv.bgz"
annot_file <- file.path(gwas_dir, "variants.tsv.bgz")

# --- Dynamic EBI harmonised folder URL helper ---
get_ebi_url <- function(gcst_id) {
  num_id    <- as.numeric(gsub("GCST", "", gcst_id))
  lower     <- ((num_id - 1) %/% 1000) * 1000 + 1
  range_str <- paste0("GCST", sprintf("%06d", lower),
                      "-GCST", sprintf("%06d", lower + 999))
  paste0("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/",
         range_str, "/", gcst_id, "/harmonised/")
}

all_traits <- list(
  ApoA    = "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30630_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz",
  ApoB    = "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30640_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz",
  BMI     = "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz",
  HDL     = "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30760_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz",
  LDL     = "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30780_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz",
  SBP     = "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/4080_irnt.gwas.imputed_v3.both_sexes.tsv.bgz",
  TG      = "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30870_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz",
  AF      = "",    # Manual: place AF*.* in gwas_data/
  HF      = get_ebi_url("GCST009541"),
  CES     = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90104001-GCST90105000/GCST90104541/harmonised/GCST90104541.h.tsv.gz",
  PA      = "",    # Manual: place PA*.* in gwas_data/
  CAD     = "",    # Manual: place CAD*.* in gwas_data/
  PAD     = "",    # Manual: place PAD*.* in gwas_data/
  SMK     = "",    # Manual: place SMK*.* in gwas_data/
  T2D     = ""     # Manual: place T2D*.* in gwas_data/
)

trait_sample_info <- list(
  ApoA = list(type = "quant",  N = 317756, OLS = TRUE),
  ApoB = list(type = "quant",  N = 317756, OLS = TRUE),
  BMI  = list(type = "quant",  N = 359983, OLS = TRUE),
  HDL  = list(type = "quant",  N = 317788, OLS = TRUE),
  LDL  = list(type = "quant",  N = 317778, OLS = TRUE),
  SBP  = list(type = "quant",  N = 317754, OLS = TRUE),
  TG   = list(type = "quant",  N = 318024, OLS = TRUE),
  PA   = list(type = "quant",  N = 608595, OLS = TRUE),
  SMK  = list(type = "quant",  N = 462690, OLS = TRUE),

  T2D  = list(type = "binary", cases = 74124,  controls = 824006,
              population_prevalence = 0.105),
  AF   = list(type = "binary", cases = 60620,  controls = 970216,
              population_prevalence = 0.02),
  CAD  = list(type = "binary", cases = 71602,  controls = 260875,
              population_prevalence = 0.017),
  HF   = list(type = "binary", cases = 47309,  controls = 930014,
              population_prevalence = 0.017),
  PAD  = list(type = "binary", cases = 31307,  controls = 211753,
              population_prevalence = 0.038),
  CES  = list(type = "binary", cases = 10804,  controls = 1234808,
              population_prevalence = 0.003)
)

Exposures <- c("AF", "CAD", "CES", "HF", "PAD")
Outcomes  <- c("ApoA", "ApoB", "BMI", "HDL", "LDL", "PA", "SMK", "SBP", "T2D", "TG")

exposure_order <- NULL   # NULL = automatic (recommended)
outcome_order  <- NULL   # NULL = automatic (recommended)

# --- Startup validation ---
stopifnot(
  "Some Exposures are not in all_traits" = all(Exposures %in% names(all_traits)),
  "Some Outcomes are not in all_traits"  = all(Outcomes  %in% names(all_traits)),
  "Trait names must be unique"           = !anyDuplicated(names(all_traits))
)

manual_traits <- names(Filter(function(u) nchar(trimws(u)) == 0,
                              all_traits[names(all_traits) %in% c(Exposures, Outcomes)]))
if (length(manual_traits) > 0) {
  message("  [INFO] Manual traits declared: ", paste(manual_traits, collapse = ", "))
}

beta_X_cols <- paste0(tolower(Exposures), "_beta")
se_X_cols   <- paste0(tolower(Exposures), "_se")
se_Y_cols   <- paste0(tolower(Outcomes),  "_se")
beta_Y_cols <- paste0(tolower(Outcomes),  "_beta")
pval_X_cols <- paste0(tolower(Exposures), "_pval")

