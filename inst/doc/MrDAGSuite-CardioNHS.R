## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")


## ----chunk01, eval = FALSE----------------------------------------------------
# setwd("path/to/project/")   # project directory


## ----chunk02, eval = FALSE----------------------------------------------------
# # install.packages("pak")
# pak::pak("lb664/MrDAGSuite")
# library(MrDAGSuite)


## ----chunk03, eval = FALSE----------------------------------------------------
# # --- Working directory: the directory with the pipeline scripts and the data (by default the current one) ---
# work_dir   <- getwd()
# 
# # --- Directories ---
# gwas_dir   <- "gwas_data/CardioNHS"         # where the raw GWAS files are stored
# ref_dir    <- "ld_reference"                # where the 1000 Genomes LD reference panel is stored
# output_dir <- "Output/CardioNHS/Forward"    # where all output files are written
# 
# ldsc_munged_dir <- file.path(gwas_dir, "ldsc_munged")    # harmonised files used by LD Score Regression
# ldsc_ref_dir    <- file.path(ref_dir, "ldsc_reference")  # LD scores and HapMap3 SNP list
# 
# # --- Analysis parameters ---
# clump_r2  <- 0.01    # LD clumping r2 threshold
# p_thresh  <- 5e-8    # genome-wide significance threshold for an instrument
# clump_kb  <- 10000   # clumping window (kb)


## ----chunk04, eval = FALSE----------------------------------------------------
# all_traits <- list(
#   # Broad UKB / Neale Lab (downloaded automatically)
#   ApoA = "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30630_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz",
#   ApoB = "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30640_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz",
#   BMI  = "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/21001_irnt.gwas.imputed_v3.both_sexes.tsv.bgz",
#   HDL  = "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30760_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz",
#   LDL  = "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30780_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz",
#   SBP  = "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/4080_irnt.gwas.imputed_v3.both_sexes.tsv.bgz",
#   TG   = "https://broad-ukb-sumstats-us-east-1.s3.amazonaws.com/round2/additive-tsvs/30870_irnt.gwas.imputed_v3.both_sexes.varorder.tsv.bgz",
# 
#   # GWAS Catalog harmonised files (downloaded automatically)
#   HF   = get_ebi_url("GCST009541"),
#   CES  = "http://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90104001-GCST90105000/GCST90104541/harmonised/GCST90104541.h.tsv.gz",
# 
#   # Manual traits: place the file in gwas_dir, with a name starting with the trait name
#   AF   = "",
#   PA   = "",
#   CAD  = "",
#   PAD  = "",
#   SMK  = "",
#   T2D  = ""
# )


## ----chunk05, eval = FALSE----------------------------------------------------
# trait_sample_info <- list(
#   ApoA = list(type = "quant",  N = 317756, OLS = TRUE),
#   ApoB = list(type = "quant",  N = 317756, OLS = TRUE),
#   BMI  = list(type = "quant",  N = 359983, OLS = TRUE),
#   HDL  = list(type = "quant",  N = 317788, OLS = TRUE),
#   LDL  = list(type = "quant",  N = 317778, OLS = TRUE),
#   SBP  = list(type = "quant",  N = 317754, OLS = TRUE),
#   TG   = list(type = "quant",  N = 318024, OLS = TRUE),
#   PA   = list(type = "quant",  N = 608595, OLS = TRUE),
#   SMK  = list(type = "quant",  N = 462690, OLS = TRUE),
# 
#   T2D  = list(type = "binary", cases = 74124, controls = 824006, population_prevalence = 0.07),
#   AF   = list(type = "binary", cases = 60620, controls = 970216, population_prevalence = 0.024),
#   CAD  = list(type = "binary", cases = 71602, controls = 260875, population_prevalence = 0.037),
#   HF   = list(type = "binary", cases = 47309, controls = 930014, population_prevalence = 0.017),
#   PAD  = list(type = "binary", cases = 31307, controls = 211753, population_prevalence = 0.015),
#   CES  = list(type = "binary", cases = 10804, controls = 1234808, population_prevalence = 0.003)
# )


## ----chunk06, eval = FALSE----------------------------------------------------
# Exposures <- c("ApoA", "ApoB", "BMI", "HDL", "LDL", "PA", "SMK", "SBP", "T2D", "TG")
# Outcomes  <- c("AF", "CAD", "CES", "HF", "PAD")


## ----chunk04b, eval = FALSE---------------------------------------------------
# scripts <- system.file("scripts", package = "MrDAGSuite")
# file.copy(list.files(scripts, full.names = TRUE), "path/to/working/directory")   # any directory chosen by the user
# 
# setwd("path/to/working/directory")
# source("MrDAG_data_preprocessing_V0.1.0_CardioNHS_Forward.R")   # loads the configuration file and the pipeline


## ----chunk08, eval = FALSE----------------------------------------------------
# load(file.path(work_dir, output_mrdag_input_rdata))   # MrDAG_input.RData in the output folder
# ls()
# #> [1] "gencov"            "gencovDiagnostics" "MrDAG_data"
# #> [4] "MrDAGcheck"        "rs"
# 
# dim(MrDAG_data)               # instruments x traits
# #> [1] 1540   15
# colnames(MrDAG_data)          # outcomes first, then exposures
# #>  [1] "AF"   "CAD"  "CES"  "HF"   "PAD"  "ApoA" "ApoB" "BMI"  "HDL"  "LDL"
# #> [11] "PA"   "SMK"  "SBP"  "T2D"  "TG"
# MrDAGcheck                    # which columns are outcomes (Y_idx) and which are exposures (X_idx)
# #> $Y_idx
# #> [1] 1 2 3 4 5
# #> $X_idx
# #>  [1]  6  7  8  9 10 11 12 13 14 15
# length(gencov)                # one 15 x 15 sample-overlap matrix per instrument
# #> [1] 1540


## ----chunk09, eval = FALSE----------------------------------------------------
# pp_edge <- 2 / (ncol(MrDAG_data) - 1)    # edge prior, 0.1429 for 15 traits


## ----chunk10, eval = FALSE----------------------------------------------------
# output <- MrDAG(data = MrDAG_data, gencov = gencov,
#                 niter = 1200000, burnin = 200000, thin = 100,
#                 pp = pp_edge, tempMax = 10,
#                 MRconstraint = MrDAGcheck,
#                 fileName = "MrDAG_CardioNHS", filePath = file.path(work_dir, output_dir),
#                 seed = 28061971)


## ----chunk11, eval = FALSE----------------------------------------------------
# gencovCorrect <- get_gencovCorrect(MrDAG_data, gencov = gencov)
# gencovCorrect


## ----chunk12, eval = FALSE----------------------------------------------------
# ord <- c(c(2, 5, 4, 1, 3), c(6, 9, 7, 10, 11, 12, 8, 13, 14, 15))


## ----chunk13, eval = FALSE----------------------------------------------------
# PEIP <- get_postEdgeProb(output, ord = ord)
# round(PEIP$postEdgeProb, 2)


## ----chunk14, eval = FALSE----------------------------------------------------
# TPM <- get_TPM(output, ord = ord, threshold = 0.5)
# 
# MAP <- get_MAP(output, ord = ord, model = "CPDAG", sortby = "logMargLik")


## ----chunk15, eval = FALSE----------------------------------------------------
# # effect of the exposure in column 7 (ApoB) on the outcome in column 2 (CAD)
# effect <- get_genPredEffect(output, 2, 7)


## ----chunk16, eval = FALSE----------------------------------------------------
# effects <- get_genPredEffects(output, CI = 0.95, ord = ord, progress = TRUE)
# 
# round(effects$genPredEffects, 2)      # posterior mean
# round(effects$genPredEffects_LL, 2)   # lower limit of the interval
# round(effects$genPredEffects_UL, 2)   # upper limit of the interval


## ----chunk17, eval = FALSE----------------------------------------------------
# est <- effects$genPredEffects
# big <- which(abs(est) >= 0.05, arr.ind = TRUE)
# data.frame(exposure = rownames(est)[big[, 1]],
#            outcome  = colnames(est)[big[, 2]],
#            mean     = round(est[big], 2),
#            lower    = round(effects$genPredEffects_LL[big], 2),
#            upper    = round(effects$genPredEffects_UL[big], 2))


## ----chunk18, eval = FALSE----------------------------------------------------
# plot_convDiagnostic(output, fileName = "MrDAG_CardioNHS_convDiagnostic",
#                     filePath = "Figures/CardioNHS/")


## ----chunk19, eval = FALSE----------------------------------------------------
# # all PEIPs, in the order given by ord
# plot_postEdgeProb(output, ord = ord, fileName = "MrDAG_CardioNHS_postEdgeProb_I",
#               filePath = "Figures/CardioNHS/")
# 
# # only the edges with PEIP at least 0.5
# plot_postEdgeProb(output, ord = ord, threshold = 0.5, fileName = "MrDAG_CardioNHS_postEdgeProb_III",
#               filePath = "Figures/CardioNHS/")
# 
# # genetically predicted effects with 95% credible intervals, effects below 0.01 in absolute value left blank
# plot_genPredEffects(output, ord = ord, CI = 0.95, threshold = 0.01,
#                    fileName = "MrDAG_CardioNHS_genPredEffects_III",
#                    filePath = "Figures/CardioNHS/")


## ----chunk20, eval = FALSE----------------------------------------------------
# qgraph_postEdgeProb(output, threshold = 0.5, legend = TRUE,
#                 fileName = "MrDAG_CardioNHS_qgraph_postEdgeProb",
#                 filePath = "Figures/CardioNHS/")
# 
# qgraph_genPredEffects(output, threshold = 0.01, legend = TRUE,
#                      fileName = "MrDAG_CardioNHS_qgraph_genPredEffects",
#                      filePath = "Figures/CardioNHS/")


## ----chunk21, eval = FALSE----------------------------------------------------
# library(ggplot2)
# library(ggpubr)
# 
# fig_edges   <- plot_postEdgeProb(output, donotplot = TRUE, ord = ord,
#                              plotPar = list(cex.axis = 0.8, cex.lab = 0.8, geom_text.size = 1.8)) +
#                theme(legend.key.height = unit(0.25, "cm"))
# 
# fig_effects <- plot_genPredEffects(output, donotplot = TRUE, ord = ord, CI = 0.95, threshold = 0.01,
#                                   plotPar = list(cex.axis = 0.8, cex.lab = 0.8, geom_text.size = 2.2)) +
#                theme(legend.key.height = unit(0.25, "cm"))
# 
# pdf("Figures/CardioNHS/Fig_CardioNHS.pdf", width = 7.5, height = 4.8)
# print(ggarrange(ggarrange(fig_edges,   ncol = 1, legend = "bottom"),
#                 ggarrange(fig_effects, ncol = 1, legend = "bottom"),
#                 ncol = 2, widths = c(1.6, 1)))
# dev.off()


## ----chunk22, eval = FALSE----------------------------------------------------
# invIV_output <- invIV4MrDAG(data = MrDAG_data, output, gencov = gencov)
# invIV_output$invalidIV    # names of the invalid instruments


## ----chunk23, eval = FALSE----------------------------------------------------
# boot <- boot4MrDAG(data = MrDAG_data, gencov = gencov, nboot = 100,
#                    niter = 200000, burnin = 100000, thin = 100,
#                    pp = pp_edge, tempMax = 10, MRconstraint = MrDAGcheck,
#                    fileName = NULL, seed = 28061971)
# 
# plot_bootEdgeFreq(boot, ord = ord, fileName = "boot4MrDAG_CardioNHS_bootEdgeFreq",
#                   filePath = "Figures/CardioNHS/")
# plot_bootGenPredEffects(boot, ord = ord, CI = 0.95, threshold = 0.01,
#                        fileName = "boot4MrDAG_CardioNHS_bootGenPredEffects",
#                        filePath = "Figures/CardioNHS/")
# 
# analysis <- analyzer4MrDAG(data = MrDAG_data, gencov = gencov, ninvIV = 5, nboot = 100,
#                            niter = 200000, burnin = 100000, thin = 100,
#                            pp = pp_edge, tempMax = 10, MRconstraint = MrDAGcheck,
#                            fileName = NULL, seed = 28061971)


## ----chunk24, eval = FALSE----------------------------------------------------
# simData <- simMrDAG(nRep = 1, q = 5, p = 15, nSubject_X = 50000,
#                     ratioSubject_XY = 1,   # two-sample design, no overlap between the samples (default)
#                     nIV = 100,             # number of instruments (default)
#                     fileName = NULL, seed = 28061971)


## ----fig2, eval = TRUE, echo = FALSE, out.width = "75%", fig.align = "center"----
knitr::include_graphics("figures/Fig_2_CardioNHS_Forward.png")

