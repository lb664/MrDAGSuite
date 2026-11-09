# =============================================================================
# MrDAG_data_preprocessing.R
# Version: 0.1.0
#
# Not meant to be sourced directly: source the run's *_config.R first, which
# declares every variable this file expects (work_dir, gwas_dir, output_dir,
# all_traits, trait_sample_info, Exposures, Outcomes, etc.).
# =============================================================================

# =============================================================================
# § 1 — ENVIRONMENT & DEPENDENCIES
# =============================================================================

step_time <- proc.time()

setwd(work_dir)

## Resolve every path from *_config.R to absolute now, once, while cwd is
## guaranteed correct (= work_dir): everything else in this script -- and
## anything called from it, e.g. estimate_ldsc_gencov()'s munge_dir/
## hm3_snplist arguments -- then stays correct regardless of whether the
## working directory drifts later (e.g. if a later pipeline stage is
## re-invoked directly in a fresh session without re-running from the top of
## this script; this is what stranded munge()'s per-trait .sumstats.gz output
## in the wrong directory instead of ldsc_munged/ on a later PsychCog rerun).
## ldsc_munged_dir/ldsc_ref_dir are normalised explicitly, not just
## gwas_dir/ref_dir, since *_config.R already computed them via
## file.path(gwas_dir, ...)/file.path(ref_dir, ...) BEFORE this line runs --
## reassigning gwas_dir/ref_dir here would not retroactively fix them.
work_dir        <- normalizePath(work_dir,        mustWork = TRUE)
gwas_dir        <- normalizePath(gwas_dir,        mustWork = FALSE)
ref_dir         <- normalizePath(ref_dir,         mustWork = FALSE)
output_dir      <- normalizePath(output_dir,      mustWork = FALSE)
ldsc_munged_dir <- normalizePath(ldsc_munged_dir, mustWork = FALSE)
ldsc_ref_dir    <- normalizePath(ldsc_ref_dir,    mustWork = FALSE)

if (!dir.exists(output_dir))
  dir.create(output_dir, recursive = TRUE)

## Internally, trait-derived column names (beta_X_cols, se_Y_cols, etc.) are
## lowercase throughout this script -- left as-is below, since that
## convention is used in every column-detection/extraction step. Saved
## output uses the proper case declared in Exposures/Outcomes instead, for
## consistency across files; standardise_trait_names() does that renaming
## on a local copy right before each save, never touching the working
## lowercase objects the rest of the script relies on.
trait_case_map     <- setNames(c(Outcomes, Exposures), tolower(c(Outcomes, Exposures)))
trait_case_map_inv <- setNames(names(trait_case_map), trait_case_map)   # proper case -> lowercase

standardise_trait_names <- function(nm, case_map = trait_case_map) {
  for (lc in names(case_map)) {
    hit <- grepl(paste0("^", lc, "_"), nm)
    if (any(hit)) nm[hit] <- sub(paste0("^", lc, "_"), paste0(case_map[[lc]], "_"), nm[hit])
  }
  nm
}

# Same idea for bare trait-name vectors (no "_suffix"), e.g. data$exposures.
standardise_trait_values <- function(x, case_map = trait_case_map) {
  hit    <- x %in% names(case_map)
  x[hit] <- unname(case_map[x[hit]])
  x
}

## Rounds for on-screen display: round(x, digits) is more readable than signif() for typical-magnitude
## numbers, but falls back to signif(x, digits) wherever round() would show a small-but-real value as a
## misleading 0. Never affects anything saved to disk -- only console/message output.
.display_round <- function(x, digits = 3) {
  r <- round(x, digits)
  zeroed <- !is.na(x) & x != 0 & r == 0
  if (any(zeroed)) r[zeroed] <- signif(x[zeroed], digits)
  r
}

# --- Banner ---
message(
  "\n",
  "##\n",
  "## Bayesian Mendelian randomization with Directed Acyclic Graphs exploration\n",
  "## MrDAGSuite Data Preprocessing\n",
  "## Version: 0.1.0\n",
  "##\n",
  "## Copyright (C) 2026 L. Bottolo and V. Zuber\n",
  "##\n",
  "##"
)

# --- Manual file pre-flight checklist ----------------------------------------
if (length(manual_traits) > 0) {
  .chk_exts <- c(".tsv.gz", ".gz", ".bgz", ".tsv.bgz", ".tsv", ".txt", ".csv", ".tbl")

  message("")
  message("[MANUAL FILE CHECKLIST]")
  message("  Verify name, extension, and size for each manually placed file.")

  for (.tr in manual_traits) {
    .cands <- list.files(gwas_dir, pattern = paste0("^", .tr, "[^a-zA-Z0-9]"),
                         full.names = TRUE)
    .cands <- .cands[grepl(paste(gsub("\\.", "\\\\.", .chk_exts), collapse = "|"),
                           .cands, ignore.case = TRUE)]

    if (length(.cands) == 0L) {
      message(sprintf("  [MISS] %-8s  No file found in '%s/' starting with '%s'",
                      .tr, gwas_dir, .tr))
      message(sprintf("         Accepted extensions: %s",
                      paste(.chk_exts, collapse = "  ")))
      next
    }

    .src   <- .cands[which.max(file.info(.cands)$size)]
    .fname <- basename(.src)
    .sz    <- file.info(.src)$size
    .p     <- strsplit(.fname, "\\.")[[1L]]
    .n     <- length(.p)
    .ext   <- if (.n >= 3L && .p[.n] %in% c("gz", "bgz", "bz2", "xz", "zip"))
                 paste0(".", .p[.n - 1L], ".", .p[.n])
               else if (.n >= 2L) paste0(".", .p[.n])
               else ""

    .ok_size   <- .sz >= 1024^2
    .status    <- if (.ok_size) "OK" else "WARN"

    .clean_name <- sub(paste0("\\", .ext, "$"), "", .fname)
    if (nchar(.clean_name) > 45)
      .clean_name <- paste0(substr(.clean_name, 1, 42), "...")

    message(sprintf("  [%-4s] %-8s  File: %-45s  Ext: %-10s  Size: %6.1f MB",
                    .status, .tr, .clean_name, .ext, .sz / 1024^2))
    if (!.ok_size)
      message("         [WARNING] File < 1 MB — verify this is the correct file.")
  }

  rm(.chk_exts, .tr, .cands, .src, .fname, .sz, .p, .n, .ext,
     .ok_size, .status, .clean_name)

  # --- Column synonym reminder -------------------------------------------------
  message("")
  message("[REMINDER — COLUMN SYNONYM LIST (§4)]")
  message("  If any manually placed file uses non-standard column names, add them")
  message("  to the relevant synonym vector in §4. No other code changes are needed.")
  message("    .rsid_syns  .beta_syns  .se_syns  .pval_syns")
  message("    .ea_syns    .oa_syns    .eaf_syns")
}

message("\n", strrep("=", 70))
message("[Step 1] Verifying and installing dependencies...")
message(strrep("=", 70))

cran_packages <- c("data.table", "dplyr", "httr", "remotes", "R.utils",
                   "vroom", "tidyr")
for (pkg in cran_packages) {
  if (!suppressPackageStartupMessages(requireNamespace(pkg, quietly = TRUE))) {
    message("  [INSTALL] ", pkg); install.packages(pkg, quiet = TRUE)
  } else { message("  [OK] ", pkg) }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
for (pkg in c("VariantAnnotation")) {
  if (!suppressPackageStartupMessages(requireNamespace(pkg, quietly = TRUE))) {
    message("  [INSTALL] ", pkg)
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else { message("  [OK] ", pkg) }
}

github_packages <- list(ieugwasr          = "MRCIEU/ieugwasr",
                        genetics.binaRies = "MRCIEU/genetics.binaRies",
                        MVMR              = "WSpiller/MVMR",
                        GenomicSEM        = "GenomicSEM/GenomicSEM")
for (pkg_name in names(github_packages)) {
  # requireNamespace() (even with quietly = TRUE) loads the namespace to check it -- for
  # GenomicSEM that pulls in readr, which prints "Registered S3 methods overwritten by
  # 'readr'" the moment its namespace loads, regardless of quietly. Only wrapping the call
  # itself in suppressPackageStartupMessages() catches that.
  if (!suppressPackageStartupMessages(requireNamespace(pkg_name, quietly = TRUE))) {
    message("  [INSTALL] ", pkg_name)
    remotes::install_github(github_packages[[pkg_name]], quiet = TRUE)
  } else { message("  [OK] ", pkg_name) }
}

# GenomicSEM depends on OpenMx (CRAN), which is not installed via
# install_github() -- checked and installed separately here if needed.
if (!suppressPackageStartupMessages(requireNamespace("OpenMx", quietly = TRUE))) {
  message("  [INSTALL] OpenMx")
  install.packages("OpenMx", quiet = TRUE)
} else { message("  [OK] OpenMx") }

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(httr); library(remotes)
  library(vroom); library(tidyr); library(ieugwasr); library(genetics.binaRies)
  library(MVMR); library(GenomicSEM)
})

plink_path <- tryCatch(
  genetics.binaRies::get_plink_binary(),
  error = function(e) stop("FATAL: Could not locate PLINK binary. ", e$message)
)
message("  [OK] PLINK: ", plink_path)

ref_base              <- file.path(ref_dir, "1kg.v3")
ref_path_for_clumping <- file.path(ref_base, "EUR")
if (!dir.exists(ref_dir)) dir.create(ref_dir, recursive = TRUE)

if (!file.exists(paste0(ref_path_for_clumping, ".bed"))) {
  message("  [!] 1KG reference panel missing. Downloading (~2 GB)...")
  dest_tar  <- file.path(ref_dir, "1kg.v3.tgz")
  dl_result <- download.file("http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz",
                              destfile = dest_tar, mode = "wb")
  if (dl_result != 0) stop("FATAL: 1KG reference download failed.")
  untar(dest_tar, exdir = ref_base)
  if (!file.exists(paste0(ref_path_for_clumping, ".bed")))
    stop("FATAL: EUR.bed not found after extraction.")
  message("  [OK] 1KG reference panel installed.")
} else {
  message("  [OK] 1KG reference panel found.")
}

ld_dir      <- file.path(ldsc_ref_dir, "eur_w_ld_chr")
hm3_snplist <- file.path(ldsc_ref_dir, "w_hm3.snplist")
if (!dir.exists(ldsc_ref_dir)) dir.create(ldsc_ref_dir, recursive = TRUE)

if (!dir.exists(ld_dir) || length(list.files(ld_dir, pattern = "\\.l2\\.ldscore")) == 0) {
  message("  [!] LD score reference missing. Downloading (~33 MB)...")
  dest_ld <- file.path(ldsc_ref_dir, "eur_w_ld_chr.tar.gz")
  dl_ld   <- download.file("https://zenodo.org/records/8182036/files/eur_w_ld_chr.tar.gz?download=1",
                            destfile = dest_ld, mode = "wb")
  if (dl_ld != 0) stop("FATAL: LD score reference download failed.")
  md5_ld <- tools::md5sum(dest_ld)
  if (md5_ld != "e2f16343c4cfaa76caa7d0c03d26b489")
    stop("FATAL: eur_w_ld_chr.tar.gz MD5 mismatch (got ", md5_ld,
         ") -- download is corrupted or truncated. Delete and retry.")
  untar(dest_ld, exdir = ldsc_ref_dir)
  if (!dir.exists(ld_dir))
    stop("FATAL: eur_w_ld_chr/ not found after extraction.")
  message("  [OK] LD score reference installed.")
} else {
  message("  [OK] LD score reference found.")
}

if (!file.exists(hm3_snplist)) {
  message("  [!] HapMap3 SNP list missing. Downloading...")
  dest_hm3 <- file.path(ldsc_ref_dir, "w_hm3.snplist.gz")
  dl_hm3   <- download.file("https://zenodo.org/records/7773502/files/w_hm3.snplist.gz?download=1",
                             destfile = dest_hm3, mode = "wb")
  if (dl_hm3 != 0) stop("FATAL: HapMap3 SNP list download failed.")
  md5_hm3 <- tools::md5sum(dest_hm3)
  if (md5_hm3 != "153ecc2bcfa740afafe656e6a384d769")
    stop("FATAL: w_hm3.snplist.gz MD5 mismatch (got ", md5_hm3,
         ") -- download is corrupted or truncated. Delete and retry.")
  R.utils::gunzip(dest_hm3, hm3_snplist, remove = FALSE, overwrite = TRUE)
  if (!file.exists(hm3_snplist))
    stop("FATAL: w_hm3.snplist not found after extraction.")
  message("  [OK] HapMap3 SNP list installed.")
} else {
  message("  [OK] HapMap3 SNP list found.")
}

ldsc_script <- file.path(work_dir, "estimate_ldsc_gencov.R")
if (!file.exists(ldsc_script))
  stop("FATAL: estimate_ldsc_gencov.R not found in work_dir: ", ldsc_script)
source(ldsc_script)
message("  [OK] estimate_ldsc_gencov.R loaded.")
elapsed <- proc.time() - step_time
message(sprintf("  [DONE] Step 1 completed in %.1f seconds.", elapsed["elapsed"]))

# Synonym lists — add new column names here to support additional file formats
.beta_syns <- c("hm_beta","beta","Beta","logOR","LogOR","OR","Effect_A2",
                "Effect","BETA","b","effect")          # effect: EvangelouE
.se_syns   <- c("standard_error","se","SE","se_gc","StdErr","stderr",
                "StdErrLogOR","hm_se",
                "SE_SLEEPDURATION")                   # SLEEP.zip compound col
.pval_syns <- c("p_value","pval","Pvalue","pvalue","p-value_gc","P.value",
                "P","PVAL","Pval","hm_p_val","Pval",
                "P_SLEEPDURATION")                    # SLEEP.zip compound col
.ea_syns   <- c("hm_effect_allele","effect_allele","EA","EFFECT_ALLELE",
                "A1","ALLELE1","ALT","Allele1","minor_allele","allele1",
                "effectAllele")
.oa_syns   <- c("hm_other_allele","other_allele","noneffect_allele",
                "NEA","OTHER_ALLELE","A2","ALLELE0","REF","Allele2","allele2",
                "otherAllele","allele0")
.eaf_syns  <- c("hm_effect_allele_frequency","effect_allele_frequency",
                "EAF","Freq1","FRQ_A","minor_AF","A1FREQ","freq1",
                "a1freq","ImpMAF")
.rsid_syns <- c("hm_rsid","rs_dbSNP147","rs_dbsnp147","snptestid","rsid",
                "ID","SNP","variant_id","variant","snp","rsids")

.all_synonyms <- unique(c(
  .rsid_syns, .beta_syns, .se_syns, .pval_syns,
  .ea_syns,   .oa_syns,   .eaf_syns,
  # chrpos columns required by detect_columns() / build_chrpos()
  "hm_chrom","hm_pos",
  "chromosome","base_pair_location",
  "Chr","Pos","CHROM","POS","CHR","BP","chr","bp_hg19","Chrom",
  "MarkerName","markername","variant",
  # compound SLEEP.zip column (beta)
  "BETA_SLEEPDURATION",
  # needed when MarkerName contains rsids — CHR+POS used for chrpos instead
  "CHR","POS","Chr","Pos"
))

.pick <- function(syns, cn) {
  h <- syns[syns %in% cn]; if (length(h) > 0) h[1] else NA_character_
}

detect_columns <- function(cn, trait_name = "") {
  prefix <- if (nchar(trait_name) > 0) paste0("[", trait_name, "] ") else ""

  if      (all(c("hm_chrom","hm_pos") %in% cn))                                   chrpos_method <- "hm_cols"
  else if (all(c("chromosome","base_pair_location") %in% cn) && !"variant" %in% cn) chrpos_method <- "base_cols"
  else if (all(c("Chr","Pos") %in% cn))                                           chrpos_method <- "Chr_Pos"
  else if (all(c("Chrom","Pos") %in% cn))                                         chrpos_method <- "Chrom_Pos"
  else if (all(c("CHROM","POS") %in% cn))                                         chrpos_method <- "CHROM_POS"
  else if (all(c("CHR","BP") %in% cn))                                            chrpos_method <- "CHR_BP"
  else if (all(c("CHR","POS") %in% cn))                                           chrpos_method <- "CHR_POS"
  else if (all(c("chr","bp_hg19") %in% cn))                                       chrpos_method <- "chr_bp_hg19"
  else if ("variant" %in% cn)                                                     chrpos_method <- "variant_field"
  else if ("MarkerName" %in% cn)                                                  chrpos_method <- "MarkerName_col"
  else if ("markername" %in% cn)                                                  chrpos_method <- "markername_col"
  else                                                                            chrpos_method <- "rsid_only"

  # Pass 2
  beta_col <- .pick(.beta_syns, cn); se_col   <- .pick(.se_syns,   cn)
  pval_col <- .pick(.pval_syns, cn); ea_col   <- .pick(.ea_syns,   cn)
  oa_col   <- .pick(.oa_syns,   cn); eaf_col  <- .pick(.eaf_syns,  cn)
  rsid_col <- .pick(.rsid_syns, cn)

  fatal <- FALSE
  for (field in c("beta_col","se_col","pval_col","ea_col")) {
    val <- get(field)
    if (is.na(val)) {
      message(prefix, "[ERROR] Required field '", field, "' not detected."); fatal <- TRUE
    } else if (!val %in% cn) {
      message(prefix, "[ERROR] Column '", val, "' not in data."); fatal <- TRUE
    }
  }
  if (fatal) {
    message(prefix, "[ERROR] Available: ", paste(cn, collapse=", "))
    message(prefix, "[ERROR] Skipping this trait."); return(NULL)
  }

  message(prefix, "[COLUMNS] chrpos=", chrpos_method,
          " beta=", beta_col, " se=", se_col, " pval=", pval_col,
          " ea=", ea_col,
          " oa=",  if(is.na(oa_col))  "NA" else oa_col,
          " eaf=", if(is.na(eaf_col)) "NA" else eaf_col,
          " rsid=",if(is.na(rsid_col))"NA" else rsid_col)

  list(chrpos_method=chrpos_method, rsid_col=rsid_col,
       beta_col=beta_col, se_col=se_col, pval_col=pval_col,
       ea_col=ea_col, oa_col=oa_col, eaf_col=eaf_col)
}

build_chrpos <- function(dt, chrpos_method) {
  switch(chrpos_method,
    "hm_cols" = {
      chr <- as.character(dt$hm_chrom); pos <- as.character(dt$hm_pos)
      na_r <- is.na(chr)|is.na(pos)|chr==""|pos==""
      if (any(na_r) && all(c("chromosome","base_pair_location") %in% names(dt))) {
        # Only fill where raw columns are actually populated
        raw_chr <- as.character(dt$chromosome[na_r])
        raw_pos <- as.character(dt$base_pair_location[na_r])
        ok <- !is.na(raw_chr) & !is.na(raw_pos) & raw_chr != "" & raw_pos != ""
        chr[na_r][ok] <- raw_chr[ok]
        pos[na_r][ok] <- raw_pos[ok]
      }
      paste(chr, pos, sep=":")
    },
    "base_cols"     = paste(as.character(dt$chromosome),
                            as.character(dt$base_pair_location), sep=":"),
    "Chr_Pos"       = paste(as.character(dt$Chr), as.character(dt$Pos), sep=":"),
    "Chrom_Pos"     = paste(as.character(dt$Chrom), as.character(dt$Pos), sep=":"),
    "chr_bp_hg19"   = paste(as.character(dt$chr), as.character(dt$bp_hg19), sep=":"),
    "MarkerName_col"= gsub("[_/][^:]*$", "",
                          iconv(as.character(dt$MarkerName),
                                to = "UTF-8", sub = "")),
    # EvangelouE: lowercase markername, format 5:29439275 (no allele suffix)
    "markername_col"= gsub("[_/][^:]*$", "",
                          iconv(as.character(dt$markername),
                                to = "UTF-8", sub = "")),
    "CHROM_POS"     = paste(as.character(dt$CHROM), as.character(dt$POS), sep=":"),
    "CHR_BP"        = paste(as.character(dt$CHR),   as.character(dt$BP),  sep=":"),
    "CHR_POS"       = paste(as.character(dt$CHR),   as.character(dt$POS), sep=":"),
    "variant_field" = gsub("^([^:]+:[^:]+):.*","\\1", dt$variant),
    "rsid_only"     = as.character(dt$SNP),
    stop("Unknown chrpos_method: ", chrpos_method)
  )
}

extract_rsid <- function(dt, rsid_col) {
  if (is.na(rsid_col) || !rsid_col %in% names(dt))
    return(rep(NA_character_, nrow(dt)))
  # iconv to UTF-8 to prevent regex failure on invalid bytes
  raw  <- iconv(as.character(dt[[rsid_col]]), to="UTF-8", sub="")
  rsid <- sub(".*?(rs[0-9]+).*","\\1", raw, perl=TRUE)
  rsid[!grepl("^rs[0-9]+$", rsid)] <- NA_character_
  rsid
}

get_trait_path <- function(trait) {
  if (!is.null(manual_paths[[trait]]))
    return(manual_paths[[trait]])
  if (!is.null(auto_paths[[trait]]) && file.exists(auto_paths[[trait]]))
    return(auto_paths[[trait]])
  # Pattern fallback: find largest file starting with TraitAcronym_
  cands <- list.files(gwas_dir, pattern = paste0("^", trait, "_"),
                       full.names = TRUE)
  cands <- cands[file.info(cands)$size >= 1024^2]
  if (length(cands) > 0)
    return(cands[which.max(file.info(cands)$size)])
  file.path(gwas_dir, paste0(trait, ".tsv.gz"))  # legacy fallback
}

count_rows_plain <- function(path) {
  # Fast row count for uncompressed files only. Returns (lines - 1) to
  # exclude the header, or NA if the command fails.
  tryCatch({
    raw <- system(paste0("wc -l < '", path, "'"),
                  intern = TRUE, ignore.stderr = TRUE)
    n   <- suppressWarnings(as.integer(trimws(raw[1L])))
    if (!is.na(n) && n >= 1L) n - 1L else NA_integer_
  }, error   = function(e) NA_integer_,
     warning = function(w) NA_integer_)
}

sort_traits <- function(traits) {
  classify <- function(t) {
    u <- trimws(all_traits[[t]])
    if (nchar(u) == 0)                          return(3L)  # manual
    if (grepl("broad-ukb|nealelab", u))          return(1L)  # UKB
    return(2L)                                               # GWAS Catalog
  }
  get_nrows <- function(t) {
    p <- get_trait_path(t)
    if (!file.exists(p)) return(0L)
    is_compressed <- grepl("\\.(gz|bgz)$", p, ignore.case = TRUE)
    if (is_compressed) {
      # File-size heuristic for compressed files — no decompression
      return(as.integer(file.info(p)$size / 100L))
    }
    # Uncompressed manually placed file — exact line count
    n <- count_rows_plain(p)
    if (is.na(n)) as.integer(file.info(p)$size / 100L) else n
  }
  groups <- sapply(traits, classify)
  counts <- sapply(traits, get_nrows)
  traits[order(groups, -counts)]
}

smart_fread <- function(path, trait_name = "") {
  prefix        <- if (nchar(trait_name) > 0) paste0("[", trait_name, "] ") else ""
  col_spec      <- vroom::cols(MarkerName = vroom::col_character(),
                               .default   = vroom::col_guess())
  is_compressed <- grepl("\\.(gz|bgz)$", path, ignore.case = TRUE)
  if (!is_compressed) {
    # Pass 1: header only
    hdr <- tryCatch(
      suppressMessages(suppressWarnings(
        vroom::vroom(path, n_max = 0L,
                     show_col_types = FALSE, progress = TRUE,
                     comment = "##")
      )),
      error = function(e) {
        message(prefix, "[ERROR] Cannot read file header: ", e$message)
        NULL
      }
    )
    if (is.null(hdr)) return(NULL)

    # Synonym match
    keep     <- intersect(.all_synonyms, names(hdr))
    keep_idx <- match(keep, names(hdr))   # integer positions — faster skip

    # Guard
    if (length(keep) < 4L) {
      message(prefix, "[ERROR] Fewer than 4 recognised columns found in header.")
      message(prefix, "[ERROR] Columns in file : ", paste(names(hdr), collapse = ", "))
      message(prefix, "[ERROR] Recognised       : ",
              if (length(keep) == 0L) "(none)" else paste(keep, collapse = ", "))
      message(prefix, "[ERROR] Add missing column names to the synonym lists in §4.")
      return(NULL)
    }

    # Pass 2 — fast path: integer positional col_select
    result <- tryCatch(
      as.data.table(suppressWarnings(
        vroom::vroom(path, col_select = keep_idx,
                     show_col_types = FALSE, progress = TRUE,
                     comment = "##")
      )),
      error = function(e) NULL
    )
    # Restore column names (positional selection drops them to V1, V2, ...)
    if (!is.null(result) && nrow(result) > 0L) {
      setnames(result, keep)
      message(prefix, "[READ] positional-select  ",
              length(keep), "/", ncol(hdr), " cols  \u2192 ", basename(path))
      return(result)
    }

    # Pass 2 — fallback: named col_select
    result <- tryCatch(
      as.data.table(suppressWarnings(
        vroom::vroom(path, col_select = dplyr::all_of(keep),
                     show_col_types = FALSE, progress = TRUE,
                     comment = "##")
      )),
      error = function(e) NULL
    )
    if (!is.null(result) && nrow(result) > 0L) {
      message(prefix, "[READ] named-select       ",
              length(keep), "/", ncol(hdr), " cols  \u2192 ", basename(path))
      return(result)
    }

    # Pass 2 — last resort: Latin-1 encoding
    message(prefix, "[READ] latin1-select      trying encoding fallback...")
    return(tryCatch(
      as.data.table(suppressWarnings(
        vroom::vroom(path, col_select = dplyr::all_of(keep),
                     show_col_types = FALSE, progress = TRUE,
                     comment = "##",
                     locale = vroom::locale(encoding = "latin1"))
      )),
      error = function(e) {
        message(prefix, "[ERROR] Cannot read file: ", e$message); NULL
      }
    ))
  }

  # Compressed files (.gz, .bgz, .tsv.gz): fread() is fastest.
  result <- tryCatch(
    fread(path, fill = TRUE, quote = "",
          nThread = max(1L, data.table::getDTthreads() - 2L)),
    error = function(e) NULL
  )
  if (!is.null(result) && nrow(result) > 0 &&
      !all(grepl("^V[0-9]+$", names(result)))) {
    if ("MarkerName" %in% names(result))
      result[, MarkerName := as.character(MarkerName)]
    return(result)
  }
  result <- tryCatch(
    as.data.table(suppressWarnings(
      vroom::vroom(path, show_col_types = FALSE, progress = TRUE,
                   col_types = col_spec, comment = "##")
    )), error = function(e) NULL
  )
  if (!is.null(result) && nrow(result) > 0) return(result)
  message(prefix, "[RETRY] Trying vroom() with Latin-1 encoding...")
  tryCatch(
    as.data.table(suppressWarnings(
      vroom::vroom(path, show_col_types = FALSE, progress = TRUE,
                   col_types = col_spec,
                   locale    = vroom::locale(encoding = "latin1"))
    )), error = function(e) {
      message(prefix, "[ERROR] Cannot read file: ", e$message); NULL
    })
}

to_num_matrix <- function(df, cols) {
  m <- vapply(df[, cols, drop=FALSE],
              function(x) suppressWarnings(as.numeric(x)),
              numeric(nrow(df)))
  matrix(m, nrow=nrow(df), ncol=length(cols),
         dimnames=list(NULL, cols))
}

harmonise_allele <- function(df, trait, use_toupper = FALSE, remove_unmatched = TRUE) {
  tl <- tolower(trait)

  # Guard: EA and alt must exist — use [[ to avoid $ partial-matching NULL returns
  if (!"EA" %in% names(df))
    stop("[HARMONISE] EA column missing for trait ", trait,
         ". Columns: ", paste(names(df), collapse = ", "))
  if (!"alt" %in% names(df))
    stop("[HARMONISE] alt column missing for trait ", trait)

  ea  <- toupper(as.character(df[["EA"]]))
  nea <- if (!"NEA" %in% names(df)) NULL
         else toupper(as.character(df[["NEA"]]))
  alt <- toupper(as.character(df[["alt"]]))

  n_match_ea  <- sum(alt == ea,  na.rm = TRUE)
  n_match_nea <- if (is.null(nea)) 0L else sum(alt == nea, na.rm = TRUE)
  n_odd       <- if (is.null(nea)) 0L
                 else sum(alt != ea & alt != nea, na.rm = TRUE)

  message("  [HARMONISE] ", trait,
          ": EA=", n_match_ea, " NEA=", n_match_nea, " Unmatched=", n_odd)

  if (n_odd > 0 && remove_unmatched) {
    keep <- is.na(alt) | alt == ea | (!is.null(nea) & alt == nea)
    df   <- df[keep, , drop = FALSE]
    message("  [DROP] Removing ", n_odd,
            " variants where neither EA nor NEA matches alt.")
    ea  <- toupper(as.character(df[["EA"]]))
    nea <- if (is.null(nea)) NULL
           else toupper(as.character(df[["NEA"]]))
    alt <- toupper(as.character(df[["alt"]]))
  }

  beta_col_raw <- paste0(tl, "_beta_not_aligned")
  beta_raw <- suppressWarnings(as.numeric(df[[beta_col_raw]]))

  df[[paste0(tl, "_beta")]] <- ifelse(alt == ea, beta_raw, -beta_raw)

  drop_cols <- intersect(c("EA", "NEA", beta_col_raw), names(df))
  df <- df[, !names(df) %in% drop_cols, drop = FALSE]
  df
}

manual_paths <- list()
auto_paths   <- list()

# Resolves a trait's actual readable file path -- extracting from a ZIP
# archive to a temp file first if needed. Previously this logic lived only
# inside read_trait_lazy() (used for column detection); estimate_ldsc_gencov.R
# called vroom directly on trait_paths[[tr]] without it, so for any trait
# whose source file is actually a ZIP archive (confirmed directly: SP's
# file is a real .zip despite a misleadingly plain .txt-style basename).
# it was reading raw zip bytes as if they were the delimited text table --
# explaining "Elements CHR, BP, ALLELE1... don't exist" (vroom parsed
# garbage from compressed bytes, not the real header). Both column
# detection and the actual LDSC reads now resolve through this same
# function, so they can never disagree about which path to read.
resolve_trait_path <- function(path, trait_name = "") {
  prefix <- if (nchar(trait_name) > 0) paste0("[", trait_name, "] ") else ""
  if (!grepl("\\.zip$", path, ignore.case = TRUE)) return(path)
  tmp_dir  <- tempdir()
  zip_list <- utils::unzip(path, list = TRUE)
  zip_list <- zip_list[!grepl("README|readme", zip_list$Name), ]
  target   <- zip_list$Name[which.max(zip_list$Length)]
  utils::unzip(path, files = target, exdir = tmp_dir, overwrite = TRUE)
  resolved <- file.path(tmp_dir, target)
  message(prefix, "[ZIP] Extracted: ", basename(target))
  resolved
}

# vroom's own delimiter auto-guess is unreliable for some GWAS files: it
# fails outright for some (explicit "Could not guess the delimiter" errors)
# and silently mis-splits others (the whole header read as one garbled
# column name, rather than erroring). Confirmed directly: MDD's file is
# genuinely space-delimited ("MarkerName A1 A2 Freq LogOR StdErrLogOR P"),
# while AN/BD/SCZ use runs of MULTIPLE spaces for column alignment (e.g.
# "CHROM     POS     ID      REF     ALT BETA SE PVAL"), which a literal
# single-space delim would parse as many empty columns rather than the
# real ones -- there is no single delimiter safe to assume across all
# traits, and "whitespace" itself isn't one delimiter but two different
# formats. Returns "\t", "," or "WHITESPACE" (one-or-more run of spaces,
# to be read via readr::read_table() rather than vroom's fixed-delim mode).
detect_file_delim <- function(path) {
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con), add = TRUE)
  # Skip leading '##' metadata lines (PGC's VCF-style sumstats format, among
  # others) -- confirmed directly: AN/BD/SCZ's first line is literally
  # "##fileFormat=PGCsumstatsVCFv1.0", which has no delimiter at all, so
  # sniffing it directly always fell through to the "space" default. These
  # files are actually plain tab-delimited; the "multi-space" look in past
  # error dumps was tabs rendered as wide whitespace in the console, not a
  # real whitespace-padded format.
  first_line <- ""
  repeat {
    ln <- tryCatch(readLines(con, n = 1), error = function(e) character(0))
    if (length(ln) == 0) break
    if (!startsWith(ln, "##")) { first_line <- ln; break }
  }
  if (grepl("\t", first_line, fixed = TRUE)) return("\t")
  if (grepl(" {2,}", first_line)) return("WHITESPACE")
  if (grepl(",", first_line, fixed = TRUE)) return(",")
  " "
}

# Reads a delimited GWAS file, transparently switching to readr::read_table()
# (which collapses runs of whitespace into one delimiter) when detect_file_delim()
# identified column-aligned whitespace rather than a single consistent
# character -- vroom::vroom(delim = " ") would otherwise produce many empty
# fields between padded columns instead of the real ones.
read_gwas_delim <- function(path, delim, ...) {
  if (identical(delim, "WHITESPACE")) {
    if (!requireNamespace("readr", quietly = TRUE))
      stop("FATAL: readr is required to read whitespace-column-aligned ",
           "files (install.packages('readr')).")
    readr::read_table(path, ...)
  } else {
    vroom::vroom(path, delim = delim, ...)
  }
}

read_trait_lazy <- function(path, trait_name = "") {
  prefix <- if (nchar(trait_name) > 0) paste0("[", trait_name, "] ") else ""
  path   <- resolve_trait_path(path, trait_name)

  # ── vroom lazy load — force MarkerName/markername to character ────────────
  # Prevents chr:pos strings like "5:85928892" being parsed as HH:MM:SS times
  file_delim <- detect_file_delim(path)
  dt <- tryCatch(
    suppressMessages(suppressWarnings(
      read_gwas_delim(path, file_delim, show_col_types = FALSE, progress = TRUE,
                   comment = "##",
                   col_types = vroom::cols(
                     MarkerName = vroom::col_character(),
                     markername = vroom::col_character(),
                     .default   = vroom::col_guess()
                   ))
    )),
    error = function(e) {
      message(prefix, "[ERROR] Cannot read file: ", e$message); NULL
    }
  )
  if (is.null(dt) || nrow(dt) == 0L) return(NULL)

  # ── Strip leading '#' from column names (VCF #CHROM header) ──────────────
  names(dt) <- sub("^#", "", names(dt))

  # ── hm_beta all-NA fallback: use raw beta column instead ─────────────────
  if ("hm_beta" %in% names(dt) && "beta" %in% names(dt)) {
    hm_vals <- suppressWarnings(as.numeric(dt[["hm_beta"]]))
    if (all(is.na(hm_vals))) {
      message(prefix, "[NOTE] hm_beta is all NA — using raw beta column instead.")
      dt[["hm_beta"]] <- dt[["beta"]]
    }
  }

  # ── SLEEP.zip compound column renaming ────────────────────────────────────
  if ("BETA_SLEEPDURATION" %in% names(dt) && !"beta" %in% names(dt))
    names(dt)[names(dt) == "BETA_SLEEPDURATION"] <- "beta"
  if ("SE_SLEEPDURATION"   %in% names(dt) && !"SE"   %in% names(dt))
    names(dt)[names(dt) == "SE_SLEEPDURATION"]   <- "SE"
  if ("P_SLEEPDURATION"    %in% names(dt) && !"P"    %in% names(dt))
    names(dt)[names(dt) == "P_SLEEPDURATION"]    <- "P"

  # ── REGENIE log10p -> P (p = 10^(-log10p)) ──────────────────────────────
  if (!any(.pval_syns %in% names(dt))) {
    log10p_col <- intersect(c("log10p", "LOG10P"), names(dt))
    if (length(log10p_col) > 0) {
      dt[["P"]] <- 10^(-as.numeric(dt[[log10p_col[1]]]))
      message(prefix, "[NOTE] Converted '", log10p_col[1], "' -> 'P' (p = 10^(-log10p)).")
    }
  }

  # ── FRQ_A_* prefix → EAF (e.g. ADHD: FRQ_A_19099) ───────────────────────
  if (!any(.eaf_syns %in% names(dt))) {
    frq_a_col <- grep("^FRQ_A_", names(dt), value = TRUE)
    if (length(frq_a_col) > 0) {
      names(dt)[names(dt) == frq_a_col[1]] <- "EAF"
      message(prefix, "[NOTE] Renamed '", frq_a_col[1], "' -> 'EAF' (FRQ_A_* prefix match).")
    }
  }

  # ── MarkerName-contains-rsid detection ───────────────────────────────────
  if ("MarkerName" %in% names(dt)) {
    sample_vals <- head(as.character(dt[["MarkerName"]]), 200L)
    sample_vals <- sample_vals[!is.na(sample_vals)]
    pct_rs      <- mean(grepl("^rs", sample_vals))
    if (pct_rs > 0.8) {
      has_pos <- any(c("CHR","POS","Chr","Pos","CHROM") %in% names(dt))
      if (has_pos) {
        message(prefix, "[NOTE] MarkerName is ", round(pct_rs * 100),
                "% rs-prefixed and CHR/POS present — rerouting via CHR+POS.")
        names(dt)[names(dt) == "MarkerName"] <- "rsid_from_col"
        if (!"SNP" %in% names(dt))
          names(dt)[names(dt) == "rsid_from_col"] <- "SNP"
      } else {
        message(prefix, "[NOTE] MarkerName is ", round(pct_rs * 100),
                "% rs-prefixed and no CHR/POS columns found",
                " — renaming to SNP for rsid-based join.")
        if (!"SNP" %in% names(dt))
          names(dt)[names(dt) == "MarkerName"] <- "SNP"
      }
    }
  }

  # ── Lazy select: only retain synonym-matched columns ─────────────────────
  keep <- intersect(.all_synonyms, names(dt))
  if (length(keep) < 4L) {
    message(prefix, "[ERROR] Fewer than 4 recognised columns: ",
            paste(names(dt), collapse = ", "))
    return(NULL)
  }
  dt <- dt[, keep, drop = FALSE]
  message(prefix, "[READ] lazy-vroom  ",
          length(keep), " cols  \u2192 ", basename(path))
  dt
}

step_time <- proc.time()
message("\n", strrep("=", 70))
message("[Checkpoint] Checking whether clumping has already been completed...")
message(strrep("=", 70))

checkpoint_valid <- FALSE
if (file.exists(output_iv_rdata) && file.exists(output_clumped_rdata)) {
  checkpoint_load_ok <- tryCatch({
    load(output_iv_rdata)       # restores `IV`
    load(output_clumped_rdata)  # restores `IV_clump`
    TRUE
  }, error = function(e) {
    message("  [CHECKPOINT] Failed to load existing IV/IV_clump: ", e$message)
    FALSE
  })

  if (checkpoint_load_ok) {
    # Saved IV/IV_clump have proper-case trait columns (standardised at save
    # time, below); translate back to the lowercase this script uses
    # internally before the checkpoint check and everything downstream.
    names(IV)       <- standardise_trait_names(names(IV),       trait_case_map_inv)
    names(IV_clump) <- standardise_trait_names(names(IV_clump), trait_case_map_inv)
    expected_cols_clump <- c(beta_X_cols, se_X_cols, beta_Y_cols, se_Y_cols)
    missing_cols  <- c(setdiff(expected_cols_clump, names(IV_clump)),
                       setdiff(pval_X_cols,          names(IV)))
    if (length(missing_cols) > 0) {
      message("  [CHECKPOINT] Existing IV/IV_clump is missing columns expected ",
              "by the current Exposures/Outcomes configuration: ",
              paste(missing_cols, collapse = ", "),
              " -- treating checkpoint as invalid, running full computation.")
    } else {
      checkpoint_valid <- TRUE
    }
  }
}

if (checkpoint_valid) {
  message("  [OK] Valid IV_clump found (", format(nrow(IV_clump), big.mark = ","),
          " rows) -- skipping download, census, exposure/outcome merging, ",
          "instrument selection, and LD clumping entirely.")
  message("  [OK] Proceeding directly to F-statistic/LDSC (§8) using the ",
          "existing IV_clump.")

  IV <- IV %>% dplyr::rename(pval = min.pval)
} else {
  message("  [OK] No valid checkpoint found -- running full computation ",
          "from §2 (download) onward.")
}

if (!checkpoint_valid) {

# =============================================================================
# § 2 — URL RESOLUTION & DOWNLOAD
# =============================================================================

step_time <- proc.time()
message("\n", strrep("=", 70))
message("[Step 2] Resolving URLs and downloading GWAS files...")
message(strrep("=", 70))

if (!dir.exists(gwas_dir)) dir.create(gwas_dir, recursive = TRUE)

manual_exts <- c(".tsv.gz", ".gz", ".bgz", ".tsv.bgz", ".tsv", ".txt", ".csv", ".tbl")

check_size <- function(path, label) {
  if (!file.exists(path)) {
    message("  [ERROR] ", label, " — file not found."); return(FALSE)
  }
  sz <- file.info(path)$size
  if (sz < 1024^2) {
    message("  [WARNING] ", label, " — small (", round(sz/1024), " KB).")
    return(FALSE)
  }
  message("  [OK] ", label, " (", round(sz/1024^2, 1), " MB).")
  return(TRUE)
}

recompress_to_gz <- function(src, dest, label) {
  message("  [CONVERT] ", label, ": ", basename(src), " → ", basename(dest))
  # If already .gz, copy directly — do not re-wrap
  if (grepl("\\.gz$", src, ignore.case = TRUE)) {
    message("  [CONVERT] Source is already .gz — copying directly.")
    ok <- file.copy(src, dest, overwrite = TRUE)
    if (!ok) message("  [ERROR] file.copy() failed.")
    return(ok)
  }
  tryCatch({
    dt <- suppressWarnings(
      vroom::vroom(src, show_col_types = FALSE, progress = TRUE,
                   col_types = vroom::cols(.default = vroom::col_character()))
    )
    vroom::vroom_write(dt, dest, delim = "\t", progress = TRUE)
    message("  [CONVERT] Written: ",
            format(nrow(dt), big.mark = ","), " rows, ",
            ncol(dt), " columns.")
    TRUE
  }, error = function(e) {
    message("  [ERROR] Conversion failed: ", e$message)
    if (file.exists(dest)) unlink(dest); FALSE
  })
}

safe_content_length <- function(response) {
  if (is.null(response)) return(NA_real_)
  val <- httr::headers(response)[["content-length"]]
  if (is.null(val) || length(val) == 0L || is.na(val)) return(NA_real_)
  as.numeric(val)
}

zip_sourced   <- character(0)
active_traits <- c(Exposures, Outcomes)
download_list <- all_traits[names(all_traits) %in% active_traits]
download_list[["_ANNOT"]] <- annot_url

for (n in names(download_list)) {
  url <- download_list[[n]]

  # ── Annotation file: fixed destination, handled separately ──────────────────
  if (n == "_ANNOT") {
    dest <- annot_file
    if (file.exists(dest)) {
      local_size  <- file.info(dest)$size
      response    <- tryCatch(httr::HEAD(url), error = function(e) NULL)
      server_size <- safe_content_length(response)
      size_ok     <- local_size >= 1024^2 &&
                     (is.na(server_size) ||
                      abs(local_size - server_size) <= server_size * 0.01)
      if (size_ok) {
        message("  [SKIP] _ANNOT (", round(local_size/1024^2, 1), " MB)."); next
      }
      message("  [!] _ANNOT incomplete. Re-downloading..."); unlink(dest)
    }
    message("  [FETCH] Downloading annotation file...")
    dl_result <- tryCatch(download.file(url, dest, method = "libcurl", quiet = TRUE),
                          error = function(e) { message("  [ERROR] ", e$message); -1L })
    check_size(dest, "_ANNOT")
    next
  }

  # ── Manual trait — read directly from original file, no conversion ──────────
  if (nchar(trimws(url)) == 0) {
    candidates <- list.files(gwas_dir, pattern = paste0("^", n, "[^a-zA-Z0-9]"),
                              full.names = TRUE)
    candidates <- candidates[
      grepl(paste(gsub("\\.", "\\\\.", manual_exts), collapse = "|"),
            candidates, ignore.case = TRUE)
    ]
    if (length(candidates) == 0) {
      message("  [MANUAL] ", n, " — no file found starting with '", n,
              "' in '", gwas_dir, "/'.")
      next
    }
    src <- candidates[which.max(file.info(candidates)$size)]
    manual_paths[[n]] <- src
    message("  [MANUAL] ", n, " — ", basename(src),
            " (", round(file.info(src)$size/1024^2, 1), " MB).")
    next
  }

  # ── Resolve EBI harmonised folder ───────────────────────────────────────────
  if (grepl("/harmonised/$", url)) {
    message("  [RESOLVE] ", n)
    page <- tryCatch(httr::content(httr::GET(url), "text", encoding = "UTF-8"),
                     error = function(e) { message("  [ERROR] ", e$message); NULL })
    if (is.null(page)) next
    fname <- regmatches(page, regexpr("[^\"'>]+\\.h(\\.tsv)?\\.gz", page))
    if (length(fname) == 0 || is.na(fname[1]))
      fname <- regmatches(page, regexpr("[^\"'>]+\\.h\\.bgz", page))
    if (length(fname) == 0 || is.na(fname[1])) {
      message("  [ERROR] Cannot resolve filename for ", n); next
    }
    url <- paste0(url, fname[1])
    message("  [RESOLVED] -> ", url)
  }

  # ── ZIP download ─────────────────────────────────────────────────────────────
  if (grepl("\\.zip$", url, ignore.case = TRUE)) {
    zip_dest    <- file.path(gwas_dir, paste0(n, "_download.zip"))
    extract_dir <- file.path(gwas_dir, paste0(n, "_zip_extract"))
    # Pre-check: if a correctly named output file already exists, skip the download
    prev_out <- list.files(gwas_dir, pattern = paste0("^", n, "_"), full.names = TRUE)
    prev_out <- prev_out[!grepl("_download\\.zip$|_zip_extract", prev_out) &
                          file.info(prev_out)$size >= 1024^2]
    if (length(prev_out) > 0) {
      best <- prev_out[which.max(file.info(prev_out)$size)]
      message("  [SKIP] ", n, " (ZIP-sourced, ", basename(best), ", ",
              round(file.info(best)$size/1024^2, 1), " MB).")
      auto_paths[[n]] <- best
      next
    }
    message("  [FETCH] Downloading ", n, " (.zip)...")
    dl_result <- tryCatch(download.file(url, zip_dest, method = "libcurl", quiet = TRUE),
                          error = function(e) { message("  [ERROR] ", e$message); -1L })
    if (!file.exists(zip_dest) || file.info(zip_dest)$size < 1024^2) {
      message("  [ERROR] ZIP download failed for ", n); next
    }
    if (dir.exists(extract_dir)) unlink(extract_dir, recursive = TRUE)
    dir.create(extract_dir)
    unzip_result <- tryCatch(unzip(zip_dest, exdir = extract_dir),
                             error = function(e) { message("  [ERROR] ", e$message); NULL })
    if (is.null(unzip_result)) { unlink(zip_dest); next }
    all_extracted <- list.files(extract_dir, full.names = TRUE, recursive = TRUE)
    data_files    <- all_extracted[
      !grepl("readme", basename(all_extracted), ignore.case = TRUE)]
    if (length(data_files) == 0) {
      message("  [ERROR] No data files after extraction for ", n)
      unlink(zip_dest); unlink(extract_dir, recursive = TRUE); next
    }
    message("  [ZIP] Candidates:")
    for (f in data_files)
      message("        ", basename(f), " (", round(file.info(f)$size/1024^2,1), " MB)")
    target_file <- data_files[which.max(file.info(data_files)$size)]
    message("  [ZIP] Selected: ", basename(target_file))
    # Build dest: TraitAcronym_OriginalFilename (add .gz suffix if not compressed)
    out_base <- basename(target_file)
    if (!grepl("\\.(gz|bgz)$", out_base, ignore.case = TRUE))
      out_base <- paste0(out_base, ".gz")
    dest <- file.path(gwas_dir, paste0(n, "_", out_base))
    ok <- recompress_to_gz(target_file, dest, n)
    unlink(zip_dest); unlink(extract_dir, recursive = TRUE)
    if (!ok) next
    zip_sourced <- c(zip_sourced, n)
    check_size(dest, n)
    auto_paths[[n]] <- dest
    next
  }

  dest <- file.path(gwas_dir, paste0(n, "_", basename(url)))

  # Integrity check
  if (file.exists(dest)) {
    local_size  <- file.info(dest)$size
    response    <- tryCatch(httr::HEAD(url), error = function(e) NULL)
    server_size <- safe_content_length(response)
    size_ok     <- local_size >= 1024^2 &&
                   (is.na(server_size) ||
                    abs(local_size - server_size) <= server_size * 0.01)
    if (size_ok) {
      message("  [SKIP] ", n, " (", basename(dest), ", ",
              round(local_size/1024^2, 1), " MB).")
      auto_paths[[n]] <- dest
      next
    }
    message("  [!] ", n, " incomplete. Re-downloading..."); unlink(dest)
  }

  message("  [FETCH] Downloading ", n, " -> ", basename(dest), "...")
  dl_result <- tryCatch(download.file(url, dest, method = "libcurl", quiet = TRUE),
                        error = function(e) { message("  [ERROR] ", e$message); -1L })
  check_size(dest, n)
  auto_paths[[n]] <- dest
}

elapsed <- proc.time() - step_time
message(sprintf("  [DONE] Step 2 completed in %.1f seconds.", elapsed["elapsed"]))

# =============================================================================
# § 3 — CENSUS CHECK
# =============================================================================

step_time <- proc.time()
message("\n", strrep("=", 70))
message("[Step 3] Census check...")
message(strrep("=", 70))

missing_auto <- missing_manual <- character(0)

# Only check traits that are declared as Exposures or Outcomes
for (n in c(Exposures, Outcomes)) {
  is_manual <- nchar(trimws(all_traits[[n]])) == 0
  if (is_manual) {
    if (!is.null(manual_paths[[n]])) {
      src_size <- file.info(manual_paths[[n]])$size
      message("  [MANUAL] ", n, " — ", basename(manual_paths[[n]]),
              " (", round(src_size/1024^2, 1), " MB)")
    } else {
      missing_manual <- c(missing_manual, n)
      message("  [MISSING-MANUAL] ", n)
    }
  } else {
    ap <- auto_paths[[n]]
    if (!is.null(ap) && file.exists(ap) && file.info(ap)$size >= 1000) {
      message("  [PRESENT] ", n, " — ", basename(ap),
              " (", round(file.info(ap)$size/1024^2, 1), " MB)")
    } else {
      missing_auto <- c(missing_auto, n)
      message("  [MISSING-AUTO]   ", n)
    }
  }
}

if (!file.exists(annot_file) || file.info(annot_file)$size < 1024^2) {
  missing_auto <- c(missing_auto, "_ANNOT")
  message("  [MISSING-AUTO] Annotation file: ", annot_file)
} else {
  message("  [PRESENT] Annotation file (",
          round(file.info(annot_file)$size/1024^2, 1), " MB)")
}

if (length(missing_auto) > 0) {
  cat("\n", strrep("*", 60), "\n")
  cat("STOP: AUTOMATED DOWNLOAD FAILURES\n\n")
  for (m in missing_auto) cat("  - ", m, "\n")
  cat(strrep("*", 60), "\n\n")
  stop("Census Check Failed: ", length(missing_auto), " file(s) missing.")
}
if (length(missing_manual) > 0) {
  cat("\n", strrep("-", 60), "\n")
  cat("ACTION REQUIRED — manual files not yet placed:\n\n")
  for (m in missing_manual)
    cat("  - Place file starting with '", m, "' in '", gwas_dir, "/'\n", sep="")
  cat("Accepted extensions: .txt .tsv .csv .gz .bgz .tsv.gz\n")
  cat("Processing will continue for available traits.\n")
  cat(strrep("-", 60), "\n\n")
}

message("  [OK] Census complete.")
elapsed <- proc.time() - step_time
message(sprintf("  [DONE] Step 3 completed in %.1f seconds.", elapsed["elapsed"]))

# =============================================================================
# § 5 — LOAD AND MERGE ALL EXPOSURES (wide format)
# =============================================================================

step_time <- proc.time()
message("\n", strrep("=", 70))
message("[Step 4] Loading and merging exposures...")
message(strrep("=", 70))

# Load UKB annotation file
message("  [READ] UKB annotation file...")
annot_dt <- tryCatch(
  as.data.table(suppressMessages(suppressWarnings(
    vroom::vroom(annot_file,
                 col_select    = c("variant","chr","pos","rsid","ref","alt"),
                 altrep = FALSE,
                 show_col_types = FALSE, progress = TRUE)
  ))),
  error = function(e) stop("FATAL: Cannot read annotation file: ", e$message)
)
message("  [OK] Annotation: ", format(nrow(annot_dt), big.mark=","), " variants.")

data_merge <- as.data.frame(annot_dt) %>%
  tidyr::unite("MarkerName", chr:pos, sep=":")
message("  [OK] Backbone: ", format(nrow(data_merge), big.mark=","), " variants.")

ordered_exposures <- if (!is.null(exposure_order)) {
  stopifnot("exposure_order must be a permutation of Exposures" =
              setequal(exposure_order, Exposures))
  exposure_order
} else {
  sort_traits(Exposures)
}
message("  [ORDER] Exposure processing order: ",
        paste(ordered_exposures, collapse=" → "))

withCallingHandlers({
for (trait in ordered_exposures) {
  message("\n  --- Exposure: ", trait, " ---")
  trait_path <- get_trait_path(trait)
  if (!file.exists(trait_path)) { message("  [SKIP] File not found: ", trait_path); next }

  dt <- read_trait_lazy(trait_path, trait)
  if (is.null(dt) || nrow(dt) == 0) { message("  [SKIP] Empty file."); next }
  message("  [OK] ", format(nrow(dt), big.mark=","), " rows loaded.")

  cols <- detect_columns(names(dt), trait)
  if (is.null(cols)) next

  tl <- tolower(trait)

  if (cols$chrpos_method == "variant_field") {
    # ---- UKB / Neale Lab: merge via annotation file by variant key ----------
    # Beta is NOT re-oriented here, deliberately (A6FIX conclusion, Article_6, 2026-08-29):
    # Neale Lab's beta is unconditionally the effect of the ALT allele (4th field of
    # `variant` = chr:pos:ref:alt), whichever allele is the minor one. `minor_allele` is
    # not this file's effect allele, so flipping beta by it (the 2026-08-26 "fix") corrupted
    # ~27% of SNPs; it was reverted. The A6FIX re-run reproduced this pipeline's MrDAG_data exactly.
    # Lazy tibble subset (preserves ALTREP) — only the 4 needed columns
    dt_sel <- as.data.frame(dt[, c("variant", cols$beta_col,
                                    cols$se_col, cols$pval_col),
                                drop = FALSE])
    names(dt_sel)[names(dt_sel) == cols$beta_col] <- paste0(tl, "_beta")
    names(dt_sel)[names(dt_sel) == cols$se_col]   <- paste0(tl, "_se")
    names(dt_sel)[names(dt_sel) == cols$pval_col] <- paste0(tl, "_pval")

    message("  [MERGE] ", trait, " → backbone via variant key (chr:pos:ref:alt)")
    data_merge <- merge(data_merge, dt_sel, by = "variant", sort = FALSE)
    message("  [MERGE] After ", trait, ": ",
            format(nrow(data_merge), big.mark = ","), " rows.")

  } else {

    # Effect_A2 correction — before renaming
    if (!is.na(cols$beta_col) && cols$beta_col == "Effect_A2" &&
        !"Effect_A1" %in% names(dt)) {
      tmp         <- cols$ea_col
      cols$ea_col <- cols$oa_col
      cols$oa_col <- tmp
      message("  [CORRECT] Effect_A2 without Effect_A1: treating A2 as effect allele.")
    }

    if (!is.na(cols$beta_col) && cols$beta_col == "OR") {
      or_raw <- suppressWarnings(as.numeric(dt[[cols$beta_col]]))
      if (any(or_raw < 0, na.rm = TRUE)) {
        message("  [CORRECT] Column 'OR' has negative values — header error: ",
                "already log(OR). Renaming to 'logOR', no transformation applied.")
      } else {
        message("  [CORRECT] Column 'OR' is all positive — genuine OR: ",
                "log-transforming and renaming to 'logOR'.")
        dt[[cols$beta_col]] <- log(or_raw)
      }
      names(dt)[names(dt) == "OR"] <- "logOR"
      cols$beta_col <- "logOR"
    }

    has_rsid <- !is.na(cols$rsid_col) &&
                cols$rsid_col %in% names(dt) &&
                sum(grepl("^rs",
                          iconv(as.character(dt[[cols$rsid_col]]),
                                to = "UTF-8", sub = "")),
                    na.rm = TRUE) > 0

    if (has_rsid) {
      rsid_vec  <- extract_rsid(as.data.table(dt[, cols$rsid_col, drop = FALSE]),
                                cols$rsid_col)
      keep_rows <- !is.na(rsid_vec) & grepl("^rs", rsid_vec) &
                   rsid_vec %in% data_merge$rsid
      dt        <- dt[keep_rows, , drop = FALSE]
      rsid_vec  <- rsid_vec[keep_rows]
      dt        <- dt[!duplicated(rsid_vec), , drop = FALSE]
      rsid_vec  <- rsid_vec[!duplicated(rsid_vec)]

      dt_sel <- data.frame(
        rsid = rsid_vec,
        EA   = as.character(dt[[cols$ea_col]]),
        NEA  = if (!is.na(cols$oa_col) && cols$oa_col %in% names(dt))
                 as.character(dt[[cols$oa_col]]) else NA_character_,
        stringsAsFactors = FALSE
      )
      dt_sel[[paste0(tl, "_beta_not_aligned")]] <- suppressWarnings(as.numeric(dt[[cols$beta_col]]))
      dt_sel[[paste0(tl, "_se")]]               <- suppressWarnings(as.numeric(dt[[cols$se_col]]))
      dt_sel[[paste0(tl, "_pval")]]             <- suppressWarnings(as.numeric(dt[[cols$pval_col]]))
      if (!is.na(cols$eaf_col) && cols$eaf_col %in% names(dt))
        dt_sel[[paste0(tl, "_eaf")]]            <- suppressWarnings(as.numeric(dt[[cols$eaf_col]]))

      message("  [MERGE] ", trait, " → inner join rsid (via ", cols$rsid_col, ")")
      data_merge <- merge(data_merge, dt_sel, by = "rsid", sort = FALSE)

    } else {
      chrpos_vec <- build_chrpos(as.data.table(dt), cols$chrpos_method)
      keep_rows  <- chrpos_vec %in% data_merge$MarkerName
      dt         <- dt[keep_rows, , drop = FALSE]
      chrpos_vec <- chrpos_vec[keep_rows]
      dt         <- dt[!duplicated(chrpos_vec), , drop = FALSE]
      chrpos_vec <- chrpos_vec[!duplicated(chrpos_vec)]

      dt_sel <- data.frame(
        MarkerName = chrpos_vec,
        EA         = as.character(dt[[cols$ea_col]]),
        NEA        = if (!is.na(cols$oa_col) && cols$oa_col %in% names(dt))
                       as.character(dt[[cols$oa_col]]) else NA_character_,
        stringsAsFactors = FALSE
      )
      dt_sel[[paste0(tl, "_beta_not_aligned")]] <- suppressWarnings(as.numeric(dt[[cols$beta_col]]))
      dt_sel[[paste0(tl, "_se")]]               <- suppressWarnings(as.numeric(dt[[cols$se_col]]))
      dt_sel[[paste0(tl, "_pval")]]             <- suppressWarnings(as.numeric(dt[[cols$pval_col]]))
      if (!is.na(cols$eaf_col) && cols$eaf_col %in% names(dt))
        dt_sel[[paste0(tl, "_eaf")]]            <- suppressWarnings(as.numeric(dt[[cols$eaf_col]]))

      message("  [MERGE] ", trait, " → inner join MarkerName (via ", cols$chrpos_method, ")")
      data_merge <- merge(data_merge, dt_sel, by = "MarkerName", sort = FALSE)
    }
    rm(dt); gc()

    message("  [MERGE] After ", trait, ": ",
            format(nrow(data_merge), big.mark = ","), " rows.")

    need_toupper <- cols$chrpos_method %in% c("MarkerName_col", "markername_col",
                                               "chr_bp_hg19", "base_cols", "CHR_BP")
    data_merge <- harmonise_allele(data_merge, trait, use_toupper = need_toupper)
  }
}
}, warning = function(w) {
  if (grepl("parsing issues", conditionMessage(w), ignore.case = TRUE)) {
    message("  [PARSE-WARN] vroom: ", conditionMessage(w),
            " — affected cells set to NA, data otherwise unchanged.")
    invokeRestart("muffleWarning")
  }
})

rm(annot_dt); gc()

if (is.null(data_merge) || nrow(data_merge) == 0)
  stop("FATAL: No exposure data could be merged.")

id_cols  <- intersect(c("rsid","MarkerName","variant","ref","alt"),
                      names(data_merge))
exp_cols <- unlist(lapply(tolower(Exposures), function(t)
              c(paste0(t,"_beta"), paste0(t,"_se"), paste0(t,"_pval"))))
exp_cols <- exp_cols[exp_cols %in% names(data_merge)]
data_merge <- data_merge[, c(id_cols, exp_cols), drop=FALSE]

message("\n  [OK] Merged exposure dataset: ",
        format(nrow(data_merge), big.mark=","), " rows, ",
        ncol(data_merge), " columns.")

elapsed <- proc.time() - step_time
message(sprintf("  [DONE] Step 4 completed in %.1f seconds.", elapsed["elapsed"]))

# =============================================================================
# § 6 — INSTRUMENT SELECTION
# =============================================================================

step_time <- proc.time()
message("\n", strrep("=", 70))
message("[Step 5] Instrument selection (min p-value across all exposures)...")
message(strrep("=", 70))

# Detect pval columns that are actually present (some traits may have been skipped)
pval_cols_present <- paste0(tolower(Exposures), "_pval")
pval_cols_present <- pval_cols_present[pval_cols_present %in% names(data_merge)]

missing_pval <- setdiff(paste0(tolower(Exposures), "_pval"), names(data_merge))
if (length(missing_pval) > 0)
  message("  [WARNING] Missing pval columns (trait skipped?): ",
          paste(missing_pval, collapse=", "))

pmat    <- data_merge %>% dplyr::select(all_of(pval_cols_present))
min_pval <- apply(pmat, MARGIN=1, FUN=min, na.rm=TRUE)
rm(pmat)

message("  [INFO] Variants with min p < ", p_thresh, " : ",
        format(sum(min_pval < p_thresh), big.mark=","))

data_merge$min.pval <- min_pval
IV <- data_merge[min_pval < p_thresh, ]
rm(data_merge, min_pval); gc()

if (nrow(IV) == 0) stop("FATAL: No variants pass the significance threshold.")
message("  [OK] Instruments selected: ", format(nrow(IV), big.mark=","))

# Set row names to rsid only if all values are unique valid RSIDs
if (!anyDuplicated(IV$rsid) && !any(is.na(IV$rsid)) &&
    all(grepl("^rs", IV$rsid))) {
  rownames(IV) <- IV$rsid
  message("  [OK] Row names set to rsid.")
} else {
  n_dup <- sum(duplicated(IV$rsid))
  n_na  <- sum(is.na(IV$rsid))
  n_bad <- sum(!grepl("^rs", IV$rsid[!is.na(IV$rsid)]))
  message("  [WARNING] Row names not set: ",
          n_dup, " duplicates, ", n_na, " NAs, ",
          n_bad, " non-rs identifiers.")
}
local({
  IV <- IV
  names(IV) <- standardise_trait_names(names(IV))
  save(IV, file = output_iv_rdata)
})
message("  [OK] IV saved: ", output_iv_rdata)

elapsed <- proc.time() - step_time
message(sprintf("  [DONE] Step 5 completed in %.1f seconds.", elapsed["elapsed"]))

# =============================================================================
# § 7 — LOAD AND MERGE ALL OUTCOMES (wide format)
# =============================================================================

step_time <- proc.time()
message("\n", strrep("=", 70))
message("[Step 6] Loading and merging outcomes (left joins)...")
message(strrep("=", 70))

IV_full        <- IV    # data.frame — preserved for independent left joins
outcome_frames <- list()

ordered_outcomes <- if (!is.null(outcome_order)) {
  stopifnot("outcome_order must be a permutation of Outcomes" =
              setequal(outcome_order, Outcomes))
  outcome_order
} else {
  sort_traits(Outcomes)
}
message("  [ORDER] Outcome processing order: ",
        paste(ordered_outcomes, collapse=" → "))

withCallingHandlers({
for (trait in ordered_outcomes) {
  message("\n  --- Outcome: ", trait, " ---")
  trait_path <- get_trait_path(trait)
  if (!file.exists(trait_path)) { message("  [SKIP] File not found: ", trait_path); next }

  dt <- read_trait_lazy(trait_path, trait)
  if (is.null(dt) || nrow(dt) == 0) { message("  [SKIP] Empty file."); next }
  message("  [OK] ", format(nrow(dt), big.mark=","), " rows loaded.")

  cols <- detect_columns(names(dt), trait)
  if (is.null(cols)) next

  tl <- tolower(trait)

  if (cols$chrpos_method == "variant_field") {
    if (!"variant" %in% names(IV_full)) {
      message("  [WARNING] ", trait, " is UKB-format (variant key) but IV_full ",
              "has no `variant` column to merge on. Excluding.")
      next
    }

    df_sel <- as.data.frame(dt[, c("variant", cols$beta_col, cols$se_col),
                                drop = FALSE])
    names(df_sel)[names(df_sel) == cols$beta_col] <- paste0(tl, "_beta")
    names(df_sel)[names(df_sel) == cols$se_col]   <- paste0(tl, "_se")

    # Beta is already ALT-oriented (Neale Lab) -- see the exposures-loop note above.
    message("  [MERGE] ", trait, " → left join variant (direct UKB match, ",
            "no harmonisation needed)")
    merged <- merge(IV_full[, "variant", drop = FALSE], df_sel,
                    by = "variant", all.x = TRUE, sort = FALSE)

    n_matched <- sum(!is.na(merged[[paste0(tl, "_beta")]]))
    n_missing <- nrow(merged) - n_matched
    message("  [MERGE] ", trait, ": ", format(n_matched, big.mark = ","),
            " matched | ", format(n_missing, big.mark = ","), " missing.")

    if (n_matched == 0) {
      message("  [WARNING] No variants matched for ", trait, ". Excluding.")
      next
    }

    keep_out <- c("variant",
                  intersect(c(paste0(tl,"_beta"), paste0(tl,"_se")), names(merged)))
    outcome_frames[[trait]] <- merged[, keep_out, drop = FALSE]
    next
  }

  # Effect_A2 correction — before renaming (same as the exposures loop above;
  # was previously missing here, which is what silently negated AF: its GWAS
  # file reports Effect_A2, so without this swap the effect allele was taken
  # to be A1 instead of A2, inverting the sign of every harmonised beta)
  if (!is.na(cols$beta_col) && cols$beta_col == "Effect_A2" &&
      !"Effect_A1" %in% names(dt)) {
    tmp         <- cols$ea_col
    cols$ea_col <- cols$oa_col
    cols$oa_col <- tmp
    message("  [CORRECT] Effect_A2 without Effect_A1: treating A2 as effect allele.")
  }

  if (!is.na(cols$beta_col) && cols$beta_col == "OR") {
    or_raw <- suppressWarnings(as.numeric(dt[[cols$beta_col]]))
    if (any(or_raw < 0, na.rm = TRUE)) {
      message("  [CORRECT] Column 'OR' has negative values — header error: ",
              "already log(OR). Renaming to 'logOR', no transformation applied.")
    } else {
      message("  [CORRECT] Column 'OR' is all positive — genuine OR: ",
              "log-transforming and renaming to 'logOR'.")
      dt[[cols$beta_col]] <- log(or_raw)
    }
    names(dt)[names(dt) == "OR"] <- "logOR"
    cols$beta_col <- "logOR"
  }

  # ── Determine merge key BEFORE the join
  # has_rsid checks the raw rsid column for rs-prefixed values
  has_rsid <- !is.na(cols$rsid_col) &&
              cols$rsid_col %in% names(dt) &&
              sum(grepl("^rs",
                        iconv(as.character(dt[[cols$rsid_col]]),
                              to = "UTF-8", sub = "")),
                  na.rm = TRUE) > 0

  if (has_rsid) {
    # ── rsid path: extract key, pre-filter to IV rsids, join
    rsid_vec <- extract_rsid(as.data.table(dt[, cols$rsid_col, drop = FALSE]),
                             cols$rsid_col)
    keep_rows <- !is.na(rsid_vec) & grepl("^rs", rsid_vec) &
                 rsid_vec %in% IV_full$rsid
    dt        <- dt[keep_rows, , drop = FALSE]
    rsid_vec  <- rsid_vec[keep_rows]
    dt        <- dt[!duplicated(rsid_vec), , drop = FALSE]
    rsid_vec  <- rsid_vec[!duplicated(rsid_vec)]

    # Build narrow df_sel — rename before merge to prevent column collisions
    df_sel <- data.frame(
      rsid = rsid_vec,
      EA   = as.character(dt[[cols$ea_col]]),
      NEA  = if (!is.na(cols$oa_col) && cols$oa_col %in% names(dt))
               as.character(dt[[cols$oa_col]]) else NA_character_,
      stringsAsFactors = FALSE
    )
    df_sel[[paste0(tl, "_beta_not_aligned")]] <- suppressWarnings(as.numeric(dt[[cols$beta_col]]))
    df_sel[[paste0(tl, "_se")]]               <- suppressWarnings(as.numeric(dt[[cols$se_col]]))
    if (!is.na(cols$pval_col) && cols$pval_col %in% names(dt))
      df_sel[[paste0(tl, "_pval")]]           <- suppressWarnings(as.numeric(dt[[cols$pval_col]]))

    message("  [MERGE] ", trait, " → left join rsid (via ", cols$rsid_col, ")")
    merged <- merge(IV_full[, "rsid", drop = FALSE], df_sel,
                    by = "rsid", all.x = TRUE, sort = FALSE)

  } else {
    chrpos_vec <- build_chrpos(as.data.table(dt), cols$chrpos_method)
    keep_rows  <- chrpos_vec %in% IV_full$MarkerName
    dt         <- dt[keep_rows, , drop = FALSE]
    chrpos_vec <- chrpos_vec[keep_rows]
    dt         <- dt[!duplicated(chrpos_vec), , drop = FALSE]
    chrpos_vec <- chrpos_vec[!duplicated(chrpos_vec)]

    df_sel <- data.frame(
      MarkerName = chrpos_vec,
      EA         = as.character(dt[[cols$ea_col]]),
      NEA        = if (!is.na(cols$oa_col) && cols$oa_col %in% names(dt))
                     as.character(dt[[cols$oa_col]]) else NA_character_,
      stringsAsFactors = FALSE
    )
    df_sel[[paste0(tl, "_beta_not_aligned")]] <- suppressWarnings(as.numeric(dt[[cols$beta_col]]))
    df_sel[[paste0(tl, "_se")]]               <- suppressWarnings(as.numeric(dt[[cols$se_col]]))
    if (!is.na(cols$pval_col) && cols$pval_col %in% names(dt))
      df_sel[[paste0(tl, "_pval")]]           <- suppressWarnings(as.numeric(dt[[cols$pval_col]]))

    message("  [MERGE] ", trait, " → left join MarkerName (via ", cols$chrpos_method, ")")
    merged <- merge(IV_full[, "MarkerName", drop = FALSE], df_sel,
                    by = "MarkerName", all.x = TRUE, sort = FALSE)
  }
  rm(dt); gc()

  # ── All heavy ops now run on ~150K rows -----------------------------------
  # Columns are already renamed in df_sel above — no post-merge renaming needed.

  n_matched <- sum(!is.na(merged[[paste0(tl, "_beta_not_aligned")]]))
  n_missing <- nrow(merged) - n_matched
  message("  [MERGE] ", trait, ": ", format(n_matched, big.mark = ","),
          " matched | ", format(n_missing, big.mark = ","), " missing.")

  if (n_matched == 0) {
    message("  [WARNING] No variants matched for ", trait, ". Excluding.")
    next
  }

  # Allele harmonisation — attach alt from IV_full, harmonise, remove alt
  if (has_rsid) {
    merged$alt <- IV_full$alt[match(merged$rsid,       IV_full$rsid)]
  } else {
    merged$alt <- IV_full$alt[match(merged$MarkerName, IV_full$MarkerName)]
  }
  need_toupper <- cols$chrpos_method %in% c("MarkerName_col","chr_bp_hg19","base_cols")
  merged <- harmonise_allele(merged, trait, use_toupper=need_toupper)
  merged$alt <- NULL

  # Store only the merge key + beta/se columns — discard everything else
  tl_beta <- paste0(tl, "_beta"); tl_se <- paste0(tl, "_se")
  key_col <- if (has_rsid) "rsid" else "MarkerName"
  keep_out <- c(key_col,
                intersect(c(tl_beta, tl_se), names(merged)))
  outcome_frames[[trait]] <- merged[, keep_out, drop=FALSE]
}
}, warning = function(w) {
  if (grepl("parsing issues", conditionMessage(w), ignore.case = TRUE)) {
    message("  [PARSE-WARN] vroom: ", conditionMessage(w),
            " — affected cells set to NA, data otherwise unchanged.")
    invokeRestart("muffleWarning")
  }
})

if (length(outcome_frames) == 0)
  stop("FATAL: No outcomes could be merged.")

# Bind all per-outcome beta/se columns onto IV_full — pure data.frame operations
IV <- IV_full

for (trait in names(outcome_frames)) {
  tl       <- tolower(trait)
  fr       <- outcome_frames[[trait]]           # data.frame
  new_cols <- c(paste0(tl, "_beta"), paste0(tl, "_se"))
  new_cols <- new_cols[new_cols %in% names(fr)]

  if ("rsid" %in% names(fr)) {
    idx <- match(IV$rsid, fr$rsid)
  } else if ("variant" %in% names(fr)) {
    idx <- match(IV$variant, fr$variant)
  } else {
    idx <- match(IV$MarkerName, fr$MarkerName)
  }
  IV <- cbind(IV, fr[idx, new_cols, drop=FALSE])
}

# Drop variants missing from any outcome — pure data.frame
outcome_beta_cols <- paste0(tolower(names(outcome_frames)), "_beta")
outcome_beta_cols <- outcome_beta_cols[outcome_beta_cols %in% names(IV)]
complete_rows     <- complete.cases(IV[, outcome_beta_cols, drop=FALSE])
n_dropped         <- sum(!complete_rows)
IV                <- IV[complete_rows, , drop=FALSE]

message("  [OK] Dropped (missing >=1 outcome): ",
        format(n_dropped, big.mark=","))
message("  [OK] Final IV: ",
        format(nrow(IV), big.mark=","), " rows, ", ncol(IV), " columns.")

elapsed <- proc.time() - step_time
message(sprintf("  [DONE] Step 6 completed in %.1f seconds.", elapsed["elapsed"]))

# =============================================================================
# § 8 — LD CLUMPING
# =============================================================================

step_time <- proc.time()
message("\n", strrep("=", 70))
message("[Step 7] LD clumping (r2 = ", clump_r2, ", kb = ", clump_kb, ")...")
message(strrep("=", 70))

IV <- IV %>% dplyr::rename(pval = min.pval)

clump_input <- IV %>%
  dplyr::select(rsid, pval) %>%
  dplyr::filter(!is.na(rsid) & grepl("^rs", rsid)) %>%
  dplyr::distinct(rsid, .keep_all=TRUE)

n_no_rsid <- nrow(IV) - nrow(clump_input)
if (n_no_rsid > 0)
  message("  [INFO] Excluded (no valid RSID)   : ", format(n_no_rsid, big.mark=","))

if (nrow(clump_input) == 0)
  stop("FATAL: No RSID-tagged instruments for clumping. ",
       "Check that outcome MarkerName formats match IV MarkerName.")

message("\n  [CLUMP] Single pass on ",
        format(nrow(clump_input), big.mark=","), " RSIDs...")
clumped_final <- tryCatch(
  ieugwasr::ld_clump(clump_input,
                     clump_kb  = clump_kb,
                     clump_r2  = clump_r2,
                     clump_p   = clump_p,
                     plink_bin = plink_path,
                     bfile     = ref_path_for_clumping),
  error = function(e) stop("FATAL: ld_clump() failed: ", e$message)
)
gc()

message("  [OK] Independent instruments: ",
        format(nrow(clumped_final), big.mark=","))
if (nrow(clumped_final) == 0) stop("FATAL: No instruments after clumping.")

IV_clump <- IV %>% dplyr::filter(rsid %in% clumped_final$rsid)
rownames(IV_clump) <- IV_clump$rsid
IV_clump$min_pval  <- IV$pval[match(IV_clump$rsid, IV$rsid)]
rm(clumped_final, clump_input); gc()

message("  [OK] IV_clump: ",
        format(nrow(IV_clump), big.mark=","), " rows.")
local({
  IV_clump <- IV_clump
  names(IV_clump) <- standardise_trait_names(names(IV_clump))
  save(IV_clump, file = output_clumped_rdata)
})
message("  [OK] IV_clump saved: ", output_clumped_rdata)

} # end of `if (!checkpoint_valid)` -- everything from §2 (download) through

f_stat_traits <- Exposures[beta_X_cols %in% names(IV_clump) &
                           se_X_cols   %in% names(IV_clump)]

if (length(f_stat_traits) < 2) {
  message("  [SKIP] Conditional F-statistics require >= 2 exposures with beta and se columns present.")
  cond_F_named <- NULL
} else {
  f_stat_beta_cols <- beta_X_cols[Exposures %in% f_stat_traits]
  f_stat_se_cols   <- se_X_cols[Exposures %in% f_stat_traits]
  betaX_f <- to_num_matrix(IV_clump, f_stat_beta_cols)
  seX_f   <- to_num_matrix(IV_clump, f_stat_se_cols)
  colnames(betaX_f) <- f_stat_traits
  colnames(seX_f)   <- f_stat_traits

  # format_mvmr() needs a placeholder outcome — any valid vector of the
  # correct length works, since strength_mvmr() never reads it.
  placeholder_outcome <- tolower(Outcomes)[paste0(tolower(Outcomes), "_beta") %in% names(IV_clump)][1]
  placeholder_beta <- suppressWarnings(as.numeric(IV_clump[[paste0(placeholder_outcome, "_beta")]]))
  placeholder_se   <- suppressWarnings(as.numeric(IV_clump[[paste0(placeholder_outcome, "_se")]]))

  gencov_is_ldsc      <- FALSE
  gencov_ldsc_partial <- FALSE
  all_ldsc_traits <- c(f_stat_traits, Outcomes)
  ldsc_traits_ok  <- all_ldsc_traits %in% names(trait_sample_info)
  if (!all(ldsc_traits_ok)) {
    message("  [WARNING] Missing trait_sample_info for: ",
            paste(all_ldsc_traits[!ldsc_traits_ok], collapse = ", "),
            " — falling back to gencov = 0 for all exposures.")
    gencov_arg <- 0
  } else {
    ldsc_result <- tryCatch(
      estimate_ldsc_gencov(
        traits        = all_ldsc_traits,
        trait_paths   = setNames(sapply(all_ldsc_traits, get_trait_path), all_ldsc_traits),
        sample_info = trait_sample_info,
        munge_dir   = ldsc_munged_dir,
        ld_dir      = ld_dir,
        hm3_snplist = hm3_snplist,
        annot_file  = annot_file
      ),
      error = function(e) {
        message("  [WARNING] LDSC sample-overlap estimation failed: ", e$message,
                " — falling back to gencov = 0.")
        NULL
      }
    )
    if (is.null(ldsc_result)) {
      gencov_arg <- 0
    } else {
      ldsc_intercept <- ldsc_result$I           # Z1*Z2 scale, raw
      ldsc_gencov    <- ldsc_result$Ns_gencov   # Ns*gencov scale
      ldsc_gencov_S  <- ldsc_result$S
      ldsc_rg        <- ldsc_result$rg
      ldsc_rg_se     <- ldsc_result$rg_se
      ldsc_rg_p      <- ldsc_result$rg_p
      ldsc_trait_N   <- ldsc_result$N
      ldsc_traits    <- ldsc_result$trait_names

      present_exposures <- intersect(f_stat_traits, ldsc_traits)
      missing_exposures <- setdiff(f_stat_traits, ldsc_traits)
      if (length(missing_exposures) > 0) {
        message("  [WARNING] LDSC could not process exposure(s): ",
                paste(missing_exposures, collapse = ", "),
                " — gencov = 0 for any pair involving ",
                if (length(missing_exposures) == 1) "this exposure" else "these exposures",
                "; LDSC-based sample-overlap is still used for the remaining ",
                length(present_exposures), " exposure(s) (all estimates ",
                "are saved in LDSC_output.RData above regardless).")
      }

      # Build the full f_stat_traits x f_stat_traits correlation matrix:
      # LDSC-based for pairs where both exposures were successfully
      # estimated, 0 (i.e. assumed independent) for any pair involving an
      # exposure LDSC could not process. Diagonal is always 1.
      exposure_corr <- matrix(0, length(f_stat_traits), length(f_stat_traits),
                              dimnames = list(f_stat_traits, f_stat_traits))
      diag(exposure_corr) <- 1
      if (length(present_exposures) >= 2) {
        exposure_corr[present_exposures, present_exposures] <-
          ldsc_result$full_matrix[present_exposures, present_exposures]
      }

      eig_min <- min(eigen(exposure_corr, symmetric = TRUE, only.values = TRUE)$values)
      if (eig_min < -1e-8) {
        message("  [WARNING] LDSC-estimated exposure correlation matrix is not ",
                "positive-semi-definite (smallest eigenvalue = ", .display_round(eig_min, 3),
                "). Projecting onto the nearest valid correlation matrix ",
                "before use.")
        exposure_corr <- as.matrix(Matrix::nearPD(exposure_corr, corr = TRUE,
                                                   keepDiag = FALSE)$mat)
      }

      gencov_arg          <- MVMR::phenocov_mvmr(exposure_corr, seX_f)
      gencov_is_ldsc      <- length(present_exposures) >= 2
      gencov_ldsc_partial <- gencov_is_ldsc && length(missing_exposures) > 0
      message("  [OK] LDSC-based sample-overlap covariance estimated for ",
              length(present_exposures), " of ", length(f_stat_traits),
              " exposures", if (length(missing_exposures) > 0) " (see above)." else ".")

      ldsc_gencov_overlap <- gencov_arg
      ldsc_intercept_se   <- ldsc_result$I_se                    # block-jackknife SE of each intercept
      ldsc_intercept_p    <- 2 * stats::pnorm(-abs(ldsc_intercept / ldsc_intercept_se))   # two-sided P
      diag(ldsc_intercept_p) <- NA
      save(ldsc_intercept, ldsc_intercept_se, ldsc_intercept_p, ldsc_gencov, ldsc_gencov_S, ldsc_rg,
           ldsc_rg_se, ldsc_rg_p, ldsc_trait_N, ldsc_traits, ldsc_gencov_overlap,
           file = output_ldsc_gencov_rdata)
      message("  [OK] LDSC intercept + Ns*gencov + genetic covariance/",
              "correlation (+ SE/p) + ldsc_gencov_overlap matrices saved: ",
              output_ldsc_gencov_rdata)
    }
  }

  cond_F_named <- tryCatch({
    mvmr_input <- MVMR::format_mvmr(
      BXGs   = as.data.frame(betaX_f),
      BYG    = placeholder_beta,
      seBXGs = as.data.frame(seX_f),
      seBYG  = placeholder_se,
      RSID   = IV_clump$rsid
    )
    invisible(capture.output(
      cond_F_raw <- suppressWarnings(
        MVMR::strength_mvmr(r_input = mvmr_input, gencov = gencov_arg)
      )
    ))
    f_vals <- as.numeric(unlist(cond_F_raw[1, ]))
    setNames(f_vals, f_stat_traits)   # columns of cond_F_raw follow BXGs order
  }, error = function(e) {
    message("  [WARNING] Conditional F-statistic computation failed: ", e$message)
    NULL
  })
}

if (!is.null(cond_F_named)) {
  weak_traits <- names(cond_F_named)[cond_F_named < 10]
  if (length(weak_traits) > 0)
    message("  [WARNING] Conditional F-statistic < 10 (weak instrument): ",
            paste(weak_traits, collapse = ", "))
  else
    message("  [OK] All conditional F-statistics >= 10.")
}

exp_pval_cols <- setNames(pval_X_cols, Exposures)
exp_pval_cols <- exp_pval_cols[exp_pval_cols %in% names(IV)]

if (length(exp_pval_cols) == 0)
  stop("FATAL: none of the current Exposures (", paste(Exposures, collapse = ", "),
       ") have a matching _pval column in IV. This usually means IV was ",
       "generated for a different Exposures configuration than the one ",
       "currently declared in §0 -- check that Exposures matches what ",
       "this IV/IV_clump was actually built from.")

total_IVs <- nrow(IV_clump)

## Own_Selected_IVs: how many of the panel's IVs (AFTER clumping, from
## IV_clump -- not IV, which is the pre-clumping GW_Sig_Variants count) are
## this exposure's OWN genome-wide-significant hits, vs. having entered the
## pooled panel only because a sibling exposure was significant there.
iv_summary <- do.call(rbind, lapply(names(exp_pval_cols), function(exp) {
  pcol  <- exp_pval_cols[[exp]]
  pvals <- suppressWarnings(as.numeric(IV[[pcol]]))
  pvals_clumped <- suppressWarnings(as.numeric(IV_clump[[pcol]]))
  data.frame(
    Exposure         = exp,
    GW_Sig_Variants  = sum(pvals < p_thresh, na.rm = TRUE),
    Independent_IVs  = total_IVs,
    Own_Selected_IVs = sum(pvals_clumped < p_thresh, na.rm = TRUE),
    Conditional_F    = if (!is.null(cond_F_named) && exp %in% names(cond_F_named))
                         cond_F_named[[exp]] else NA_real_,
    stringsAsFactors = FALSE
  )
}))
## Independent IVs is constant across every row (it's the same LD-clumped
## panel for all exposures), so it is reported once in the title instead of
## as its own column. Own-Sel. IVs and its (%) are combined into a single
## column, but the count and percentage are formatted with independent
## fixed-width sub-fields *inside* that one cell, so both stay vertically
## aligned across rows despite the combined string's total length varying
## row to row (e.g. "7 (0.5%)" vs "424 (27.5%)").
top_bot <- function(l, m, r)
  paste0("  ", l,
         strrep("\u2500", 14), m,
         strrep("\u2500", 23), m,
         strrep("\u2500", 21), m,
         strrep("\u2500", 16), r)

message("\n", top_bot("\u250c", "\u252c", "\u2510"))
message(sprintf("  \u2502%-77s\u2502",
                sprintf("  IV SUMMARY TABLE (per exposure, after LD clumping): %s independent IVs",
                        format(total_IVs, big.mark = ","))))
message(top_bot("\u251c", "\u252c", "\u2524"))
message(sprintf("  \u2502 %-12s \u2502 %-21s \u2502 %-19s \u2502 %-14s \u2502",
                "Exposure", "GW-Sig. Variants", "Own-Sel. IVs (%)", "Conditional F"))
gencov_subtitle <- if (isTRUE(gencov_ldsc_partial)) {
  "(gencov: LDSC, partial)"
} else if (isTRUE(gencov_is_ldsc)) {
  "(gencov: LDSC)"
} else {
  "(gencov: 0)"
}
## Own-Sel. IVs (%) gets no before/after-clumping subtitle -- unlike
## GW-Sig. Variants, there's no "before clumping" variant of it to
## disambiguate from.
message(sprintf("  \u2502 %-12s \u2502 %-21s \u2502 %-19s \u2502 %-14s \u2502",
                "", "(before clumping)", "", gencov_subtitle))
message(top_bot("\u251c", "\u252c", "\u2524"))
for (i in seq_len(nrow(iv_summary))) {
  f_disp <- if (is.na(iv_summary$Conditional_F[i])) "NA"
            else sprintf("%.1f", iv_summary$Conditional_F[i])
  own_pct <- 100 * iv_summary$Own_Selected_IVs[i] / total_IVs
  own_cell <- sprintf("%7s (%5s%%)",
                      format(iv_summary$Own_Selected_IVs[i], big.mark = ","),
                      sprintf("%.1f", own_pct))
  message(sprintf("  \u2502 %-12s \u2502 %21s \u2502 %19s \u2502 %14s \u2502",
                  iv_summary$Exposure[i],
                  format(iv_summary$GW_Sig_Variants[i], big.mark = ","),
                  own_cell,
                  f_disp))
}
message(top_bot("\u2514", "\u2534", "\u2518"))

save(iv_summary, file = output_iv_summary_rdata)
message("  [OK] IV summary RData saved: ", output_iv_summary_rdata)

elapsed <- proc.time() - step_time
message(sprintf("  [DONE] Step 7 completed in %.1f seconds.", elapsed["elapsed"]))

step_time <- proc.time()
message("\n", strrep("=", 70))
message("[Step 8] Automatic matrix extraction...")
message(strrep("=", 70))

# Verify all expected columns are present
missing_betaX <- setdiff(beta_X_cols, names(IV_clump))
missing_seY   <- setdiff(se_Y_cols,   names(IV_clump))
missing_betaY <- setdiff(beta_Y_cols, names(IV_clump))
all_missing   <- c(missing_betaX, missing_seY, missing_betaY)

if (length(all_missing) > 0) {
  message("  [WARNING] Missing columns (trait may have been skipped):")
  for (m in all_missing) message("    - ", m)
  # Remove missing traits from extraction vectors
  beta_X_cols <- beta_X_cols[beta_X_cols %in% names(IV_clump)]
  se_Y_cols   <- se_Y_cols[se_Y_cols     %in% names(IV_clump)]
  beta_Y_cols <- beta_Y_cols[beta_Y_cols %in% names(IV_clump)]
}

message("  [OK] betaX cols: ", paste(standardise_trait_names(beta_X_cols), collapse=", "))
message("  [OK] seY   cols: ", paste(standardise_trait_names(se_Y_cols),   collapse=", "))
message("  [OK] betaY cols: ", paste(standardise_trait_names(beta_Y_cols), collapse=", "))

rs <- IV_clump$rsid

betaX <- to_num_matrix(IV_clump, beta_X_cols)
seY   <- to_num_matrix(IV_clump, se_Y_cols)
betaY <- to_num_matrix(IV_clump, beta_Y_cols)

rownames(betaX) <- rs; rownames(seY) <- rs; rownames(betaY) <- rs

# Strip suffix to recover plain trait names for colnames
colnames(betaX) <- gsub("_beta$", "", beta_X_cols)
colnames(seY)   <- gsub("_se$",   "", se_Y_cols)
colnames(betaY) <- gsub("_beta$", "", beta_Y_cols)

message("  [OK] betaX : ", nrow(betaX), " x ", ncol(betaX),
        " (", paste(standardise_trait_values(colnames(betaX)), collapse=", "), ")")
message("  [OK] seY   : ", nrow(seY),   " x ", ncol(seY),
        " (", paste(standardise_trait_values(colnames(seY)),   collapse=", "), ")")
message("  [OK] betaY : ", nrow(betaY), " x ", ncol(betaY),
        " (", paste(standardise_trait_values(colnames(betaY)), collapse=", "), ")")

elapsed <- proc.time() - step_time
message(sprintf("  [DONE] Step 8 completed in %.1f seconds.", elapsed["elapsed"]))

# =============================================================================
# § 10 — IVW WEIGHTING AND FINAL SAVE
# =============================================================================

step_time <- proc.time()
message("\n", strrep("=", 70))
message("[Step 9] IVW weighting and saving final output...")
message(strrep("=", 70))

mean_se   <- apply(seY, MARGIN=1, FUN=mean, na.rm=TRUE)
betaY_ivw <- betaY / mean_se
betaX_ivw <- betaX / mean_se

rownames(betaX_ivw) <- rs
rownames(betaY_ivw) <- rs
message("  [OK] Row names set to rsid in betaX_ivw and betaY_ivw.")

message("  [OK] IVW weights computed (mean seY across ",
        ncol(seY), " outcomes).")
message("  [INFO] Weight range: [",
        round(min(mean_se), 5), ", ", round(max(mean_se), 5), "]")

data <- list(
  betaY_ivw = betaY_ivw,
  betaX_ivw = betaX_ivw,
  rs        = rs,
  exposures = colnames(betaX),
  outcomes  = colnames(betaY)
)

local({
  data <- data
  colnames(data$betaX_ivw) <- standardise_trait_values(colnames(data$betaX_ivw))
  colnames(data$betaY_ivw) <- standardise_trait_values(colnames(data$betaY_ivw))
  data$exposures <- standardise_trait_values(data$exposures)
  data$outcomes  <- standardise_trait_values(data$outcomes)
  save(data, file = output_matrix_rdata)
})
message("  [OK] data saved: ", output_matrix_rdata)

message("\n  [SUMMARY]")
message("    Independent instruments : ", length(rs))
message("    Exposures               : ", paste(standardise_trait_values(colnames(betaX)), collapse=", "))
message("    Outcomes                : ", paste(standardise_trait_values(colnames(betaY)), collapse=", "))

elapsed <- proc.time() - step_time
message(sprintf("  [DONE] Step 9 completed in %.1f seconds.", elapsed["elapsed"]))

# =============================================================================
# § 11 — MrDAG-READY OUTPUT (data + gencov, positionally aligned)
# =============================================================================
# The object MrDAG() actually consumes: `data` with outcomes first then
# exposures (MrDAG's documented column order) and a matching per-IV gencov
# list, built the same way as Article_18's build_gencov_for_roles() -- LDSC
# intercept x per-IV SE_i,r x SE_i,s off-diagonal, SE_i,r^2 diagonal (stored; not used by the current
# algorithms), with
# IVs missing an SE for a trait imputed with that trait's mean SE rather
# than dropped. SEs are rescaled by mean_se (§10) to match betaY_ivw/betaX_ivw,
# which are on that same rescaled units -- gencov must match the scale of
# the data it corrects (MrDAG.R's own gencov formula multiplies by se_j*se_k
# in the units of `data`, not the raw GWAS-file SE).

step_time <- proc.time()
message("\n", strrep("=", 70))
message("[Step 10] Building MrDAG-ready data + gencov...")
message(strrep("=", 70))


# beta_*_cols and se_*_cols are filtered independently for presence in §8-9,
# so betaX/betaY's surviving columns are not guaranteed to match seX/seY's --
# reconcile by name, warn and drop anything without both, rather than assume
# positional alignment.
se_X_cols_final <- paste0(tolower(colnames(betaX)), "_se")
have_seX <- se_X_cols_final %in% names(IV_clump)
if (!all(have_seX))
  message("  [WARNING] No SE column for: ", paste(standardise_trait_values(colnames(betaX)[!have_seX]), collapse = ", "),
          " -- excluded from the MrDAG-ready output.")
mrdag_exposures <- colnames(betaX)[have_seX]
seX <- to_num_matrix(IV_clump, se_X_cols_final[have_seX])
rownames(seX) <- rs; colnames(seX) <- mrdag_exposures

mrdag_outcomes <- intersect(colnames(betaY), colnames(seY))
if (length(mrdag_outcomes) < ncol(betaY))
  message("  [WARNING] No SE column for: ",
          paste(standardise_trait_values(setdiff(colnames(betaY), mrdag_outcomes)), collapse = ", "),
          " -- excluded from the MrDAG-ready output.")

seX_ivw <- seX / mean_se
seY_ivw <- (seY / mean_se)[, mrdag_outcomes, drop = FALSE]
betaX_ivw_final <- betaX_ivw[, mrdag_exposures, drop = FALSE]
betaY_ivw_final <- betaY_ivw[, mrdag_outcomes, drop = FALSE]

gencov_traits <- c(mrdag_outcomes, mrdag_exposures)   # outcomes-then-exposures, matches MrDAG_data below
se_by_trait <- setNames(
  c(asplit(seY_ivw, 2), asplit(seX_ivw, 2)),
  gencov_traits
)

# gencov_traits are lowercase (matching betaX/betaY's column-naming
# convention); trait_sample_info/get_trait_path() expect the proper-case
# names declared in Exposures/Outcomes -- translate for anything that
# leaves this script's lowercase convention (trait_case_map is defined once,
# near the top of this file).
traits_proper <- unname(trait_case_map[gencov_traits])

mrdag_gencov <- tryCatch({
  # Reuses §8's LDSC result instead of re-running estimate_ldsc_gencov() (same GenomicSEM::ldsc()
  # regression, same trait_sample_info/munge_dir/ld_dir/hm3_snplist/annot_file): LDSC's cross-trait
  # intercepts are estimated pairwise, so a trait pair's I[.,.] entry doesn't depend on which other
  # traits happen to be batched into the same call -- §8's all_ldsc_traits (f_stat_traits + Outcomes)
  # is always a superset of gencov_traits here (same beta/se-presence filter for exposures; outcomes
  # only ever shrink via mrdag_outcomes's intersect()), so nothing needed here is missing from it.
  if (!exists("ldsc_result") || is.null(ldsc_result))
    stop("no LDSC result available from the F-statistic/LDSC step (§8) to reuse")

  missing_from_step8 <- setdiff(traits_proper, ldsc_result$trait_names)
  if (length(missing_from_step8) > 0)
    stop("§8's LDSC result does not cover: ", paste(missing_from_step8, collapse = ", "))

  ldsc_mrdag <- ldsc_result

  se_trait_mean <- setNames(
    vapply(gencov_traits, function(tr) mean(se_by_trait[[tr]], na.rm = TRUE), numeric(1)),
    gencov_traits
  )

  for (tr in gencov_traits) {
    missing_tr <- is.na(se_by_trait[[tr]])
    if (any(missing_tr)) {
      fallback <- se_trait_mean[[tr]]
      message("  [IMPUTE] ", sum(missing_tr), " of ", length(rs),
              " IV(s) missing an SE for ", standardise_trait_values(tr), "; imputed with its mean SE (",
              .display_round(fallback, 3), ").")
      se_by_trait[[tr]][missing_tr] <- fallback
    }
  }

  m <- length(gencov_traits)
  gencov_list <- vector("list", length(rs))
  for (i in seq_along(rs)) {
    Xi <- matrix(0, m, m, dimnames = list(gencov_traits, gencov_traits))
    for (a in gencov_traits) for (b in gencov_traits) {
      if (a == b) {
        Xi[a, a] <- se_by_trait[[a]][i]^2   # main diagonal always stored; the algorithms (MrDAG, MR2x) decide not to use it
      } else if (match(a, gencov_traits) < match(b, gencov_traits)) {
        # ldsc_mrdag$I, not full_matrix: per gencov_derivation.pdf Sec. 6(b),
        # the raw bivariate cross-trait LDSC intercept already equals
        # N_o,rs*rho_rs/sqrt(N_r*N_s) directly -- full_matrix instead divides
        # by the univariate (per-trait) intercepts, a different, unrelated
        # normalisation meant for building a bounded correlation matrix
        # elsewhere in this script, not for this per-IV Xi construction.
        gi <- ldsc_mrdag$I[trait_case_map[[a]], trait_case_map[[b]]]
        Xi[a, b] <- Xi[b, a] <- gi * se_by_trait[[a]][i] * se_by_trait[[b]][i]
      }
    }
    gencov_list[[i]] <- Xi
  }
  gencov_list
}, error = function(e) {
  ## Loud, not silent: a NULL gencov means MrDAG would later run WITHOUT any sample-overlap
  ## correction and nothing downstream would say so.
  message("\n", strrep("!", 70))
  message("  [WARNING] MrDAG-ready gencov could NOT be built: ", e$message)
  message("  [WARNING] MrDAG_input.RData will be saved with gencov = NULL, i.e. NO sample-overlap")
  message("  [WARNING] correction. Do not use it for the corrected analysis.")
  message(strrep("!", 70), "\n")
  warning("MrDAG-ready gencov could not be built (", e$message, "); saved with gencov = NULL.",
          call. = FALSE, immediate. = TRUE)
  NULL
})

q <- length(mrdag_outcomes); p <- length(mrdag_exposures)
MrDAGcheck <- list(Y_idx = seq_len(q), X_idx = seq_len(p) + q)
MrDAG_data <- cbind(betaY_ivw_final, betaX_ivw_final)
rownames(MrDAG_data) <- rs
colnames(MrDAG_data) <- standardise_trait_values(colnames(MrDAG_data))   # proper case, matching other saved files

gencov <- mrdag_gencov
if (!is.null(gencov))
  gencov <- lapply(gencov, function(Xi) {
    proper <- standardise_trait_values(rownames(Xi))
    dimnames(Xi) <- list(proper, proper)
    Xi
  })

# ── gencov correction diagnostics ───────────────────────────────────────────
# Magnitude of the sample-overlap (off-diagonal) correction relative to the raw sufficient statistic. The main
# diagonal (each trait's own measurement error variance) is stored in gencov but is not used by the current algorithms.
gencovDiagnostics <- NULL
if (!is.null(gencov)) {
  Xc            <- scale(MrDAG_data, center = TRUE, scale = FALSE)
  tXX_no_gencov <- crossprod(Xc)
  correction    <- Reduce(`+`, lapply(gencov, function(Xi) { diag(Xi) <- 0; Xi }))
  frob_no_gencov <- norm(tXX_no_gencov, "F")
  overall <- list(
    frobenius_no_gencov   = frob_no_gencov,
    frobenius_with_gencov = norm(tXX_no_gencov - correction, "F"),
    frobenius_correction  = norm(correction, "F"),
    pct_of_signal         = 100 * norm(correction, "F") / frob_no_gencov
  )
  mean_gencov  <- Reduce(`+`, gencov) / length(gencov)
  se_diag      <- sqrt(diag(mean_gencov))
  mean_gencov_correlation <- mean_gencov / outer(se_diag, se_diag)
  diag(mean_gencov_correlation) <- 1
  dimnames(mean_gencov_correlation) <- dimnames(mean_gencov)
  if (exists("digits", inherits = TRUE)) {
    overall <- lapply(overall, round, digits = digits)
    mean_gencov_correlation <- round(mean_gencov_correlation, digits = digits)
  }
  gencovDiagnostics <- list(overall = overall, mean_gencov_correlation = mean_gencov_correlation)
  message("")
  message("  [INFO] Sample-overlap (off-diagonal) correction: ", round(overall$pct_of_signal, 1),
          "% of the raw sufficient statistic's Frobenius norm.")
}
save(MrDAG_data, MrDAGcheck, gencov, rs, gencovDiagnostics, file = output_mrdag_input_rdata)
message("  [OK] MrDAG-ready input (data, MrDAGcheck, gencov, gencovDiagnostics) saved: ", output_mrdag_input_rdata)
message("  [INFO] Use directly as: MrDAG(data = MrDAG_data, MrDAGcheck = MrDAGcheck, gencov = gencov, ...)")

elapsed <- proc.time() - step_time
message(sprintf("  [DONE] Step 10 completed in %.1f seconds.", elapsed["elapsed"]))

message("\n", strrep("=", 70))
message("[COMPLETE] All steps finished successfully.")
message(strrep("=", 70), "\n")

gc()