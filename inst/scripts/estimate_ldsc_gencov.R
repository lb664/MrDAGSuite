
# =============================================================================
# estimate_ldsc_gencov.R -- companion to MrDAG_data_preprocessing
# Version 0.1.0
#
# LD Score Regression: sample-overlap intercept, genetic
# covariance/correlation (+ SE, p-value). Sourced by the main function (§1).
# =============================================================================

## .ldsc_with_intercept_se() -- a minimal fork of GenomicSEM::ldsc() (package version
## installed as of this writing) whose ONLY change is capturing `intercept.se`, which
## ldsc() already computes internally via the same block-jackknife used for every other
## LDSC SE it reports, but never returns -- only ever writes to its own log text. Every
## other line is verbatim from GenomicSEM::ldsc()'s body; changes are marked "## ADDED".
## Kept as a full fork rather than a small wrapper because intercept.se is a local
## variable computed deep inside ldsc()'s own loop, with no hook to extract it otherwise.
.ldsc_with_intercept_se <- function (traits, sample.prev, population.prev, ld, wld, trait.names = NULL,
    sep_weights = FALSE, chr = 22, n.blocks = 200, ldsc.log = NULL,
    stand = FALSE, select = FALSE, chisq.max = NA)
{
    time <- proc.time()
    begin.time <- Sys.time()
    if (is.null(ldsc.log)) {
        logtraits <- gsub(".*/", "", traits)
        log2 <- paste(logtraits, collapse = "_")
        if (object.size(log2) > 200) {
            log2 <- substr(log2, 1, 100)
        }
        log.file <- file(paste0(log2, "_ldsc.log"), open = "wt")
    }
    else {
        log.file <- file(paste0(ldsc.log, "_ldsc.log"), open = "wt")
    }
    if (!isFALSE(select) && !identical(select, "ODD") && !identical(select,
        "EVEN") && !is.numeric(select)) {
        stop("select must be one of the following values: FALSE, 'ODD', 'EVEN', or a numeric vector of chromosome numbers.")
    }
    GenomicSEM:::.LOG("Multivariate ld-score regression of ", length(traits),
        " traits ", "(", paste(traits, collapse = " "), ")",
        " began at: ", begin.time, file = log.file)
    if (identical(select, "ODD") || identical(select, "EVEN")) {
        odd <- seq(1, chr, 2)
        even <- seq(2, chr, 2)
    }
    n.traits <- length(traits)
    n.V <- n.traits * (n.traits + 1)/2
    if (n.traits > 18) {
        n.blocks <- (((n.traits + 1) * (n.traits + 2))/2) + 1
        GenomicSEM:::.LOG("     ", file = log.file, print = FALSE)
        GenomicSEM:::.LOG("Setting the number of blocks used to perform the block jacknife used to estimate the sampling covariance matrix (V) to ",
            n.blocks, file = log.file)
        GenomicSEM:::.LOG("This reflects the need to estimate V using at least one more block then their are nonredundant elements in the genetic covariance matrix that includes individual SNPs.",
            file = log.file)
        GenomicSEM:::.LOG("If the n.blocks is > 1000 you should carefully inspect output for any strange results, such as extremely significant Q_SNP estimates.",
            file = log.file)
        GenomicSEM:::.LOG("     ", file = log.file, print = FALSE)
        if (n.blocks > 1000) {
            warning("The number of blocks needed to estimate V is > 1000, which may result in sampling dependencies across the blocks used to estimate standard errors and can bias results.")
        }
    }
    if (!(is.null(trait.names))) {
        check_names <- stringr::str_detect(trait.names, "[+\\-*/^]")
        if (any(check_names))
            warning("Your trait names specified include mathematical arguments (e.g., +, -, *, /) that will be misread by lavaan. Please rename the traits using the trait.names argument.")
    }
    if (length(traits) == 1)
        warning("Our version of ldsc requires 2 or more traits. Please include an additional trait.")
    cov <- matrix(NA, nrow = n.traits, ncol = n.traits)
    V.hold <- matrix(NA, nrow = n.blocks, ncol = n.V)
    N.vec <- matrix(NA, nrow = 1, ncol = n.V)
    Liab.S <- rep(1, n.traits)
    I <- matrix(NA, nrow = n.traits, ncol = n.traits)
    I.se <- matrix(NA, nrow = n.traits, ncol = n.traits)   ## ADDED
    GenomicSEM:::.LOG("Reading in LD scores", file = log.file)
    if (isFALSE(select)) {
        x <- do.call("rbind", lapply(1:chr, function(i) {
            suppressMessages(readr::read_delim(file.path(ld, paste0(i,
                ".l2.ldscore.gz")), delim = "\t", escape_double = FALSE,
                trim_ws = TRUE, progress = FALSE))
        }))
    }
    if (identical(select, "ODD")) {
        x <- do.call("rbind", lapply(odd, function(i) {
            suppressMessages(readr::read_delim(file.path(ld, paste0(i,
                ".l2.ldscore.gz")), delim = "\t", escape_double = FALSE,
                trim_ws = TRUE, progress = FALSE))
        }))
    }
    if (identical(select, "EVEN")) {
        x <- do.call("rbind", lapply(even, function(i) {
            suppressMessages(readr::read_delim(file.path(ld, paste0(i,
                ".l2.ldscore.gz")), delim = "\t", escape_double = FALSE,
                trim_ws = TRUE, progress = FALSE))
        }))
    }
    if (is.numeric(select)) {
        x <- do.call("rbind", lapply(select, function(i) {
            suppressMessages(readr::read_delim(file.path(ld, paste0(i,
                ".l2.ldscore.gz")), delim = "\t", escape_double = FALSE,
                trim_ws = TRUE, progress = FALSE))
        }))
    }
    x$CM <- NULL
    x$MAF <- NULL
    if (sep_weights) {
        if (isFALSE(select)) {
            w <- do.call("rbind", lapply(1:chr, function(i) {
                suppressMessages(readr::read_delim(file.path(wld, paste0(i,
                  ".l2.ldscore.gz")), delim = "\t", escape_double = FALSE,
                  trim_ws = TRUE, progress = FALSE))
            }))
        }
        if (identical(select, "EVEN")) {
            w <- do.call("rbind", lapply(even, function(i) {
                suppressMessages(readr::read_delim(file.path(wld, paste0(i,
                  ".l2.ldscore.gz")), delim = "\t", escape_double = FALSE,
                  trim_ws = TRUE, progress = FALSE))
            }))
        }
        if (identical(select, "ODD")) {
            w <- do.call("rbind", lapply(odd, function(i) {
                suppressMessages(readr::read_delim(file.path(wld, paste0(i,
                  ".l2.ldscore.gz")), delim = "\t", escape_double = FALSE,
                  trim_ws = TRUE, progress = FALSE))
            }))
        }
        if (is.numeric(select)) {
            w <- do.call("rbind", lapply(select, function(i) {
                suppressMessages(readr::read_delim(file.path(wld, paste0(i,
                  ".l2.ldscore.gz")), delim = "\t", escape_double = FALSE,
                  trim_ws = TRUE, progress = FALSE))
            }))
        }
    }
    else {
        w <- x
    }
    w$CM <- NULL
    w$MAF <- NULL
    colnames(w)[ncol(w)] <- "wLD"
    if (isFALSE(select)) {
        m <- do.call("rbind", lapply(1:chr, function(i) {
            suppressMessages(readr::read_csv(file.path(ld, paste0(i,
                ".l2.M_5_50")), col_names = FALSE))
        }))
    }
    if (identical(select, "EVEN")) {
        m <- do.call("rbind", lapply(even, function(i) {
            suppressMessages(readr::read_csv(file.path(ld, paste0(i,
                ".l2.M_5_50")), col_names = FALSE))
        }))
    }
    if (identical(select, "ODD")) {
        m <- do.call("rbind", lapply(odd, function(i) {
            suppressMessages(readr::read_csv(file.path(ld, paste0(i,
                ".l2.M_5_50")), col_names = FALSE))
        }))
    }
    if (is.numeric(select)) {
        m <- do.call("rbind", lapply(select, function(i) {
            suppressMessages(readr::read_csv(file.path(ld, paste0(i,
                ".l2.M_5_50")), col_names = FALSE))
        }))
    }
    M.tot <- sum(m)
    m <- M.tot
    s <- 0
    all_y <- lapply(traits, function(chi1) {
        y1 <- suppressMessages(na.omit(readr::read_delim(chi1, delim = "\t",
            escape_double = FALSE, trim_ws = TRUE, progress = FALSE)))
        GenomicSEM:::.LOG("Read in summary statistics [", s <<- s + 1, "/",
            n.traits, "] from: ", chi1, file = log.file)
        merged <- merge(y1[, c("SNP", "N", "Z", "A1")], w[, c("SNP",
            "wLD")], by = "SNP", sort = FALSE)
        merged <- merge(merged, x, by = "SNP", sort = FALSE)
        merged <- merged[with(merged, order(CHR, BP)), ]
        GenomicSEM:::.LOG("Out of ", nrow(y1), " SNPs, ", nrow(merged), " remain after merging with LD-score files",
            file = log.file)
        if (is.na(chisq.max)) {
            chisq.max <- max(0.001 * max(merged$N), 80)
        }
        rm <- (merged$Z^2 > chisq.max)
        merged <- merged[!rm, ]
        GenomicSEM:::.LOG("Removing ", sum(rm), " SNPs with Chi^2 > ", chisq.max,
            "; ", nrow(merged), " remain", file = log.file)
        merged
    })
    s <- 1
    for (j in 1:n.traits) {
        chi1 <- traits[j]
        y1 <- all_y[[j]]
        y1$chi1 <- y1$Z^2
        for (k in j:length(traits)) {
            if (j == k) {
                GenomicSEM:::.LOG("     ", "     ", file = log.file, print = FALSE)
                GenomicSEM:::.LOG("Estimating heritability [", s, "/", n.V,
                  "] for: ", chi1, file = log.file)
                samp.prev <- sample.prev[j]
                pop.prev <- population.prev[j]
                merged <- y1
                n.snps <- nrow(merged)
                merged$intercept <- 1
                merged$x.tot <- merged$L2
                merged$x.tot.intercept <- 1
                tot.agg <- (M.tot * (mean(merged$chi1) - 1))/mean(merged$L2 *
                  merged$N)
                tot.agg <- max(tot.agg, 0)
                tot.agg <- min(tot.agg, 1)
                merged$ld <- pmax(merged$L2, 1)
                merged$w.ld <- pmax(merged$wLD, 1)
                merged$c <- tot.agg * merged$N/M.tot
                merged$het.w <- 1/(2 * (1 + (merged$c * merged$ld))^2)
                merged$oc.w <- 1/merged$w.ld
                merged$w <- merged$het.w * merged$oc.w
                merged$initial.w <- sqrt(merged$w)
                merged$weights <- merged$initial.w/sum(merged$initial.w)
                N.bar <- mean(merged$N)
                weighted.LD <- as.matrix(cbind(merged$L2, merged$intercept) *
                  merged$weights)
                weighted.chi <- as.matrix(merged$chi1 * merged$weights)
                n.annot <- 1
                select.from <- floor(seq(from = 1, to = n.snps,
                  length.out = (n.blocks + 1)))
                select.to <- c(select.from[2:n.blocks] - 1, n.snps)
                xty.block.values <- matrix(data = NA, nrow = n.blocks,
                  ncol = (n.annot + 1))
                xtx.block.values <- matrix(data = NA, nrow = ((n.annot +
                  1) * n.blocks), ncol = (n.annot + 1))
                colnames(xty.block.values) <- colnames(xtx.block.values) <- colnames(weighted.LD)
                replace.from <- seq(from = 1, to = nrow(xtx.block.values),
                  by = (n.annot + 1))
                replace.to <- seq(from = (n.annot + 1), to = nrow(xtx.block.values),
                  by = (n.annot + 1))
                for (i in 1:n.blocks) {
                  xty.block.values[i, ] <- t(t(weighted.LD[select.from[i]:select.to[i],
                    ]) %*% weighted.chi[select.from[i]:select.to[i],
                    ])
                  xtx.block.values[replace.from[i]:replace.to[i],
                    ] <- as.matrix(t(weighted.LD[select.from[i]:select.to[i],
                    ]) %*% weighted.LD[select.from[i]:select.to[i],
                    ])
                }
                xty <- as.matrix(colSums(xty.block.values))
                xtx <- matrix(data = NA, nrow = (n.annot + 1),
                  ncol = (n.annot + 1))
                colnames(xtx) <- colnames(weighted.LD)
                for (i in 1:nrow(xtx)) {
                  xtx[i, ] <- t(colSums(xtx.block.values[seq(from = i,
                    to = nrow(xtx.block.values), by = ncol(weighted.LD)),
                    ]))
                }
                reg <- solve(xtx, xty)
                intercept <- reg[2]
                coefs <- reg[1]/N.bar
                reg.tot <- coefs * m
                delete.from <- seq(from = 1, to = nrow(xtx.block.values),
                  by = ncol(xtx.block.values))
                delete.to <- seq(from = ncol(xtx.block.values),
                  to = nrow(xtx.block.values), by = ncol(xtx.block.values))
                delete.values <- matrix(data = NA, nrow = n.blocks,
                  ncol = (n.annot + 1))
                colnames(delete.values) <- colnames(weighted.LD)
                for (i in 1:n.blocks) {
                  xty.delete <- xty - xty.block.values[i, ]
                  xtx.delete <- xtx - xtx.block.values[delete.from[i]:delete.to[i],
                    ]
                  delete.values[i, ] <- solve(xtx.delete, xty.delete)
                }
                tot.delete.values <- delete.values[, 1:n.annot]
                pseudo.values <- matrix(data = NA, nrow = n.blocks,
                  ncol = length(reg))
                colnames(pseudo.values) <- colnames(weighted.LD)
                for (i in 1:n.blocks) {
                  pseudo.values[i, ] <- (n.blocks * reg) - ((n.blocks -
                    1) * delete.values[i, ])
                }
                jackknife.cov <- cov(pseudo.values)/n.blocks
                jackknife.se <- sqrt(diag(jackknife.cov))
                intercept.se <- jackknife.se[length(jackknife.se)]
                coef.cov <- jackknife.cov[1:n.annot, 1:n.annot]/(N.bar^2)
                cat.cov <- coef.cov * (m %*% t(m))
                tot.cov <- sum(cat.cov)
                tot.se <- sqrt(tot.cov)
                V.hold[, s] <- pseudo.values[, 1]
                N.vec[1, s] <- N.bar
                if (is.na(pop.prev) == F & is.na(samp.prev) ==
                  F) {
                  conversion.factor <- (pop.prev^2 * (1 - pop.prev)^2)/(samp.prev *
                    (1 - samp.prev) * dnorm(qnorm(1 - pop.prev))^2)
                  Liab.S[j] <- conversion.factor
                  GenomicSEM:::.LOG("     ", file = log.file, print = FALSE)
                  GenomicSEM:::.LOG("Please note that the results initially printed to the screen and log file reflect the NON-liability h2 and cov_g. However, a liability conversion is being used for trait ",
                    chi1, " when creating the genetic covariance matrix used as input for Genomic SEM and liability scale results are printed at the end of the log file.",
                    file = log.file)
                  GenomicSEM:::.LOG("     ", file = log.file, print = FALSE)
                }
                cov[j, j] <- reg.tot
                I[j, j] <- intercept
                I.se[j, j] <- intercept.se   ## ADDED
                lambda.gc <- median(merged$chi1)/qchisq(0.5,
                  df = 1)
                mean.Chi <- mean(merged$chi1)
                ratio <- (intercept - 1)/(mean.Chi - 1)
                ratio.se <- intercept.se/(mean.Chi - 1)
                GenomicSEM:::.LOG("Heritability Results for trait: ", chi1,
                  file = log.file)
                GenomicSEM:::.LOG("Mean Chi^2 across remaining SNPs: ", round(mean.Chi,
                  4), file = log.file)
                GenomicSEM:::.LOG("Lambda GC: ", round(lambda.gc, 4), file = log.file)
                GenomicSEM:::.LOG("Intercept: ", round(intercept, 4), " (",
                  round(intercept.se, 4), ")", file = log.file)
                GenomicSEM:::.LOG("Ratio: ", round(ratio, 4), " (", round(ratio.se,
                  4), ")", file = log.file)
                GenomicSEM:::.LOG("Total Observed Scale h2: ", round(reg.tot,
                  4), " (", round(tot.se, 4), ")", file = log.file)
                GenomicSEM:::.LOG("h2 Z: ", format(reg.tot/tot.se, digits = 3),
                  file = log.file)
            }
            if (j != k) {
                GenomicSEM:::.LOG("     ", file = log.file, print = FALSE)
                chi2 <- traits[k]
                GenomicSEM:::.LOG("Calculating genetic covariance [", s, "/",
                  n.V, "] for traits: ", chi1, " and ", chi2,
                  file = log.file)
                y2 <- all_y[[k]]
                y <- merge(y1, y2[, c("SNP", "N", "Z", "A1")],
                  by = "SNP", sort = FALSE)
                y$Z.x <- ifelse(y$A1.y == y$A1.x, y$Z.x, -y$Z.x)
                y$ZZ <- y$Z.y * y$Z.x
                y$chi2 <- y$Z.y^2
                merged <- na.omit(y)
                n.snps <- nrow(merged)
                GenomicSEM:::.LOG(n.snps, " SNPs remain after merging ", chi1,
                  " and ", chi2, " summary statistics", file = log.file)
                merged$intercept <- 1
                merged$x.tot <- merged$L2
                merged$x.tot.intercept <- 1
                tot.agg <- (M.tot * (mean(merged$chi1) - 1))/mean(merged$L2 *
                  merged$N.x)
                tot.agg <- max(tot.agg, 0)
                tot.agg <- min(tot.agg, 1)
                merged$ld <- pmax(merged$L2, 1)
                merged$w.ld <- pmax(merged$wLD, 1)
                merged$c <- tot.agg * merged$N.x/M.tot
                merged$het.w <- 1/(2 * (1 + (merged$c * merged$ld))^2)
                merged$oc.w <- 1/merged$w.ld
                merged$w <- merged$het.w * merged$oc.w
                merged$initial.w <- sqrt(merged$w)
                tot.agg2 <- (M.tot * (mean(merged$chi2) - 1))/mean(merged$L2 *
                  merged$N.y)
                tot.agg2 <- max(tot.agg2, 0)
                tot.agg2 <- min(tot.agg2, 1)
                merged$ld2 <- pmax(merged$L2, 1)
                merged$w.ld2 <- pmax(merged$wLD, 1)
                merged$c2 <- tot.agg2 * merged$N.y/M.tot
                merged$het.w2 <- 1/(2 * (1 + (merged$c2 * merged$ld))^2)
                merged$oc.w2 <- 1/merged$w.ld2
                merged$w2 <- merged$het.w2 * merged$oc.w2
                merged$initial.w2 <- sqrt(merged$w2)
                merged$weights_cov <- (merged$initial.w + merged$initial.w2)/sum(merged$initial.w +
                  merged$initial.w2)
                N.bar <- sqrt(mean(merged$N.x) * mean(merged$N.y))
                weighted.LD <- as.matrix(cbind(merged$L2, merged$intercept) *
                  merged$weights)
                weighted.chi <- as.matrix(merged$ZZ * merged$weights_cov)
                n.annot <- 1
                select.from <- floor(seq(from = 1, to = n.snps,
                  length.out = (n.blocks + 1)))
                select.to <- c(select.from[2:n.blocks] - 1, n.snps)
                xty.block.values <- matrix(data = NA, nrow = n.blocks,
                  ncol = (n.annot + 1))
                xtx.block.values <- matrix(data = NA, nrow = ((n.annot +
                  1) * n.blocks), ncol = (n.annot + 1))
                colnames(xty.block.values) <- colnames(xtx.block.values) <- colnames(weighted.LD)
                replace.from <- seq(from = 1, to = nrow(xtx.block.values),
                  by = (n.annot + 1))
                replace.to <- seq(from = (n.annot + 1), to = nrow(xtx.block.values),
                  by = (n.annot + 1))
                for (i in 1:n.blocks) {
                  xty.block.values[i, ] <- t(t(weighted.LD[select.from[i]:select.to[i],
                    ]) %*% weighted.chi[select.from[i]:select.to[i],
                    ])
                  xtx.block.values[replace.from[i]:replace.to[i],
                    ] <- as.matrix(t(weighted.LD[select.from[i]:select.to[i],
                    ]) %*% weighted.LD[select.from[i]:select.to[i],
                    ])
                }
                xty <- as.matrix(colSums(xty.block.values))
                xtx <- matrix(data = NA, nrow = (n.annot + 1),
                  ncol = (n.annot + 1))
                colnames(xtx) <- colnames(weighted.LD)
                for (i in 1:nrow(xtx)) {
                  xtx[i, ] <- t(colSums(xtx.block.values[seq(from = i,
                    to = nrow(xtx.block.values), by = ncol(weighted.LD)),
                    ]))
                }
                reg <- solve(xtx, xty)
                intercept <- reg[2]
                coefs <- reg[1]/N.bar
                reg.tot <- coefs * m
                delete.from <- seq(from = 1, to = nrow(xtx.block.values),
                  by = ncol(xtx.block.values))
                delete.to <- seq(from = ncol(xtx.block.values),
                  to = nrow(xtx.block.values), by = ncol(xtx.block.values))
                delete.values <- matrix(data = NA, nrow = n.blocks,
                  ncol = (n.annot + 1))
                colnames(delete.values) <- colnames(weighted.LD)
                for (i in 1:n.blocks) {
                  xty.delete <- xty - xty.block.values[i, ]
                  xtx.delete <- xtx - xtx.block.values[delete.from[i]:delete.to[i],
                    ]
                  delete.values[i, ] <- solve(xtx.delete, xty.delete)
                }
                tot.delete.values <- delete.values[, 1:n.annot]
                pseudo.values <- matrix(data = NA, nrow = n.blocks,
                  ncol = length(reg))
                colnames(pseudo.values) <- colnames(weighted.LD)
                for (i in 1:n.blocks) {
                  pseudo.values[i, ] <- (n.blocks * reg) - ((n.blocks -
                    1) * delete.values[i, ])
                }
                jackknife.cov <- cov(pseudo.values)/n.blocks
                jackknife.se <- sqrt(diag(jackknife.cov))
                intercept.se <- jackknife.se[length(jackknife.se)]
                coef.cov <- jackknife.cov[1:n.annot, 1:n.annot]/(N.bar^2)
                cat.cov <- coef.cov * (m %*% t(m))
                tot.cov <- sum(cat.cov)
                tot.se <- sqrt(tot.cov)
                V.hold[, s] <- pseudo.values[, 1]
                N.vec[1, s] <- N.bar
                cov[k, j] <- cov[j, k] <- reg.tot
                I[k, j] <- I[j, k] <- intercept
                I.se[k, j] <- I.se[j, k] <- intercept.se   ## ADDED
                GenomicSEM:::.LOG("Results for genetic covariance between: ",
                  chi1, " and ", chi2, file = log.file)
                GenomicSEM:::.LOG("Mean Z*Z: ", round(mean(merged$ZZ), 4),
                  file = log.file)
                GenomicSEM:::.LOG("Cross trait Intercept: ", round(intercept,
                  4), " (", round(intercept.se, 4), ")", file = log.file)
                GenomicSEM:::.LOG("Total Observed Scale Genetic Covariance (g_cov): ",
                  round(reg.tot, 4), " (", round(tot.se, 4),
                  ")", file = log.file)
                GenomicSEM:::.LOG("g_cov Z: ", format(reg.tot/tot.se, digits = 3),
                  file = log.file)
                GenomicSEM:::.LOG("g_cov P-value: ", format(2 * pnorm(abs(reg.tot/tot.se),
                  lower.tail = FALSE), digits = 5), file = log.file)
            }
            s <- s + 1
        }
    }
    v.out <- cov(V.hold)/crossprod(N.vec * (sqrt(n.blocks)/m))
    ratio <- tcrossprod(sqrt(Liab.S))
    S <- cov * ratio
    scaleO <- gdata::lowerTriangle(ratio, diag = TRUE)
    V <- v.out * tcrossprod(scaleO)
    colnames(S) <- if (is.null(trait.names))
        paste0("V", 1:ncol(S))
    else trait.names
    if (mean(Liab.S) != 1) {
        r <- nrow(S)
        SE <- matrix(0, r, r)
        SE[lower.tri(SE, diag = TRUE)] <- sqrt(diag(V))
        GenomicSEM:::.LOG(c("     ", "     "), file = log.file, print = FALSE)
        GenomicSEM:::.LOG("Liability Scale Results", file = log.file)
        for (j in 1:n.traits) {
            if (is.null(trait.names)) {
                chi1 <- traits[j]
            }
            else {
                chi1 <- trait.names[j]
            }
            for (k in j:length(traits)) {
                if (j == k) {
                  GenomicSEM:::.LOG("     ", file = log.file, print = FALSE)
                  GenomicSEM:::.LOG("Liability scale results for: ", chi1,
                    file = log.file)
                  GenomicSEM:::.LOG("Total Liability Scale h2: ", round(S[j,
                    j], 4), " (", round(SE[j, j], 4), ")", file = log.file)
                }
                if (j != k) {
                  if (is.null(trait.names)) {
                    chi2 <- traits[k]
                  }
                  else {
                    chi2 <- trait.names[k]
                  }
                  GenomicSEM:::.LOG("Total Liability Scale Genetic Covariance between ",
                    chi1, " and ", chi2, ": ", round(S[k, j],
                      4), " (", round(SE[k, j], 4), ")", file = log.file)
                  GenomicSEM:::.LOG("     ", file = log.file, print = FALSE)
                }
            }
        }
    }
    if (all(diag(S) > 0)) {
        ratio <- tcrossprod(1/sqrt(diag(S)))
        S_Stand <- S * ratio
        scaleO <- gdata::lowerTriangle(ratio, diag = TRUE)
        V_Stand <- V * tcrossprod(scaleO)
        r <- nrow(S)
        SE_Stand <- matrix(0, r, r)
        SE_Stand[lower.tri(SE_Stand, diag = TRUE)] <- sqrt(diag(V_Stand))
        GenomicSEM:::.LOG(c("     ", "     "), file = log.file, print = FALSE)
        GenomicSEM:::.LOG("Genetic Correlation Results", file = log.file)
        for (j in 1:n.traits) {
            if (is.null(trait.names)) {
                chi1 <- traits[j]
            }
            else {
                chi1 <- trait.names[j]
            }
            for (k in j:length(traits)) {
                if (j != k) {
                  if (is.null(trait.names)) {
                    chi2 <- traits[k]
                  }
                  else {
                    chi2 <- trait.names[k]
                  }
                  GenomicSEM:::.LOG("Genetic Correlation between ", chi1,
                    " and ", chi2, ": ", round(S_Stand[k, j],
                      4), " (", round(SE_Stand[k, j], 4), ")",
                    file = log.file)
                  GenomicSEM:::.LOG("     ", file = log.file, print = FALSE)
                }
            }
        }
    }
    else {
        warning("Your genetic covariance matrix includes traits estimated to have a negative heritability.")
        GenomicSEM:::.LOG("Your genetic covariance matrix includes traits estimated to have a negative heritability.",
            file = log.file, print = FALSE)
        GenomicSEM:::.LOG("Genetic correlation results could not be computed due to negative heritability estimates.",
            file = log.file)
    }
    end.time <- Sys.time()
    total.time <- difftime(time1 = end.time, time2 = begin.time,
        units = "sec")
    mins <- floor(floor(total.time)/60)
    secs <- floor(total.time - mins * 60)
    GenomicSEM:::.LOG("     ", file = log.file, print = FALSE)
    GenomicSEM:::.LOG("LDSC finished running at ", end.time, file = log.file)
    GenomicSEM:::.LOG("Running LDSC for all files took ", mins, " minutes and ",
        secs, " seconds", file = log.file)
    GenomicSEM:::.LOG("     ", file = log.file, print = FALSE)
    flush(log.file)
    close(log.file)
    if (stand) {
        list(V = V, S = S, I = I, I.se = I.se, N = N.vec, m = m, V_Stand = V_Stand,   ## ADDED I.se
            S_Stand = S_Stand)
    }
    else {
        list(V = V, S = S, I = I, I.se = I.se, N = N.vec, m = m)   ## ADDED I.se
    }
}

#' Estimate the sample-overlap correlation matrix between the traits via LDSC
estimate_ldsc_gencov <- function(traits, trait_paths,
                                  sample_info, munge_dir, ld_dir, hm3_snplist,
                                  annot_file = NULL) {

  if (!requireNamespace("GenomicSEM", quietly = TRUE))
    stop("FATAL: GenomicSEM is required. Install via ",
         "remotes::install_github('GenomicSEM/GenomicSEM').")
  if (!dir.exists(ld_dir))
    stop("FATAL: LD score reference directory not found: ", ld_dir)
  if (!file.exists(hm3_snplist))
    stop("FATAL: HapMap3 SNP list not found: ", hm3_snplist)

  missing_info <- setdiff(traits, names(sample_info))
  if (length(missing_info) > 0)
    stop("FATAL: No trait_sample_info entry for: ",
         paste(missing_info, collapse = ", "))

  if (!dir.exists(munge_dir)) dir.create(munge_dir, recursive = TRUE)

  message("  [LDSC] Preparing to munge ", length(traits), " traits: ",
          paste(traits, collapse = ", "))

  files       <- character(length(traits))
  trait_names <- traits
  N_vec       <- numeric(length(traits))

  for (i in seq_along(traits)) {
    tr   <- traits[i]
    info <- sample_info[[tr]]

    files[i] <- trait_paths[[tr]]

    if (info$type == "quant") {
      N_vec[i] <- info$N
    } else if (info$type == "binary") {
      N_vec[i] <- info$cases + info$controls
    } else {
      stop("FATAL: Unknown trait type '", info$type, "' for ", tr,
           " (expected 'quant' or 'binary').")
    }
  }

  .build_manifest <- function(tr) {
    fi <- file.info(trait_paths[[tr]])
    list(
      file_size   = unname(fi$size),
      ## Numeric epoch, not as.character(fi$mtime): the character form
      ## renders in the CURRENT session's timezone, so if two runs sit on
      ## opposite sides of a DST transition (or otherwise have a different
      ## effective TZ) the same unchanged file prints a different mtime
      ## string and looks like a stale cache even though nothing changed.
      ## The numeric epoch is timezone-independent, so it doesn't have this
      ## failure mode. Confirmed on this exact pipeline: file_size and
      ## sample_info matched exactly for all 5 traits, but file_mtime as a
      ## string was off by exactly one hour on every single one -- a false
      ## invalidation, not a real input change.
      file_mtime  = as.numeric(fi$mtime),
      sample_info = sample_info[[tr]],
      ## Bump when the pre-rename/allele logic changes, so cached munged files
      ## built by an older version are rebuilt instead of silently reused.
      ldsc_input_version = 2L
    )
  }

  .manifest_path <- function(tr) file.path(munge_dir, paste0(tr, "_manifest.RData"))

  .write_manifest <- function(tr) {
    manifest <- .build_manifest(tr)
    save(manifest, file = .manifest_path(tr))
  }

  .read_manifest <- function(tr) {
    path <- .manifest_path(tr)
    if (!file.exists(path)) return(NULL)
    e <- new.env()
    tryCatch({
      load(path, envir = e)
      if (!exists("manifest", envir = e, inherits = FALSE)) return(NULL)
      get("manifest", envir = e)
    }, error = function(err) NULL)
  }

  .manifest_matches <- function(tr) {
    old <- .read_manifest(tr)
    if (is.null(old)) return(FALSE)
    isTRUE(identical(old, .build_manifest(tr)))
  }

  manifest_ok  <- vapply(traits, .manifest_matches, logical(1))
  has_sumstats <- file.exists(file.path(munge_dir, paste0(traits, ".sumstats.gz")))
  has_interim  <- file.exists(file.path(munge_dir, paste0(traits, "_ldsc_input.tsv.gz")))

  is_cached         <- manifest_ok & has_sumstats
  has_valid_interim <- manifest_ok & has_interim & !is_cached

  if (any(is_cached))
    message("  [LDSC] Already munged, skipping: ",
            paste(traits[is_cached], collapse = ", "))
  if (any(has_valid_interim))
    message("  [LDSC] Pre-rename already done, skipping to munge(): ",
            paste(traits[has_valid_interim], collapse = ", "))

  stale_sumstats <- traits[has_sumstats & !manifest_ok]
  stale_interim  <- traits[has_interim & !has_sumstats & !manifest_ok]
  stale          <- union(stale_sumstats, stale_interim)
  if (length(stale) > 0)
    message("  [LDSC] Cache invalidated (input changed since last run), ",
            "re-munging from scratch: ", paste(stale, collapse = ", "))

  need_fresh <- traits[!is_cached & !has_valid_interim]

  need_fresh_cols <- list()
  trait_data      <- list()
  really_fresh    <- character(0)
  for (tr in need_fresh) {
    dt_hdr <- tryCatch(read_trait_lazy(trait_paths[[tr]], tr), error = function(e) NULL)
    if (is.null(dt_hdr)) {
      message("  [LDSC-SKIP] ", tr, ": could not read file header -- ",
              "excluded from sample-overlap estimation, other traits unaffected.")
      next
    }
    cols_tr <- detect_columns(names(dt_hdr), tr)
    if (is.null(cols_tr)) {
      message("  [LDSC-SKIP] ", tr, ": column detection failed -- ",
              "excluded from sample-overlap estimation, other traits unaffected.")
      next
    }
    # Effect_A2 correction -- same rule as MrDAG_data_preprocessing.R: a file that
    # reports only Effect_A2 (no Effect_A1) gives the effect of allele A2, so A2 is
    # the effect allele. Without this swap the effect was attached to the wrong allele
    # and the Z sign was reversed for every SNP (AF, Nielsen et al. 2018).
    if (!is.na(cols_tr$beta_col) && cols_tr$beta_col == "Effect_A2" &&
        !"Effect_A1" %in% names(dt_hdr)) {
      tmp_ea         <- cols_tr$ea_col
      cols_tr$ea_col <- cols_tr$oa_col
      cols_tr$oa_col <- tmp_ea
      message("  [LDSC] ", tr, ": Effect_A2 without Effect_A1 -- treating A2 as effect allele.")
    }
    need_fresh_cols[[tr]] <- cols_tr
    trait_data[[tr]]      <- dt_hdr
    really_fresh <- c(really_fresh, tr)
  }
  need_fresh <- really_fresh

  message("  [LDSC] Pre-renaming columns to munge()'s expected format...")

  rsid_content_ok <- setNames(logical(length(need_fresh)), need_fresh)
  for (tr in need_fresh) {
    cols <- need_fresh_cols[[tr]]
    if (identical(cols$chrpos_method, "variant_field") || is.na(cols$rsid_col)) {
      rsid_content_ok[tr] <- NA  # not applicable -- handled by dedicated branches
      next
    }
    # Reuses read_trait_lazy()'s already-loaded data rather than re-reading
    # from disk, so column names always match what was detected.
    rsid_vals <- tryCatch(
      iconv(as.character(trait_data[[tr]][[cols$rsid_col]]), to = "UTF-8", sub = ""),
      error = function(e) NULL
    )
    if (is.null(rsid_vals)) {
      rsid_content_ok[tr] <- FALSE
    } else {
      sample_vals <- utils::head(rsid_vals[!is.na(rsid_vals)], 200)
      rsid_content_ok[tr] <- sum(grepl("^rs", sample_vals), na.rm = TRUE) > 0
    }
    if (!isTRUE(rsid_content_ok[tr]) && !is.na(rsid_content_ok[tr]))
      message("  [LDSC] ", tr, ": rsid column '", cols$rsid_col, "' does not ",
              "contain usable rsids (sample check failed) -- will use ",
              "chr:pos translation instead.")
  }

  needs_variant_lookup <- any(vapply(need_fresh, function(tr)
    identical(need_fresh_cols[[tr]]$chrpos_method, "variant_field"), logical(1)))
  needs_chrpos_lookup  <- any(vapply(need_fresh, function(tr) {
    cols <- need_fresh_cols[[tr]]
    !identical(cols$chrpos_method, "variant_field") &&
      (is.na(cols$rsid_col) || !isTRUE(rsid_content_ok[tr]))
  }, logical(1)))

  annot_lookup <- NULL
  if (needs_variant_lookup || needs_chrpos_lookup) {
    if (is.null(annot_file) || !file.exists(annot_file))
      stop("FATAL: at least one trait has no usable rsid of its own (UKB ",
           "variant-only format, no rsid column, or an rsid column whose ",
           "content failed validation), but no valid annot_file was ",
           "supplied to translate to rsid. Cannot proceed.")
    message("  [LDSC] Loading variant/chr:pos -> rsid lookup...")
    annot_raw <- vroom::vroom(annot_file,
                              col_select = c("variant", "chr", "pos", "rsid"),
                              altrep = FALSE,
                              show_col_types = FALSE, progress = FALSE)
    annot_lookup <- data.frame(
      variant    = annot_raw$variant,
      MarkerName = paste(annot_raw$chr, annot_raw$pos, sep = ":"),
      rsid       = annot_raw$rsid,
      stringsAsFactors = FALSE
    )
    rm(annot_raw); gc()
    message("  [LDSC] Lookup table: ",
            format(nrow(annot_lookup), big.mark = ","), " variants.")
  }

  .chrpos_raw_cols <- function(method) {
    switch(method,
      "hm_cols"         = c("hm_chrom", "hm_pos"),
      "base_cols"       = c("chromosome", "base_pair_location"),
      "Chr_Pos"         = c("Chr", "Pos"),
      "Chrom_Pos"       = c("Chrom", "Pos"),
      "chr_bp_hg19"     = c("chr", "bp_hg19"),
      "MarkerName_col"  = c("MarkerName"),
      "markername_col"  = c("markername"),
      "CHROM_POS"       = c("CHROM", "POS"),
      "CHR_BP"          = c("CHR", "BP"),
      "CHR_POS"         = c("CHR", "POS"),
      stop("Unrecognised chrpos_method for rsid translation: ", method)
    )
  }

  renamed_files <- setNames(rep(NA_character_, length(need_fresh)), need_fresh)

  for (tr in need_fresh) {
    cols <- need_fresh_cols[[tr]]

    result <- tryCatch({
      if (identical(cols$chrpos_method, "variant_field")) {
        if (any(is.na(c(cols$ea_col, cols$beta_col, cols$pval_col))))
          stop(tr, " (UKB format) is missing a required column (EA/beta/pval).")

        dt_raw <- trait_data[[tr]]
        dt <- data.frame(
          variant = as.character(dt_raw[["variant"]]),
          A1      = dt_raw[[cols$ea_col]],
          effect  = dt_raw[[cols$beta_col]],
          P       = dt_raw[[cols$pval_col]],
          stringsAsFactors = FALSE
        )

        rsid_vec <- annot_lookup$rsid[match(dt$variant, annot_lookup$variant)]
        parts    <- do.call(rbind, strsplit(as.character(dt$variant), ":", fixed = TRUE))
        ref_vec  <- if (ncol(parts) >= 3) parts[, 3] else NA_character_
        alt_vec  <- if (ncol(parts) >= 4) parts[, 4] else NA_character_
        # UKB (Neale lab) effects are for the ALT allele of variant = chr:pos:ref:alt, so
        # A1 must be alt and A2 ref. The detected ea_col (e.g. minor_allele) can be the
        # ref allele, which gave the wrong allele label (and so the wrong Z sign after
        # merging with another trait) for every variant whose minor allele is ref.
        a1_vec <- iconv(as.character(alt_vec), to = "UTF-8", sub = "")
        a2_vec <- iconv(as.character(ref_vec), to = "UTF-8", sub = "")

        dt <- data.frame(SNP = rsid_vec, A1 = a1_vec, A2 = a2_vec,
                         effect = dt$effect, P = dt$P, stringsAsFactors = FALSE)
        dt <- dt[!is.na(dt$SNP), , drop = FALSE]

      } else if (is.na(cols$rsid_col) || !isTRUE(rsid_content_ok[tr])) {
        if (is.na(cols$chrpos_method) || identical(cols$chrpos_method, "rsid_only"))
          stop(tr, ": no usable rsid, and no chr:pos fallback is available.")
        if (any(is.na(c(cols$ea_col, cols$oa_col, cols$beta_col, cols$pval_col))))
          stop(tr, " (chr:pos fallback) is missing a required column (EA/OA/beta/pval).")

        raw_cols <- .chrpos_raw_cols(cols$chrpos_method)
        dt_raw   <- trait_data[[tr]]
        missing_raw <- setdiff(c(raw_cols, cols$ea_col, cols$oa_col, cols$beta_col, cols$pval_col),
                                names(dt_raw))
        if (length(missing_raw) > 0)
          stop(tr, " (chr:pos fallback): column(s) not found in the loaded data: ",
               paste(missing_raw, collapse = ", "))

        marker_vec <- build_chrpos(as.data.table(dt_raw[, raw_cols, with = FALSE]), cols$chrpos_method)
        rsid_vec   <- annot_lookup$rsid[match(marker_vec, annot_lookup$MarkerName)]

        ea_vec <- iconv(as.character(dt_raw[[cols$ea_col]]), to = "UTF-8", sub = "")
        oa_vec <- iconv(as.character(dt_raw[[cols$oa_col]]), to = "UTF-8", sub = "")

        dt <- data.frame(SNP = rsid_vec, A1 = ea_vec, A2 = oa_vec,
                         effect = dt_raw[[cols$beta_col]], P = dt_raw[[cols$pval_col]],
                         stringsAsFactors = FALSE)
        dt <- dt[!is.na(dt$SNP), , drop = FALSE]

      } else {
        need <- c(cols$rsid_col, cols$ea_col, cols$oa_col, cols$beta_col, cols$pval_col)
        if (any(is.na(need)))
          stop(tr, " is missing a required column for LDSC (rsid/EA/OA/beta/pval).")

        dt_raw <- trait_data[[tr]]
        missing_raw <- setdiff(need, names(dt_raw))
        if (length(missing_raw) > 0)
          stop(tr, ": column(s) not found in the loaded data: ",
               paste(missing_raw, collapse = ", "))

        dt <- data.frame(
          SNP    = dt_raw[[cols$rsid_col]],
          A1     = dt_raw[[cols$ea_col]],
          A2     = dt_raw[[cols$oa_col]],
          effect = dt_raw[[cols$beta_col]],
          P      = dt_raw[[cols$pval_col]],
          stringsAsFactors = FALSE
        )

        dt$SNP <- iconv(as.character(dt$SNP), to = "UTF-8", sub = "")
        dt$A1  <- iconv(as.character(dt$A1),  to = "UTF-8", sub = "")
        dt$A2  <- iconv(as.character(dt$A2),  to = "UTF-8", sub = "")

        dt$SNP <- sub(".*?(rs[0-9]+).*", "\\1", dt$SNP, perl = TRUE)
        dt$SNP[!grepl("^rs[0-9]+$", dt$SNP)] <- NA_character_
        dt <- dt[!is.na(dt$SNP), , drop = FALSE]
      }

      out_path <- file.path(munge_dir, paste0(tr, "_ldsc_input.tsv.gz"))
      vroom::vroom_write(dt, out_path, delim = "\t")
      rm(dt); trait_data[[tr]] <- NULL; gc()
      out_path
    }, error = function(e) {
      message("  [LDSC-SKIP] ", tr, ": ", conditionMessage(e),
              " -- excluded from sample-overlap estimation, other traits unaffected.")
      trait_data[[tr]] <<- NULL
      NA_character_
    })

    renamed_files[tr] <- result
    if (!is.na(result)) .write_manifest(tr)
  }

  fresh_ok <- need_fresh[!is.na(renamed_files)]

  to_munge <- c(fresh_ok, traits[has_valid_interim])
  renamed_files_all <- setNames(
    c(unname(renamed_files[fresh_ok]),
      file.path(munge_dir, paste0(traits[has_valid_interim], "_ldsc_input.tsv.gz"))),
    to_munge
  )

  if (length(is_cached[is_cached]) + length(to_munge) < 2)
    stop("FATAL: fewer than 2 traits available for LDSC (",
         length(is_cached[is_cached]), " cached + ", length(to_munge),
         " pending munge()). Cannot estimate sample-overlap covariance.")

  skipped <- setdiff(need_fresh, fresh_ok)
  if (length(skipped) > 0)
    message("  [LDSC] Excluded (pre-rename failed): ", paste(skipped, collapse = ", "))

  munge_ok <- character(0)
  # Consecutive (not total) failures: isolated trait-specific problems (a
  # malformed file, a genuinely missing column) should still be skipped one
  # at a time as before, but several in a row is a strong signal of a
  # systemic problem (wrong munge_dir/hm3_snplist, permissions, disk space)
  # rather than bad data -- abort loudly instead of silently retrying the
  # same failure for every remaining trait, which previously ran for hours
  # across an entire trait list without ever surfacing as more than one
  # easy-to-miss [LDSC-SKIP] line per trait.
  consecutive_failures <- 0
  for (tr in to_munge) {
    munge_result <- tryCatch({
      invisible(capture.output(
        GenomicSEM::munge(
          files       = unname(renamed_files_all[tr]),
          hm3         = hm3_snplist,
          trait.names = tr,
          N           = N_vec[match(tr, traits)],
          info.filter = 0.9,
          maf.filter  = 0.01,
          log.name    = file.path(munge_dir, paste0(tr, "_munge"))
        )
      ))
      # munge() writes <trait_name>.sumstats.gz into the current working
      # directory; move/verify into munge_dir for a predictable location.
      # Always relocate (overwriting any stale copy already in munge_dir):
      # reaching this loop at all means `to_munge` already decided this
      # trait needed a fresh run, so a pre-existing same-named file in
      # munge_dir is exactly what's stale and must not be silently kept
      # over the file munge() just produced.
      sumstats_file <- paste0(tr, ".sumstats.gz")
      if (file.exists(sumstats_file))
        file.rename(sumstats_file, file.path(munge_dir, sumstats_file))
      if (!file.exists(file.path(munge_dir, sumstats_file)))
        stop("munge() did not produce the expected output file")
      TRUE
    }, error = function(e) {
      message("  [LDSC-SKIP] ", tr, ": munge() failed -- ", conditionMessage(e),
              " -- excluded from sample-overlap estimation, other traits unaffected. ",
              "Its pre-renamed intermediate file is kept, so a rerun will skip ",
              "straight back to munge() without repeating the pre-rename step.")
      structure(FALSE, error_message = conditionMessage(e))
    })

    leftover_logs <- list.files(munge_dir, pattern = paste0("^", tr, "_munge"),
                                full.names = TRUE)
    if (length(leftover_logs) > 0) unlink(leftover_logs)

    if (isTRUE(munge_result)) {
      munge_ok <- c(munge_ok, tr)
      consecutive_failures <- 0
    } else {
      consecutive_failures <- consecutive_failures + 1
      if (consecutive_failures >= 3)
        stop("3 traits in a row failed munge() (most recently '", tr, "': ",
             attr(munge_result, "error_message"), ") -- this looks systemic ",
             "(wrong munge_dir/hm3_snplist, permissions, disk space), not ",
             "trait-specific data problems. Aborting rather than silently ",
             "retrying the same failure for every remaining trait.")
    }
  }

  if (length(is_cached[is_cached]) + length(munge_ok) < 2)
    stop("FATAL: fewer than 2 traits successfully available for LDSC after ",
         "munge() (", length(is_cached[is_cached]), " cached + ",
         length(munge_ok), " freshly munged). Cannot estimate sample-overlap ",
         "covariance.")

  # Final trait set: cached + freshly-munged-and-successful, preserving
  # the original `traits` order for determinism.
  final_ok    <- traits %in% c(traits[is_cached], munge_ok)
  trait_names <- traits[final_ok]
  N_vec       <- N_vec[final_ok]
  munged_paths <- file.path(munge_dir, paste0(trait_names, ".sumstats.gz"))

  sample_prev <- vapply(trait_names, function(tr) {
    info <- sample_info[[tr]]
    if (info$type == "quant") NA_real_
    else info$cases / (info$cases + info$controls)
  }, numeric(1))

  population_prev <- vapply(trait_names, function(tr) {
    info <- sample_info[[tr]]
    if (info$type == "quant") NA_real_
    else info$population_prevalence
  }, numeric(1))

  # ── ldsc(): the actual multivariable LD Score Regression ─────────────────
  message("  [LDSC] Running LD Score Regression across ",
          length(trait_names), " traits (", choose(length(trait_names), 2),
          " pairs)...")
  ldsc_log_name <- "LDSC_all_traits"

  # ldsc() itself has no cat()/message() calls, but some dependency in its call chain writes a
  # \r-updating progress indicator straight to stderr (the conventional way to show progress
  # without polluting piped/redirected stdout) -- capture.output() alone only redirects stdout, so
  # that leaked through as a wall of blank lines in consoles that can't overwrite \r in place (e.g.
  # RStudio). Sinking both stdout and stderr to a discarded temp file catches it regardless of which
  # stream or mechanism it uses.
  .quiet_sink_file <- tempfile()
  .quiet_sink_con  <- file(.quiet_sink_file, open = "wt")
  sink(.quiet_sink_con, type = "output")
  sink(.quiet_sink_con, type = "message")
  ldsc_out <- tryCatch({
    ldsc_out <- .ldsc_with_intercept_se(
      traits          = munged_paths,
      sample.prev     = sample_prev,
      population.prev = population_prev,
      ld              = ld_dir,
      wld             = ld_dir,
      trait.names     = trait_names,
      ldsc.log        = file.path(munge_dir, ldsc_log_name)
    )
    sink(type = "message"); sink(type = "output"); close(.quiet_sink_con); unlink(.quiet_sink_file)
    ldsc_out
  }, error = function(e) {
    sink(type = "message"); sink(type = "output"); close(.quiet_sink_con); unlink(.quiet_sink_file)
    stop("FATAL: GenomicSEM::ldsc() failed: ", e$message)
  })

  ## GenomicSEM::ldsc() does not always close its own log-file connection.
  ## A before/after connection-count diff can miss it when this function is
  ## called more than once per session with the same log path (e.g. once
  ## for a conditional F-statistic and once for a per-IV gencov), so close
  ## any still-open connection to this specific log file directly, by name,
  ## rather than by when it was opened.
  all_conns <- showConnections(all = TRUE)
  log_leaks <- rownames(all_conns)[grepl(ldsc_log_name, all_conns[, "description"], fixed = TRUE)]
  for (cn in log_leaks) tryCatch(close(getConnection(as.integer(cn))), error = function(e) NULL)

  I_mat <- ldsc_out$I
  rownames(I_mat) <- trait_names
  colnames(I_mat) <- trait_names

  I_se_mat <- ldsc_out$I.se
  rownames(I_se_mat) <- trait_names
  colnames(I_se_mat) <- trait_names

  if (is.null(ldsc_out$S))
    stop("FATAL: GenomicSEM::ldsc() did not return an $S (genetic ",
         "covariance) element -- check the installed GenomicSEM version's ",
         "ldsc() output structure directly before trusting this function's ",
         "S/rg/rg_se/rg_p output.")

  S_mat <- ldsc_out$S
  rownames(S_mat) <- trait_names
  colnames(S_mat) <- trait_names

  h2_diag <- diag(S_mat)
  if (any(h2_diag <= 0)) {
    message("  [WARNING] Non-positive LDSC heritability estimate for: ",
            paste(trait_names[h2_diag <= 0], collapse = ", "),
            " -- genetic correlations involving these traits will be NaN ",
            "in the returned rg matrix.")
  }
  rg_mat <- S_mat / sqrt(outer(h2_diag, h2_diag))
  diag(rg_mat) <- 1
  dimnames(rg_mat) <- list(trait_names, trait_names)

  n_traits_ldsc <- length(trait_names)

  rg_se <- matrix(NA_real_, n_traits_ldsc, n_traits_ldsc,
                   dimnames = list(trait_names, trait_names))
  rg_p  <- matrix(NA_real_, n_traits_ldsc, n_traits_ldsc,
                   dimnames = list(trait_names, trait_names))

  rg_se_result <- tryCatch({
    if (is.null(ldsc_out$V))
      stop("ldsc_out$V not present")

    V_raw <- ldsc_out$V
    expected_len <- n_traits_ldsc * (n_traits_ldsc + 1) / 2
    if (nrow(V_raw) != expected_len || ncol(V_raw) != expected_len)
      stop("ldsc_out$V is ", nrow(V_raw), "x", ncol(V_raw), ", expected ",
           expected_len, "x", expected_len, " (k*(k+1)/2 for k = ",
           n_traits_ldsc, " traits) -- element order cannot be safely ",
           "assumed from this mismatch alone")

    vech_index <- function(i, j, k) {
      if (i < j) { tmp <- i; i <- j; j <- tmp }   # lower triangle: i >= j
      offset <- if (j == 1) 0 else sum(k - seq_len(j - 1) + 1)
      offset + (i - j + 1)
    }

    for (a in seq_len(n_traits_ldsc)) {
      for (b in seq_len(n_traits_ldsc)) {
        if (a == b) next   # diagonal left NA: SE of self-correlation is not meaningful
        idx_aa <- vech_index(a, a, n_traits_ldsc)
        idx_bb <- vech_index(b, b, n_traits_ldsc)
        idx_ab <- vech_index(a, b, n_traits_ldsc)
        Sigma_sub <- V_raw[c(idx_aa, idx_bb, idx_ab), c(idx_aa, idx_bb, idx_ab)]
        S_aa <- S_mat[a, a]; S_bb <- S_mat[b, b]
        rg_ab <- rg_mat[a, b]
        # Gradient of rg_ab = S_ab / sqrt(S_aa * S_bb) w.r.t. (S_aa, S_bb, S_ab)
        grad <- c(-0.5 * rg_ab / S_aa, -0.5 * rg_ab / S_bb, 1 / sqrt(S_aa * S_bb))
        var_rg <- as.numeric(t(grad) %*% Sigma_sub %*% grad)
        if (is.finite(var_rg) && var_rg >= 0) rg_se[a, b] <- sqrt(var_rg)
      }
    }
    rg_p_tmp <- 2 * stats::pnorm(-abs(rg_mat / rg_se))
    diag(rg_p_tmp) <- NA_real_
    rg_p <- rg_p_tmp
    TRUE
  }, error = function(e) {
    message("  [WARNING] rg standard errors/p-values could not be computed ",
            "(", conditionMessage(e), ") -- rg itself, and everything the ",
            "conditional F-statistic depends on (S, I, Ns_gencov), are ",
            "unaffected; rg_se/rg_p are returned as NA.")
    FALSE
  })

  names(N_vec) <- trait_names

  Ns_gencov <- I_mat * outer(N_vec, N_vec, function(a, b) sqrt(a * b))
  dimnames(Ns_gencov) <- list(trait_names, trait_names)

  full_matrix <- I_mat / sqrt(outer(diag(I_mat), diag(I_mat)))
  diag(full_matrix) <- 1   # numerical cleanup only -- already ~1 by construction
  dimnames(full_matrix) <- list(trait_names, trait_names)

  message("  [LDSC] Done. Full ", length(trait_names), "x", length(trait_names),
          " intercept matrix computed.")

  list(
    full_matrix = full_matrix,
    I           = I_mat,
    I_se        = I_se_mat,
    Ns_gencov   = Ns_gencov,
    S           = S_mat,
    rg          = rg_mat,
    rg_se       = rg_se,
    rg_p        = rg_p,
    N           = N_vec,
    trait_names = trait_names
  )
}
