###############################################################################
# 05_sensitivity.R -- sensitivity analyses added in version 4.  Usage:
#   Rscript 05_sensitivity.R zs ss mcmc exact rotation      (any subset, in any order)
#  zs       TLGBLUP with the source scores standardized to unit variance (TLGBLUP-Zs), schemes S20 and S50
#  ss       "no pooling of genotypes": markers standardized separately within every source environment and
#           within the target environment; marker effects shared on the allele scale (TLGBLUP-SS and the
#           target-only GBLUP-SS with the target-only kernel), schemes S20 and S50
#  mcmc     chain-length check: standard chain vs two chains 5-10 times longer (partition 1 of S50) for
#           GBLUP, TLGBLUP, TLGBLUP-H, MEGBLUP-M and MEGBLUP in one trait of each data-set family;
#           effective sample sizes, Gelman-Rubin statistics and agreement of the predictions
#  exact    MEGBLUP with the exact (full-rank) kernels vs the 95%-trace kernels in LC-Optimal and G2F-2016
#  rotation G2F: every location used in turn as the target (all other locations as sources), 20% training,
#           3 partitions, models GBLUP, TLGBLUP, TLGBLUP-H, angleTL, MEGBLUP-M and MEGBLUP
# Outputs: ../outputs/results/sens_*.csv and rotation_*.csv (resumable: finished keys are skipped)
###############################################################################
source("00_functions.R")
suppressPackageStartupMessages(library(coda))
R <- file.path("..", "outputs", "results"); tmp <- file.path(tempdir(), "bglr_")
app <- function(df, f) { p <- file.path(R, f); write.table(df, p, sep = ",", row.names = FALSE,
                                                           col.names = !file.exists(p), append = file.exists(p)) }
donekeys <- function(f, cols) { p <- file.path(R, f); if (!file.exists(p)) return(character(0))
  x <- read.csv(p); do.call(paste, x[cols]) }
PARTS <- readRDS(file.path(R, "partitions.rds"))
args <- commandArgs(TRUE)
mets <- function(f) metrics(f$pred[f$test], f$obs[f$test])

## ---- zs and ss -------------------------------------------------------------------
if (any(c("zs", "ss") %in% args)) {
  done <- donekeys("sens_scores.csv", c("dataset", "trait", "scheme", "rep", "method"))
  SD <- if (nzchar(Sys.getenv("SENS_DS"))) strsplit(Sys.getenv("SENS_DS"), ",")[[1]] else names(DS)
  for (dn in SD) { D <- load_dataset(dn); srcyears <- setdiff(D$years, D$target)
    for (tr in D$traits) {
      src <- readRDS(file.path(R, sprintf("source_w_%s_%s.rds", dn, tr)))
      if ("ss" %in% args) {
        ## source-specific standardization: Z_k = scale(Z[lines of k]) equals scale(M[lines of k]) (affine
        ## invariance), so the raw marker codes are not needed.  w_k is mapped to the allele scale (w_k / sd_k)
        ## and applied to the target markers centred with the target means.
        dt <- trait_data(D, tr); tl <- dt$Line[dt$Env == D$target]
        std <- function(Z) { m <- colMeans(Z); s <- apply(Z, 2, sd); s[s < 1e-8] <- Inf
                             list(Z = sweep(sweep(Z, 2, m), 2, s, "/"), s = s) }
        set.seed(1)
        wraw <- lapply(srcyears, function(yr) { lk <- dt$Line[dt$Env == yr]; st <- std(D$Z[lk, ])
          Gk <- tcrossprod(st$Z) / ncol(st$Z)
          fm <- BGLR(y = dt[[tr]][dt$Env == yr], ETA = list(list(K = Gk, model = "RKHS")), nIter = D$nIter,
                     burnIn = D$burnIn, thin = D$thin, verbose = FALSE, saveAt = tmp)
          w <- drop(crossprod(st$Z, pinv_sym(Gk) %*% fm$ETA[[1]]$u)) / ncol(st$Z); w / st$s })
        ZTc <- sweep(D$Z[tl, ], 2, colMeans(D$Z[tl, ]))
        SS <- sapply(wraw, function(w) drop(ZTc %*% w)); colnames(SS) <- paste0("s_", srcyears); SS <- scale(SS, scale = FALSE)
        stT <- std(D$Z[tl, ]); GTs <- tcrossprod(stT$Z) / ncol(stT$Z); dimnames(GTs) <- list(tl, tl)
      }
      for (sc in c("S20", "S50")) for (r in seq_len(D$nrep)) { if (is.null(PARTS[[paste(dn, sc)]])) next
        P <- PARTS[[paste(dn, sc)]][[r]]
        todo <- c(if ("zs" %in% args) "TLGBLUP-Zs", if ("ss" %in% args) c("TLGBLUP-SS", "GBLUP-SS"))
        todo <- todo[!paste(dn, tr, sc, r, todo) %in% done]
        for (m in todo) { set.seed(1000 * r + 2); t0 <- tnow()
          f <- switch(m, `TLGBLUP-Zs` = fit_tlgblup(D, tr, P, src, tmp, "BRR", std_scores = TRUE),
                      `TLGBLUP-SS` = fit_tlgblup(D, tr, P, NULL, tmp, "BRR", S = SS, GT = GTs),
                      `GBLUP-SS` = { tt <- target_y(D, tr, P)
                        fm <- BGLR(y = tt$y, ETA = list(list(K = GTs, model = "RKHS")), nIter = D$nIter, burnIn = D$burnIn,
                                   thin = D$thin, verbose = FALSE, saveAt = tmp)
                        list(pred = fm$yHat, obs = tt$dt[[tr]], test = tt$test) })
          app(data.frame(dataset = dn, trait = tr, scheme = sc, rep = r, method = m, t(mets(f)), seconds = tnow() - t0), "sens_scores.csv")
          message(sprintf("%s %s %s %d %s r = %.3f", dn, tr, sc, r, m, mets(f)[1])) } }
    } }
}

## ---- mcmc: chain-length and convergence check ----------------------------------------
if ("mcmc" %in% args) {
  CASES <- list(c("Indica", "GY"), c("G2F_2014", "GM"))
  LONG <- list(Rice = c(nIter = 30000, burnIn = 5000, thin = 5), Maize = c(nIter = 12000, burnIn = 2000, thin = 2))
  done <- donekeys("sens_mcmc.csv", c("dataset", "trait", "method", "chain"))
  for (cs in CASES) { dn <- cs[1]; tr <- cs[2]; D0 <- load_dataset(dn); P <- PARTS[[paste(dn, "S50")]][[1]]
    src <- readRDS(file.path(R, sprintf("source_w_%s_%s.rds", dn, tr)))
    KER <- NULL; kf <- file.path("..", "outputs", "cache", sprintf("K_%s_%s.rds", dn, tr))
    KER <- if (file.exists(kf)) readRDS(kf) else me_kernels(D0, trait_data(D0, tr))
    for (m in c("GBLUP", "TLGBLUP", "TLGBLUP-H", "MEGBLUP-M", "MEGBLUP")) {
      preds <- list()
      for (ch in c("standard", "long1", "long2")) {
        D <- D0; if (ch != "standard") { L <- LONG[[DS[[dn]]$crop]]; D$nIter <- L[["nIter"]]; D$burnIn <- L[["burnIn"]]; D$thin <- L[["thin"]] }
        tg <- file.path(tempdir(), paste0("mc_", ch, "_")); unlink(paste0(tg, "*"))
        set.seed(switch(ch, standard = 1000 + match(m, CFG$methods), long1 = 11, long2 = 22)); t0 <- tnow()
        f <- switch(m, GBLUP = fit_gblup(D, tr, P, tg), TLGBLUP = fit_tlgblup(D, tr, P, src, tg, "BRR"),
                    `TLGBLUP-H` = fit_tlgblup(D, tr, P, src, tg, "BRR", hadamard = TRUE),
                    `MEGBLUP-M` = fit_megblup(D, tr, P, tg, KER, gxe = FALSE), MEGBLUP = fit_megblup(D, tr, P, tg, KER))
        el <- tnow() - t0; preds[[ch]] <- f$pred[f$test]
        ## variance-component chains saved by BGLR (one value per thinned iteration, burn-in included)
        fs <- c(varE = paste0(tg, "varE.dat"), Sys.glob(paste0(tg, "ETA_*_var[UB].dat")))
        names(fs)[-1] <- sub(".*ETA_(.*)_(var[UB]).dat", "\\2_\\1", fs[-1])
        keep <- function(x) x[-seq_len(floor(D$burnIn / D$thin))]
        ch_s <- lapply(fs, function(p) keep(scan(p, quiet = TRUE)))
        saveRDS(ch_s, file.path(tempdir(), sprintf("chains_%s_%s_%s.rds", dn, m, ch)))
        ess <- sapply(ch_s, function(x) unname(effectiveSize(mcmc(x))))
        if (!paste(dn, tr, m, ch) %in% done)
          app(data.frame(dataset = dn, trait = tr, method = m, chain = ch, nIter = D$nIter, burnIn = D$burnIn, thin = D$thin,
                         seconds = el, r = mets(f)[["r"]], parameter = names(ess), ESS = unname(ess),
                         post_mean = sapply(ch_s, mean)), "sens_mcmc.csv")
        message(sprintf("mcmc %s %s %s %s: %.0f s, r = %.3f, min ESS = %.0f", dn, tr, m, ch, el, mets(f)[["r"]], min(ess)))
      }
      ## Gelman-Rubin between the two long chains; agreement of the predictions of the standard and long chains
      c1 <- readRDS(file.path(tempdir(), sprintf("chains_%s_%s_long1.rds", dn, m)))
      c2 <- readRDS(file.path(tempdir(), sprintf("chains_%s_%s_long2.rds", dn, m)))
      psrf <- sapply(names(c1), function(k) gelman.diag(mcmc.list(mcmc(c1[[k]]), mcmc(c2[[k]])), autoburnin = FALSE)$psrf[1, 1]); names(psrf) <- names(c1)
      app(data.frame(dataset = dn, trait = tr, method = m, parameter = names(psrf), PSRF = unname(psrf),
                     cor_pred_std_long = cor(preds$standard, preds$long1), cor_pred_long_long = cor(preds$long1, preds$long2),
                     max_abs_diff_pred_sd = max(abs(preds$standard - preds$long1)) / sd(preds$long1)), "sens_mcmc_agreement.csv")
    }
    rm(KER); invisible(gc()) }
}

## ---- exact: MEGBLUP with exact kernels in two maize data sets ---------------------------
if ("exact" %in% args) {
  done <- donekeys("sens_exact.csv", c("dataset", "trait", "rep"))
  for (cs in list(c("LC_Optimal", "GY"), c("G2F_2016", "GM"))) { dn <- cs[1]; tr <- cs[2]; D <- load_dataset(dn); D$eig <- 1
    dd <- trait_data(D, tr); t0 <- tnow()
    ## exact eigen-decompositions (no truncation), computed once and reused in the five partitions
    e <- eigen(D$G, symmetric = TRUE); k <- e$values > 1e-10 * max(e$values)
    Fm <- e$vectors[match(dd$Line, rownames(D$G)), k] %*% diag(sqrt(e$values[k])); rm(e)
    sv <- svd(Fm, nv = 0); rm(Fm); kk <- sv$d^2 > 1e-10 * max(sv$d^2)
    L <- list(V = sv$u[, kk], d = sv$d[kk]^2, model = "RKHS"); rm(sv)
    blk <- lapply(split(seq_len(nrow(dd)), dd$Env), function(i) { b <- eigen(D$G[dd$Line[i], dd$Line[i]], symmetric = TRUE)
      kk <- b$values > 1e-10 * max(b$values); list(i = i, V = b$vectors[, kk], d = b$values[kk]) })
    r <- sum(sapply(blk, function(b) length(b$d))); VG <- matrix(0, nrow(dd), r); dG <- numeric(r); j <- 0
    for (b in blk) { jj <- j + seq_along(b$d); VG[b$i, jj] <- b$V; dG[jj] <- b$d; j <- max(jj) }
    KER <- list(L = L, GE = list(V = VG, d = dG, model = "RKHS")); rm(VG, blk); invisible(gc())
    message(sprintf("exact kernels %s: %.0f s (ranks %d / %d)", dn, tnow() - t0, length(L$d), r))
    for (rr in 1:3) { if (paste(dn, tr, rr) %in% done) next
      P <- PARTS[[paste(dn, "S50")]][[rr]]; set.seed(1000 * rr + 5); t0 <- tnow()
      f <- fit_megblup(D, tr, P, tmp, KER)
      app(data.frame(dataset = dn, trait = tr, scheme = "S50", rep = rr, method = "MEGBLUP-exact", t(mets(f)),
                     rank_L = length(KER$L$d), rank_GE = length(KER$GE$d), seconds = tnow() - t0), "sens_exact.csv")
      message(sprintf("exact %s rep %d: r = %.3f (%.0f s)", dn, rr, mets(f)[1], tnow() - t0)) }
    rm(KER, L); invisible(gc()) }
}

## ---- rotation: every G2F location as the target ------------------------------------------
if ("rotation" %in% args) {
  RM <- c("GBLUP", "TLGBLUP", "TLGBLUP-H", "angleTL", "MEGBLUP-M")
  done <- donekeys("rotation_metrics.csv", c("dataset", "target", "rep", "method"))
  for (dn in c("G2F_2014", "G2F_2015", "G2F_2016")) { D0 <- load_dataset(dn); tr <- D0$traits[1]
    sf <- file.path(R, sprintf("rotation_sources_%s.rds", dn))
    if (file.exists(sf)) SRC <- readRDS(sf) else { set.seed(1)
      SRC <- setNames(lapply(D0$years, function(yr) fit_source(D0, tr, yr, tmp)), D0$years); saveRDS(SRC, sf) }
    kf <- file.path("..", "outputs", "cache", sprintf("K_%s_%s.rds", dn, tr)); KER <- readRDS(kf)
    for (ti in seq_along(D0$years)) { tg <- D0$years[ti]; D <- D0; D$target <- tg
      dt <- trait_data(D, tr); lines <- dt$Line[dt$Env == tg]
      parts <- make_partitions(lines, 0.20, 3, seed = 50000 + 100 * match(dn, names(DS)) + ti)
      src <- SRC[setdiff(D0$years, tg)]
      for (rr in 1:3) for (m in RM) { if (paste(dn, tg, rr, m) %in% done) next
        set.seed(1000 * rr + match(m, CFG$methods)); t0 <- tnow(); P <- parts[[rr]]
        f <- switch(m, GBLUP = fit_gblup(D, tr, P, tmp), TLGBLUP = fit_tlgblup(D, tr, P, src, tmp, "BRR"),
                    `TLGBLUP-H` = fit_tlgblup(D, tr, P, src, tmp, "BRR", hadamard = TRUE),
                    angleTL = fit_angletl(D, tr, P, src, seed = rr),
                    `MEGBLUP-M` = fit_megblup(D, tr, P, tmp, KER, gxe = FALSE), MEGBLUP = fit_megblup(D, tr, P, tmp, KER))
        app(data.frame(dataset = dn, trait = tr, target = tg, rep = rr, method = m, n_train = length(P$train),
                       n_test = sum(f$test), t(mets(f)), seconds = tnow() - t0, share = if (is.null(f$share)) NA else f$share),
            "rotation_metrics.csv")
        message(sprintf("rotation %s %s rep %d %-9s r = %.3f", dn, tg, rr, m, mets(f)[1])) } }
    rm(KER); invisible(gc()) }
}
message("SENSITIVITY DONE: ", paste(args, collapse = " "))
