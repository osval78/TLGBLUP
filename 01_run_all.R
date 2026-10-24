###############################################################################
# 01_run_all.R -- fits GBLUP, TLGBLUP (BRR weights), TLGBLUP-F (fixed weights),
# TLGBLUP-H (TLGBLUP + Hadamard interaction S c x g), MEGBLUP and MEGBLUP-M for every
# data set x trait x scheme x partition.  Resumable: results are appended after
# every fit and the (data set, trait, scheme, partition, method) keys already in
# metrics.csv are skipped.  Usage:  Rscript 01_run_all.R [data set names]
# Outputs (../outputs/results): metrics.csv, predictions.csv.gz, transfer_weights.csv,
#   shares.csv, varcomp.csv, sources.csv, timing.csv, partitions.rds, source_w_*.rds
###############################################################################
source("00_functions.R")
R <- file.path("..", "outputs", "results"); tmp <- file.path(tempdir(), "bglr_"); dir.create(R, showWarnings = FALSE)
rd <- function(f) { p <- file.path(R, f); if (file.exists(p)) read.csv(p) else NULL }
app <- function(df, f) { p <- file.path(R, f); write.table(df, p, sep = ",", row.names = FALSE,
                                                           col.names = !file.exists(p), append = file.exists(p)) }
args <- commandArgs(TRUE); RUN <- if (length(args)) args else names(DS)
## optional selection of schemes and methods: SCHEMES="S20,S10" METHODS="angleTL" Rscript 01_run_all.R
SCHEMES <- if (nzchar(Sys.getenv("SCHEMES"))) strsplit(Sys.getenv("SCHEMES"), ",")[[1]] else names(CFG$schemes)
METHODS <- if (nzchar(Sys.getenv("METHODS"))) strsplit(Sys.getenv("METHODS"), ",")[[1]] else CFG$methods
## optional time budget (seconds): no new fit is started if it would not finish in time,
## so the script can be run in chunks (e.g. BUDGET=280 Rscript 01_run_all.R)
BUDGET <- as.numeric(Sys.getenv("BUDGET", "Inf")); T0 <- tnow(); EST <- c(MEGBLUP = 150, `MEGBLUP-M` = 40, other = 15)
over <- function(m) tnow() - T0 + (if (m %in% names(EST)) EST[[m]] else EST[["other"]]) > BUDGET
done <- rd("metrics.csv"); done <- if (is.null(done)) character(0) else paste(done$dataset, done$trait, done$scheme, done$rep, done$method)
partf <- file.path(R, "partitions.rds"); PARTS <- if (file.exists(partf)) readRDS(partf) else list()

for (dn in RUN) {
  D <- load_dataset(dn); srcyears <- setdiff(D$years, D$target)
  tl <- D$Y$Line[D$Y$Env == D$target & complete.cases(D$Y)]
  ## CV1 schemes only for data sets in which the target lines were evaluated in the source environments
  seen <- mean(unique(tl) %in% D$Y$Line[D$Y$Env != D$target]) > 0
  SCH <- SCHEMES[!is_cv1(SCHEMES) | seen]
  for (sc in SCH) {                                  # the same partitions for all traits and methods
    key <- paste(dn, sc)
    if (is_cv1(sc)) { PARTS[[key]] <- PARTS[[paste(dn, base_scheme(sc))]]; next }   # CV1 uses the CV2 partitions
    if (is.null(PARTS[[key]])) { PARTS[[key]] <- make_partitions(tl, CFG$schemes[[sc]], D$nrep,
                                   seed = 1000 * match(dn, names(DS)) + 10 * match(sc, names(CFG$schemes))); saveRDS(PARTS, partf) }
  }
  if (any(is_cv1(SCH)) && any(sapply(PARTS[paste(dn, SCH)], is.null))) stop("run the base scheme partitions first")
  for (tr in D$traits) {
    ## optional subset of partitions (REPS="1" or "1:3"): lets the runs proceed partition by partition
    REPS <- if (nzchar(Sys.getenv("REPS"))) eval(parse(text = Sys.getenv("REPS"))) else seq_len(D$nrep)
    REPS <- intersect(REPS, seq_len(D$nrep))
    todo <- expand.grid(rep = REPS, scheme = SCH, method = METHODS, stringsAsFactors = FALSE)
    todo <- todo[!paste(dn, tr, todo$scheme, todo$rep, todo$method) %in% done, ]
    if (!nrow(todo)) next
    srcf <- file.path(R, sprintf("source_w_%s_%s.rds", dn, tr))
    if (file.exists(srcf)) src <- readRDS(srcf) else {
      set.seed(1); src <- lapply(srcyears, function(yr) fit_source(D, tr, yr, tmp)); saveRDS(src, srcf)
      app(data.frame(dataset = dn, trait = tr, source_year = srcyears, n = sapply(src, `[[`, "n"),
                     varU = sapply(src, `[[`, "varU"), varE = sapply(src, `[[`, "varE"), h2 = sapply(src, `[[`, "h2")), "sources.csv")
    }
    kf <- file.path("..", "outputs", "cache", sprintf("K_%s_%s.rds", dn, tr))
    KER <- if (!any(c("MEGBLUP", "MEGBLUP-M") %in% todo$method)) NULL else if (file.exists(kf)) readRDS(kf) else { t0 <- tnow(); k <- me_kernels(D, trait_data(D, tr))
      if (D$eig < 1) saveRDS(k, kf, compress = FALSE)
      message(sprintf("%s %s: MEGBLUP kernels in %.0f s (ranks %s)", dn, tr, tnow() - t0, paste(k$rank, collapse = "/"))); k }
    for (sc in SCH) for (r in REPS) {
      ms <- todo$method[todo$scheme == sc & todo$rep == r]; if (!length(ms)) next
      P <- PARTS[[paste(dn, sc)]][[r]]; tgt <- trait_data(D, tr); tgt <- tgt$Line[tgt$Env == D$target]
      cv1 <- is_cv1(sc)
      ## CV1: source models refitted without the testing lines of this partition (no information on them leaves
      ## the source environments); their heritabilities are stored in sources_cv1.csv
      srcP <- if (!cv1 || !any(ms %in% c("TLGBLUP", "TLGBLUP-F", "TLGBLUP-H", "angleTL"))) src else {
        set.seed(r); ex <- cv1_test(D, tr, P); s1 <- lapply(srcyears, function(yr) fit_source(D, tr, yr, tmp, exclude = ex))
        app(data.frame(dataset = dn, trait = tr, scheme = sc, rep = r, source_year = srcyears, n = sapply(s1, `[[`, "n"),
                       h2 = sapply(s1, `[[`, "h2")), "sources_cv1.csv"); s1 }
      for (m in ms) {
        if (over(m)) { message("time budget reached"); quit(save = "no", status = 0) }
        set.seed(1000 * r + match(m, CFG$methods)); t0 <- tnow()
        f <- switch(m, GBLUP = fit_gblup(D, tr, P, tmp),
                    TLGBLUP = fit_tlgblup(D, tr, P, srcP, tmp, "BRR"),
                    `TLGBLUP-F` = fit_tlgblup(D, tr, P, srcP, tmp, "FIXED"),
                    `TLGBLUP-H` = fit_tlgblup(D, tr, P, srcP, tmp, "BRR", hadamard = TRUE),
                    MEGBLUP = fit_megblup(D, tr, P, tmp, KER, cv1 = cv1),
                    `MEGBLUP-M` = fit_megblup(D, tr, P, tmp, KER, gxe = FALSE, cv1 = cv1),
                    angleTL = fit_angletl(D, tr, P, srcP, seed = r))
        if (!is.null(f$tuning)) app(data.frame(dataset = dn, trait = tr, scheme = sc, rep = r, t(f$tuning)), "angletl_tuning.csv")
        el <- tnow() - t0; if (m %in% names(EST)) EST[[m]] <- 1.1 * el
        met <- data.frame(dataset = dn, trait = tr, scheme = sc, rep = r, method = m, n_train = length(P$train),
                          n_test = sum(f$test), t(metrics(f$pred[f$test], f$obs[f$test])))
        app(met, "metrics.csv"); app(data.frame(dataset = dn, trait = tr, scheme = sc, rep = r, method = m, seconds = el), "timing.csv")
        if (!is.null(f$varcomp)) app(data.frame(dataset = dn, trait = tr, scheme = sc, rep = r, method = m, component = names(f$varcomp),
                       value = unname(f$varcomp)), "varcomp.csv")
        if (!is.null(f$c)) {
          app(data.frame(dataset = dn, trait = tr, scheme = sc, rep = r, method = m, source_year = srcyears,
                         c = f$c, share = f$share), "transfer_weights.csv")
          app(data.frame(dataset = dn, trait = tr, scheme = sc, rep = r, method = m, share = f$share,
                         share_H = if (is.null(f$share_H)) NA else f$share_H), "shares.csv") }
        pr <- data.frame(dataset = dn, trait = tr, scheme = sc, rep = r, method = m, Line = tgt[f$test],
                         obs = f$obs[f$test], pred = f$pred[f$test])
        p <- file.path(R, "predictions.csv.gz"); con <- gzfile(p, ifelse(file.exists(p), "a", "w"))
        write.table(pr, con, sep = ",", row.names = FALSE, col.names = !file.exists(p) || file.size(p) == 0); close(con)
        message(sprintf("%s %s %s rep %d %-9s %5.0f s | r = %.3f  NRMSE = %.3f", dn, tr, sc, r, m, el, met$r, met$NRMSE))
      }
    }
    rm(KER); gc(verbose = FALSE)
  }
}
message("ALL DONE")
