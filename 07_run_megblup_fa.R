###############################################################################
# 07_run_megblup_fa.R -- fits MEGBLUP-FA (factor-analytic multi-environment GBLUP, REML; 06_fa_reml.R)
# with exactly the same data, partitions and masking as MEGBLUP and MEGBLUP-M (01_run_all.R):
#   CV2 schemes (S10, S20, S50, S80): the testing lines of the target environment are set as missing;
#   CV1 schemes (CV1_S20, CV1_S50): their records in the source environments are also set as missing.
# Number of factors: k = 2, or k = 1 when FA(2) is not identified ((J - k)^2 < J + k, i.e. J = 4 environments).
# Resumable.  Usage: Rscript 07_run_megblup_fa.R [data set names]      (REPS / SCHEMES as in 01_run_all.R)
# Outputs (../outputs/results): fa_metrics.csv, fa_predictions.csv.gz, fa_timing.csv, fa_details.csv, fa_rcor.csv
###############################################################################
source("00_functions.R"); source("06_fa_reml.R")
R <- file.path("..", "outputs", "results")
app <- function(df, f) { p <- file.path(R, f); write.table(df, p, sep = ",", row.names = FALSE,
                                                           col.names = !file.exists(p), append = file.exists(p)) }
args <- commandArgs(TRUE); RUN <- if (length(args)) args else c("G2F_2014", "G2F_2015", "G2F_2016", "EYT_4", "EYT_5", "EYT_6", "EYT_7")
SCHEMES <- if (nzchar(Sys.getenv("SCHEMES"))) strsplit(Sys.getenv("SCHEMES"), ",")[[1]] else c("S10", "S20", "S50", "S80", "CV1_S20", "CV1_S50")
PARTS <- readRDS(file.path(R, "partitions.rds"))
done <- if (file.exists(file.path(R, "fa_metrics.csv"))) with(read.csv(file.path(R, "fa_metrics.csv")), paste(dataset, trait, scheme, rep)) else character(0)
for (dn in RUN) {
  D <- load_dataset(dn)
  for (tr in D$traits) {
    dd <- trait_data(D, tr); J <- length(unique(dd$Env)); k <- fa_k(J); init <- NULL
    REPS <- if (nzchar(Sys.getenv("REPS"))) eval(parse(text = Sys.getenv("REPS"))) else seq_len(D$nrep)
    for (sc in SCHEMES) for (r in intersect(REPS, seq_len(D$nrep))) {
      if (paste(dn, tr, sc, r) %in% done) next
      P <- PARTS[[paste(dn, base_scheme(sc))]][[r]]; y <- dd[[tr]]
      te <- dd$Env == D$target & !(dd$Line %in% P$train); test <- dd$Line[te]; y[te] <- NA
      if (is_cv1(sc)) y[dd$Line %in% cv1_test(D, tr, P)] <- NA
      t0 <- tnow()
      f <- tryCatch(fa_reml(dd[, c("Env", "Line")], y, D$G, D$target, test, k = k, init = init), error = function(e) { message(conditionMessage(e)); NULL })
      if (is.null(f) || !f$converged) {           # restart from default starting values
        f2 <- tryCatch(fa_reml(dd[, c("Env", "Line")], y, D$G, D$target, test, k = k), error = function(e) NULL)
        if (!is.null(f2) && (is.null(f) || f2$ll > f$ll)) f <- f2 }
      el <- tnow() - t0
      init <- list(Lam = f$Lambda, psi = unname(f$psi), s2 = unname(f$s2))
      obs <- dd[[tr]][te]; pred <- unname(f$pred[test])
      met <- data.frame(dataset = dn, trait = tr, scheme = sc, rep = r, method = "MEGBLUP-FA", n_train = length(P$train),
                        n_test = length(test), t(metrics(pred, obs)))
      app(met, "fa_metrics.csv"); app(data.frame(dataset = dn, trait = tr, scheme = sc, rep = r, method = "MEGBLUP-FA", seconds = el), "fa_timing.csv")
      envs <- rownames(f$rcor); src <- setdiff(envs, D$target)
      app(data.frame(dataset = dn, trait = tr, scheme = sc, rep = r, k = k, J = J, iter = f$iter, converged = f$converged, logL = f$ll,
                     max_score = f$max_score, rank = f$rank, n_lines = f$n_lines,
                     h2_target = f$Sigma_g[D$target, D$target] / (f$Sigma_g[D$target, D$target] + f$s2[[D$target]]),
                     rg_mean = mean(f$rcor[D$target, src]), rg_min = min(f$rcor[D$target, src]), rg_max = max(f$rcor[D$target, src])), "fa_details.csv")
      app(data.frame(dataset = dn, trait = tr, scheme = sc, rep = r, source = src, rg = f$rcor[D$target, src],
                     var_target = f$Sigma_g[D$target, D$target], var_source = diag(f$Sigma_g)[src]), "fa_rcor.csv")
      pr <- data.frame(dataset = dn, trait = tr, scheme = sc, rep = r, method = "MEGBLUP-FA", Line = test, obs = obs, pred = pred)
      p <- file.path(R, "fa_predictions.csv.gz"); con <- gzfile(p, ifelse(file.exists(p), "a", "w"))
      write.table(pr, con, sep = ",", row.names = FALSE, col.names = !file.exists(p) || file.size(p) == 0); close(con)
      message(sprintf("%s %s %-7s rep %2d  k=%d  %3d it  %6.1f s | r = %.3f  conv = %s", dn, tr, sc, r, k, f$iter, el, met$r, f$converged))
    }
  }
  rm(D); .XCACHE <- new.env(); gc(verbose = FALSE)
}
message("ALL DONE")
