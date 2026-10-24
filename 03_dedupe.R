###############################################################################
# 03_dedupe.R -- housekeeping for chunked/resumed runs: removes duplicated keys
# (a fit repeated after an interrupted chunk) keeping the first occurrence, and
# optionally drops whole data sets (Rscript 03_dedupe.R drop <data sets>)
###############################################################################
suppressPackageStartupMessages(library(data.table)); R <- file.path("..", "outputs", "results")
a <- commandArgs(TRUE); drop <- if (length(a) && a[1] == "drop") a[-1] else character(0)
keys <- list(metrics.csv = c("dataset","trait","scheme","rep","method"), timing.csv = c("dataset","trait","scheme","rep","method"),
             shares.csv = c("dataset","trait","scheme","rep","method"), varcomp.csv = c("dataset","trait","scheme","rep","method","component"),
             transfer_weights.csv = c("dataset","trait","scheme","rep","method","source_year"), sources.csv = c("dataset","trait","source_year"),
             predictions.csv.gz = c("dataset","trait","scheme","rep","method","Line"))
for (f in names(keys)) { p <- file.path(R, f); if (!file.exists(p)) next
  x <- fread(p); n0 <- nrow(x); x <- x[!dataset %in% drop]; x <- unique(x, by = keys[[f]])
  fwrite(x, p); message(sprintf("%s: %d -> %d rows", f, n0, nrow(x))) }
if (length(drop)) { pf <- file.path(R, "partitions.rds"); P <- readRDS(pf); P <- P[!sub(" .*", "", names(P)) %in% drop]; saveRDS(P, pf)
  unlink(file.path(R, sprintf("source_w_%s_*.rds", drop))); unlink(file.path("..", "outputs", "cache", sprintf("*%s*", drop))) }
