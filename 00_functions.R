###############################################################################
# TLGBLUP: multi-source transfer-learning GBLUP -- core functions (BGLR)
# Seven real data sets: rice indica and japonica (four traits each), maize LC
# (grain yield under optimal and drought management) and maize G2F 2014, 2015
# and 2016 (grain moisture).  Models: GBLUP, TLGBLUP (BRR weights), TLGBLUP-F
# (fixed weights), TLGBLUP-H (TLGBLUP + Hadamard interaction between the
# transferred term Sc and the target genomic effect g), MEGBLUP (G + GxE) and
# MEGBLUP-M (MEGBLUP without the GxE term: main effects only).
###############################################################################
suppressPackageStartupMessages(library(BGLR))
suppressPackageStartupMessages(library(data.table))

## ---- 0. configuration ---------------------------------------------------------
## nIter/burnIn/thin: MCMC settings of every model fitted to the data set; nrep: number of
## random partitions per scheme.
## coding: "012" = 0/1/2 allele counts; "01" = 0/1 codes; "std" = markers supplied
## already standardized (G2F).  eig: share of the trace of the MEGBLUP kernels kept
## in their eigen-representation (1 = exact).  nIter/burnIn apply to every model
## fitted to that data set (so all comparisons within a data set use the same chain).
DS <- list(
  G2F_2014   = list(X = "X_G2F_2014.csv", Y = "Y_G2F_2014.csv", line = "Hybrid",      coding = "std",
                    traits = c(GM = "Grain_Moisture_BLUE"),                  crop = "Maize", eig = .95, nIter = 12000, burnIn = 2000, thin = 2, nrep = 10),
  G2F_2015   = list(X = "X_G2F_2015.csv", Y = "Y_G2F_2015.csv", line = "Hybrid",      coding = "std",
                    traits = c(GM = "Grain_Moisture_BLUE"),                  crop = "Maize", eig = .95, nIter = 12000, burnIn = 2000, thin = 2, nrep = 10),
  G2F_2016   = list(X = "X_G2F_2016.csv", Y = "Y_G2F_2016.csv", line = "Hybrid",      coding = "std",
                    traits = c(GM = "Grain_Moisture_BLUE"),                  crop = "Maize", eig = .95, nIter = 12000, burnIn = 2000, thin = 2, nrep = 10),
  EYT_4      = list(X = "X_EYT_4.csv",    Y = "Y_EYT_4.csv",    line = "Line",        coding = "012",
                    traits = c(GY = "GY", DTHD = "DTHD", DTMT = "DTMT", PH = "Height"), crop = "Wheat", eig = .85, nIter = 12000, burnIn = 2000, thin = 2, nrep = 5),
  EYT_5      = list(X = "X_EYT_5.csv",    Y = "Y_EYT_5.csv",    line = "Line",        coding = "012",
                    traits = c(GY = "GY", DTHD = "DTHD", DTMT = "DTMT", PH = "Height"), crop = "Wheat", eig = .85, nIter = 12000, burnIn = 2000, thin = 2, nrep = 5),
  EYT_6      = list(X = "X_EYT_6.csv",    Y = "Y_EYT_6.csv",    line = "Line",        coding = "012",
                    traits = c(GY = "GY", DTHD = "DTHD", DTMT = "DTMT", PH = "Height"), crop = "Wheat", eig = .85, nIter = 12000, burnIn = 2000, thin = 2, nrep = 5),
  EYT_7      = list(X = "X_EYT_7.csv",    Y = "Y_EYT_7.csv",    line = "Line",        coding = "012",
                    traits = c(GY = "GY", DTHD = "DTHD", DTMT = "DTMT", PH = "Height"), crop = "Wheat", eig = .85, nIter = 1200, burnIn = 200, thin = 2, nrep = 5))
## Schemes: S10, S20, S50, S80 = proportion of the target lines used for training (CV2-type: the testing
## lines may have records in the source environments).  CV1_S20 and CV1_S50 use the same partitions as S20
## and S50, but the testing lines are removed from the source environments (source models refitted in every
## partition; their records are also masked in MEGBLUP and MEGBLUP-M): prediction of new lines.
## Method order is fixed (it defines the random seed of every fit); angleTL was added in version 4.
CFG <- list(schemes = c(S50 = 0.50, S80 = 0.80, S20 = 0.20, S10 = 0.10, CV1_S20 = 0.20, CV1_S50 = 0.50), maf = 0.05,
            methods = c("GBLUP", "TLGBLUP", "TLGBLUP-F", "TLGBLUP-H", "MEGBLUP", "MEGBLUP-M", "angleTL"))
is_cv1 <- function(sc) grepl("^CV1_", sc)
## testing lines of a partition for a trait: target lines with a record of the trait that are not in training
cv1_test <- function(D, trait, part) { dt <- trait_data(D, trait); setdiff(dt$Line[dt$Env == D$target], part$train) }
base_scheme <- function(sc) sub("^CV1_", "", sc)
RAW <- file.path("..", "data", "raw")

## ---- 1. data: markers, genomic relationship matrix, phenotypes -------------
.XCACHE <- new.env()
load_markers <- function(file, idcol, coding, maf = CFG$maf) {
  if (!is.null(.XCACHE[[file]])) return(.XCACHE[[file]])
  X <- fread(file.path(RAW, file)); ids <- as.character(X[[idcol]])
  mcols <- setdiff(names(X), c(idcol, "V1", "")); M <- as.matrix(X[, ..mcols]); rownames(M) <- ids; rm(X)
  v <- apply(M, 2, var); keep <- v > 0
  if (coding == "012") { p <- colMeans(M) / 2; keep <- keep & pmin(p, 1 - p) >= maf }
  if (coding == "01")  { p <- colMeans(M);     keep <- keep & pmin(p, 1 - p) >= maf }
  if (coding == "std") {                     # standardized codes: frequency of the most common genotype
    fmax <- apply(M, 2, function(m) max(tabulate(match(round(m, 6), unique(round(m, 6))))) / length(m))
    keep <- keep & fmax <= 1 - maf }
  Z <- scale(M[, keep], center = TRUE, scale = TRUE)        # standardized markers (all lines of the data set)
  G <- tcrossprod(Z) / ncol(Z)                              # G = Z Z' / p
  out <- list(Z = Z, G = G, n_raw = ncol(M), n_used = ncol(Z)); .XCACHE[[file]] <- out; out
}
load_dataset <- function(name) {
  cf <- file.path("..", "outputs", "cache", paste0("D_", name, ".rds"))
  if (file.exists(cf)) return(readRDS(cf))
  s <- DS[[name]]; mk <- load_markers(s$X, s$line, s$coding)
  Y <- fread(file.path(RAW, s$Y)); Y <- as.data.frame(Y)
  Y <- data.frame(Env = as.character(Y$Env), Line = as.character(Y[[s$line]]), Y[, unname(s$traits), drop = FALSE])
  names(Y)[-(1:2)] <- names(s$traits)
  Y <- Y[Y$Line %in% rownames(mk$G), ]; Y <- Y[order(Y$Env, Y$Line), ]; rownames(Y) <- NULL
  yrs <- sort(unique(Y$Env))
  D <- list(name = name, Z = mk$Z, G = mk$G, Y = Y, years = yrs, target = tail(yrs, 1), traits = names(s$traits),
       n_markers_raw = mk$n_raw, n_markers = mk$n_used, nIter = s$nIter, burnIn = s$burnIn, thin = s$thin, eig = s$eig, nrep = s$nrep)
  dir.create(dirname(cf), showWarnings = FALSE); saveRDS(D, cf, compress = FALSE); D     # cache (speeds up resumed runs)
}
trait_data <- function(D, trait) D$Y[!is.na(D$Y[[trait]]), c("Env", "Line", trait)]
bglr <- function(D, y, ETA, tag) BGLR(y = y, ETA = ETA, nIter = D$nIter, burnIn = D$burnIn, thin = D$thin,
                                      verbose = FALSE, saveAt = tag)

## ---- 2. partitions of the target environment ---------------------------------
make_partitions <- function(lines, prop, nrep, seed) {
  set.seed(seed); n <- length(lines); ntr <- max(3, round(prop * n))
  lapply(seq_len(nrep), function(r) { tr <- sort(sample.int(n, ntr)); list(train = lines[tr], test = lines[-tr]) })
}

## ---- 3. source models: GBLUP in each source environment, marker effects w_k ---
## y_k = mu_k 1 + g_k + e,  g_k ~ N(0, s2 G_k).  The posterior mean g_hat_k is mapped
## to marker effects  w_k = Z_k' G_k^+ g_hat_k / p  (GBLUP / ridge-regression equivalence;
## G_k^+ = Moore-Penrose inverse, which equals G_k^{-1} when G_k has full rank).
pinv_sym <- function(A, tol = 1e-8) { e <- eigen(A, symmetric = TRUE); k <- e$values > tol * max(e$values)
  e$vectors[, k, drop = FALSE] %*% (t(e$vectors[, k, drop = FALSE]) / e$values[k]) }
fit_source <- function(D, trait, year, tag, exclude = NULL) {
  dk <- trait_data(D, trait); dk <- dk[dk$Env == year & !(dk$Line %in% exclude), ]   # CV1: testing lines removed
  Zk <- D$Z[dk$Line, , drop = FALSE]; Gk <- D$G[dk$Line, dk$Line]
  fm <- bglr(D, dk[[trait]], list(list(K = Gk, model = "RKHS")), tag)
  g <- fm$ETA[[1]]$u; alpha <- pinv_sym(Gk) %*% g
  w <- drop(crossprod(Zk, alpha)) / ncol(Zk)
  list(w = w, year = year, n = nrow(dk), varU = fm$ETA[[1]]$varU, varE = fm$varE,
       h2 = fm$ETA[[1]]$varU / (fm$ETA[[1]]$varU + fm$varE))
}

## ---- 4. models for the target environment ------------------------------------
target_y <- function(D, trait, part) {
  dt <- trait_data(D, trait); dt <- dt[dt$Env == D$target, ]; y <- dt[[trait]]
  y[!(dt$Line %in% part$train)] <- NA; list(dt = dt, y = y, test = !(dt$Line %in% part$train))
}
source_scores <- function(D, lines, sources) {
  S <- sapply(sources, function(s) drop(D$Z[lines, ] %*% s$w)); S <- matrix(S, nrow = length(lines))
  colnames(S) <- paste0("s_", sapply(sources, `[[`, "year")); scale(S, center = TRUE, scale = FALSE)   # centred scores
}
## (a) GBLUP: target data only.  y = mu 1 + g + e
fit_gblup <- function(D, trait, part, tag) {
  t <- target_y(D, trait, part)
  fm <- bglr(D, t$y, list(list(K = D$G[t$dt$Line, t$dt$Line], model = "RKHS")), tag)
  list(pred = fm$yHat, obs = t$dt[[trait]], test = t$test, varcomp = c(varG = fm$ETA[[1]]$varU, varE = fm$varE))
}
## (b) TLGBLUP: y = mu 1 + S c + g + e,  S = [Z w_1, ..., Z w_K] (source genomic scores),
##     c ~ N(0, s2_c I) ("BRR") or flat ("FIXED", TLGBLUP-F)
## (c) TLGBLUP-H: y = mu 1 + S c + g + h + e,  h ~ N(0, s2_h K_SG),  K_SG = K_S o G_T,
##     K_S = S~ S~' / K (linear kernel of the standardized source scores; the covariance of S c
##     up to scale) and o the Hadamard product: h models the interaction between the
##     transferred component S c and the target genomic effect g.
## Optional arguments used in the sensitivity analyses (05_sensitivity.R): std_scores = TRUE standardizes the
## source scores to unit variance before fitting; S and GT supply externally computed scores and kernel.
fit_tlgblup <- function(D, trait, part, sources, tag, weights = c("BRR", "FIXED"), hadamard = FALSE,
                        std_scores = FALSE, S = NULL, GT = NULL) {
  weights <- match.arg(weights); t <- target_y(D, trait, part)
  if (is.null(S)) S <- source_scores(D, t$dt$Line, sources)
  if (std_scores) { S <- scale(S); S[is.na(S)] <- 0 }
  if (is.null(GT)) GT <- D$G[t$dt$Line, t$dt$Line]
  ETA <- list(list(X = S, model = weights), list(K = GT, model = "RKHS"))
  if (hadamard) {
    Ss <- scale(S); Ss[is.na(Ss)] <- 0; KS <- tcrossprod(Ss) / ncol(Ss)
    ETA[[3]] <- list(K = KS * GT, model = "RKHS") }
  fm <- bglr(D, t$y, ETA, tag)
  cvec <- fm$ETA[[1]]$b; tl <- drop(S %*% cvec); vg <- fm$ETA[[2]]$varU
  vh <- if (hadamard) fm$ETA[[3]]$varU * mean(diag(ETA[[3]]$K)) else 0
  list(pred = fm$yHat, obs = t$dt[[trait]], test = t$test, c = cvec,
       varcomp = c(varG = vg, varH = if (hadamard) fm$ETA[[3]]$varU else NA, varE = fm$varE),
       share = var(tl) / (var(tl) + vg + vh), share_H = vh / (var(tl) + vg + vh))
}
## (d) MEGBLUP: all environments pooled.
##     y = mu + E_j (fixed) + g_i (Z_L G Z_L') + ge_ij ((Z_L G Z_L') o (Z_E Z_E')) + e
## Large data sets: the two kernels are passed to BGLR through their eigenvectors and
## eigenvalues (V, d) computed from the factor F = Z_L U D^{1/2} (G = U D U'), and only
## the leading components accounting for a share D$eig of the trace are kept.
trunc_vd <- function(sv, prop) { d <- sv$d^2; k <- d > 1e-10 * max(d); d <- d[k]
  if (prop < 1) { cs <- cumsum(d) / sum(d); k2 <- seq_len(which(cs >= prop)[1]) } else k2 <- seq_along(d)
  list(V = sv$u[, which(k)[k2], drop = FALSE], d = d[k2]) }
me_kernels <- function(D, dd) {
  if (D$eig >= 1) {                                          # exact kernels (rice data sets)
    ZL <- model.matrix(~ 0 + factor(Line, levels = rownames(D$G)), dd)
    KL <- ZL %*% D$G %*% t(ZL); KE <- tcrossprod(model.matrix(~ 0 + factor(Env), dd))
    return(list(L = list(K = KL, model = "RKHS"), GE = list(K = KL * KE, model = "RKHS"))) }
  ev <- function(K) { e <- eigen(K, symmetric = TRUE); list(u = e$vectors, d = sqrt(pmax(e$values, 0))) }
  if (!anyDuplicated(dd$Line)) L <- trunc_vd(ev(D$G[dd$Line, dd$Line]), D$eig) else {
    e <- eigen(D$G, symmetric = TRUE); k <- e$values > 1e-10 * max(e$values)
    Fm <- e$vectors[match(dd$Line, rownames(D$G)), k] %*% diag(sqrt(e$values[k])); rm(e)
    L <- trunc_vd(svd(Fm, nv = 0), D$eig); rm(Fm) }                  # Z_L G Z_L' = F F'
  blk <- lapply(split(seq_len(nrow(dd)), dd$Env), function(i) {        # (Z_L G Z_L') o (Z_E Z_E') is block diagonal
    b <- trunc_vd(ev(D$G[dd$Line[i], dd$Line[i]]), D$eig); list(i = i, V = b$V, d = b$d) })
  r <- sum(sapply(blk, function(b) length(b$d))); VG <- matrix(0, nrow(dd), r); dG <- numeric(r); j <- 0
  for (b in blk) { jj <- j + seq_along(b$d); VG[b$i, jj] <- b$V; dG[jj] <- b$d; j <- max(jj) }
  list(L = list(V = L$V, d = L$d, model = "RKHS"), GE = list(V = VG, d = dG, model = "RKHS"),
       rank = c(L = length(L$d), GE = r))
}
## (e) MEGBLUP-M: particular case of MEGBLUP without the GxE term (main effects only),
##     y = mu + E_j (fixed) + g_i (Z_L G Z_L') + e
## CV1: the records of the testing lines in the source environments are also set to missing.  With BGLR the
## missing records are imputed at every iteration, which is equivalent to removing them while keeping the
## same (pre-computed) kernels.
fit_megblup <- function(D, trait, part, tag, KER, gxe = TRUE, cv1 = FALSE) {
  dd <- trait_data(D, trait); y <- dd[[trait]]
  y[dd$Env == D$target & !(dd$Line %in% part$train)] <- NA
  if (cv1) y[dd$Line %in% cv1_test(D, trait, part)] <- NA
  XE <- model.matrix(~ factor(Env), dd)[, -1, drop = FALSE]
  ETA <- list(list(X = XE, model = "FIXED"), KER$L); if (gxe) ETA[[3]] <- KER$GE
  fm <- bglr(D, y, ETA, tag)
  it <- which(dd$Env == D$target)
  list(pred = fm$yHat[it], obs = dd[[trait]][it], test = !(dd$Line[it] %in% part$train),
       varcomp = c(varG = fm$ETA[[2]]$varU, varGE = if (gxe) fm$ETA[[3]]$varU else NA, varE = fm$varE))
}

## (f) angleTL (Gu et al., 2025) tuned by 5-fold cross-validation, the frequentist benchmark of TLGBLUP.
##     beta = argmin ||y - X beta||^2/n + lambda ||beta||^2 - 2 eta w'beta.  With c = eta/lambda and
##     delta = n lambda / p, its predictions are  y_hat = mu + c s + G_.,tr (G_tr,tr + delta I)^{-1} (y_tr - mu - c s_tr),
##     i.e. kernel ridge regression around the offset c s (s = Z_T w).  The grid covers delta (16 values, 0.01 to 1000 times the mean diagonal of G), c in
##     [-1, 3] (41 values; c = 0 is target-only ridge/GBLUP, c = 1 the distance-based estimator) and the source
##     direction w (the average of the K source estimates or any single source), all chosen by 5-fold CV.
fit_angletl <- function(D, trait, part, sources, nfold = 5, seed = 1) {
  t <- target_y(D, trait, part); S <- source_scores(D, t$dt$Line, sources)
  S <- cbind(mean = rowMeans(S), S); GT <- D$G[t$dt$Line, t$dt$Line]; tr <- which(!t$test); te <- which(t$test)
  cg <- seq(-1, 3, by = .1); dg <- mean(diag(GT)) * 10^seq(-2, 3, length.out = 16)
  pred_fun <- function(itr, ite, s, delta, cc) {           # returns matrix length(ite) x length(cc)
    e <- eigen(GT[itr, itr], symmetric = TRUE); Ki <- e$vectors %*% (t(e$vectors) / (e$values + delta))
    H <- GT[ite, itr, drop = FALSE] %*% Ki; y0 <- t$y[itr]
    sapply(cc, function(c0) { r0 <- y0 - c0 * s[itr]; mu <- mean(r0); mu + c0 * s[ite] + drop(H %*% (r0 - mu)) }) }
  set.seed(seed); fold <- sample(rep_len(seq_len(nfold), length(tr)))
  best <- c(err = Inf, k = 1, d = 1, c = 0)
  for (d in seq_along(dg)) {
    ## eigen-decomposition shared by all c and sources for this fold and delta
    err <- matrix(0, ncol(S), length(cg))
    for (f in seq_len(nfold)) { itr <- tr[fold != f]; ite <- tr[fold == f]
      e <- eigen(GT[itr, itr], symmetric = TRUE); Ki <- e$vectors %*% (t(e$vectors) / (e$values + dg[d]))
      H <- GT[ite, itr, drop = FALSE] %*% Ki
      for (k in seq_len(ncol(S))) for (j in seq_along(cg)) {
        r0 <- t$y[itr] - cg[j] * S[itr, k]; mu <- mean(r0)
        p <- mu + cg[j] * S[ite, k] + drop(H %*% (r0 - mu)); err[k, j] <- err[k, j] + sum((t$y[ite] - p)^2) } }
    m <- which(err == min(err), arr.ind = TRUE)[1, ]
    if (err[m[1], m[2]] < best["err"]) best <- c(err = err[m[1], m[2]], k = unname(m[1]), d = d, c = cg[m[2]]) }
  pr <- rep(NA_real_, nrow(t$dt)); pr[te] <- pred_fun(tr, te, S[, best[["k"]]], dg[best[["d"]]], best[["c"]])[, 1]
  list(pred = pr, obs = t$dt[[trait]], test = t$test,
       tuning = c(source = colnames(S)[best[["k"]]], delta = dg[best[["d"]]], c = best[["c"]]))
}

## ---- 5. metrics ----------------------------------------------------------------
## NRMSE = RMSE / mean (as in version 3); NRMSE_SD = RMSE / SD of the observed values (scale and location free)
metrics <- function(pred, obs) { rmse <- sqrt(mean((obs - pred)^2))
  c(r = cor(pred, obs), NRMSE = rmse / mean(obs), NRMSE_SD = rmse / sd(obs)) }
tnow <- function() unname(proc.time()[3])
