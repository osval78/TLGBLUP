###############################################################################
# 02_figures_tables.R -- figures (PNG 600 dpi + PDF) and tables (CSV + Markdown), version 5
# (version 5 adds the four wheat EYT data sets: Figure 3b with their phenotypes; crop = Rice, Maize or Wheat)
# Main-text figures: 1 workflow, 2-3 phenotypes, 4 accuracy vs training size, 5 gains vs GBLUP (corrected CIs),
# 6 CV1 vs CV2, 7 TLGBLUP vs angleTL, 8 transfer weights, 9 transfer share, 10 G2F target rotation.
# Supplementary figures S1-S4.  Statistical tests: corrected resampled t test (Nadeau & Bengio, 2003).
###############################################################################
suppressPackageStartupMessages({ library(ggplot2); library(data.table); library(grid) })
source("00_functions.R")
## version 6: data sets can be excluded (DROP_DATASETS="Indica,Japonica,LC_Optimal,LC_Drought") and the outputs written
## to tagged folders (OUT_TAG="v6" -> outputs/figures_v6, outputs/tables_v6); every results file is filtered on read.
DROP <- strsplit(Sys.getenv("DROP_DATASETS"), ",")[[1]]; TAG <- Sys.getenv("OUT_TAG")
DS <- DS[setdiff(names(DS), DROP)]
fread0 <- data.table::fread
fread <- function(...) { d <- fread0(...); if (length(DROP) && "dataset" %in% names(d)) d <- d[!dataset %in% DROP]; d }
R <- file.path("..", "outputs", "results")
FD <- file.path("..", "outputs", paste0("figures", if (nzchar(TAG)) paste0("_", TAG) else ""))
TD <- file.path("..", "outputs", paste0("tables", if (nzchar(TAG)) paste0("_", TAG) else ""))
dir.create(FD, showWarnings = FALSE); dir.create(TD, showWarnings = FALSE)
M <- fread(file.path(R, "metrics.csv")); W <- fread(file.path(R, "transfer_weights.csv")); SH <- fread(file.path(R, "shares.csv"))
S <- fread(file.path(R, "sources.csv")); TM <- fread(file.path(R, "timing.csv"))
METH <- c("GBLUP", "MEGBLUP-M", "MEGBLUP", "angleTL", "TLGBLUP-F", "TLGBLUP", "TLGBLUP-H")
COL <- c(GBLUP = "#8C8C8C", `MEGBLUP-M` = "#F0E442", MEGBLUP = "#E69F00", angleTL = "#009E73", `TLGBLUP-F` = "#56B4E9", TLGBLUP = "#0072B2", `TLGBLUP-H` = "#CC79A7")
SCH <- c(S10 = "CV2, 10% training", S20 = "CV2, 20% training", S50 = "CV2, 50% training", S80 = "CV2, 80% training",
         CV1_S20 = "CV1, 20% training", CV1_S50 = "CV1, 50% training")
PROP <- c(S10 = 10, S20 = 20, S50 = 50, S80 = 80, CV1_S20 = 20, CV1_S50 = 50)
DSL <- c(Indica = "Indica", Japonica = "Japonica", LC_Optimal = "LC-Optimal", LC_Drought = "LC-Drought",
         G2F_2014 = "G2F-2014", G2F_2015 = "G2F-2015", G2F_2016 = "G2F-2016",
         EYT_4 = "EYT-4", EYT_5 = "EYT-5", EYT_6 = "EYT-6", EYT_7 = "EYT-7")
CROPCOL <- c(Maize = "#D55E00", Rice = "#0072B2", Wheat = "#009E73")
CASES <- unlist(lapply(names(DS), function(d) paste(DSL[d], names(DS[[d]]$traits))))
lab <- function(x) { x[, case := factor(paste(DSL[dataset], trait), CASES)]; x[, crop := sapply(as.character(dataset), function(d) DS[[d]]$crop)]
  if ("scheme" %in% names(x)) x[, schemeL := factor(SCH[scheme], SCH)]; if ("scheme" %in% names(x)) x[, design := ifelse(grepl("^CV1", scheme), "CV1", "CV2")]; if ("method" %in% names(x)) x[, method := factor(method, METH)]; x }
M <- lab(M); W <- lab(W); SH <- lab(SH)
th <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        strip.background = element_rect(fill = "grey95"), strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 40, hjust = 1))
save_fig <- function(p, f, w, h) { ggsave(file.path(FD, paste0(f, ".png")), p, width = w, height = h, dpi = 600)
                                   ggsave(file.path(FD, paste0(f, ".pdf")), p, width = w, height = h) }
ci <- function(x) { m <- mean(x); s <- sd(x) / sqrt(length(x)); q <- qt(.975, length(x) - 1); c(mean = m, lo = m - q * s, hi = m + q * s) }

## ---- Figure 1: workflow ------------------------------------------------------
png(file.path(FD, "Figure1_TLGBLUP_workflow.png"), width = 10, height = 5.4, units = "in", res = 600)
grid.newpage(); pushViewport(viewport(xscale = c(0, 100), yscale = c(0, 54)))
box <- function(x, y, w, h, fill) grid.roundrect(unit(x, "native"), unit(y, "native"), unit(w, "native"), unit(h, "native"),
  r = unit(2, "mm"), gp = gpar(fill = fill, col = "grey30"))
txt <- function(x, y, lab, cex = .78, ...) grid.text(lab, unit(x, "native"), unit(y, "native"), gp = gpar(cex = cex, ...))
arr <- function(x0, y0, x1, y1) grid.lines(unit(c(x0, x1), "native"), unit(c(y0, y1), "native"),
  arrow = arrow(length = unit(2, "mm"), type = "closed"), gp = gpar(fill = "grey30", col = "grey30"))
ys <- c(47, 36, 25, 14)
for (k in 1:4) {
  box(11, ys[k], 19, 8, "#FDE9C9"); txt(11, ys[k] + 1.6, if (k < 4) paste("Source environment", k) else "Source environment K", fontface = "bold")
  txt(11, ys[k] - 1.8, if (k < 4) bquote(bold(y)[.(k)] == mu[.(k)] * bold("1") + bold(g)[.(k)] + bold(e)[.(k)]) else bquote(bold(y)[K] == mu[K] * bold("1") + bold(g)[K] + bold(e)[K]))
  box(36, ys[k], 19, 8, "#FFF7E6"); txt(36, ys[k] + 1.6, "Marker effects (summary)", fontface = "bold")
  txt(36, ys[k] - 1.8, if (k < 4) bquote(hat(bold(w))[.(k)] == bold(Z)[.(k)]^T * bold(G)[.(k)]^"+" * hat(bold(g))[.(k)] / p) else bquote(hat(bold(w))[K] == bold(Z)[K]^T * bold(G)[K]^"+" * hat(bold(g))[K] / p))
  arr(20.5, ys[k], 26.5, ys[k]); arr(45.5, ys[k], 54.5, 30.5) }
box(63, 30.5, 17, 12, "#DCEBF7"); txt(63, 33.5, "Source genomic scores", fontface = "bold")
txt(63, 29, bquote(bold(s)[k] == bold(Z)[T] * hat(bold(w))[k]), cex = .9); txt(63, 25.8, "(all target lines)", cex = .7)
box(87, 43, 24, 16, "#C6DBEF"); txt(87, 48.5, "Target environment (training lines)", fontface = "bold")
txt(87, 44.5, bquote("TLGBLUP: " * bold(y) == mu * bold("1") + bold(S) * bold(c) + bold(g) + bold(e)), cex = .74)
txt(87, 40.8, bquote("TLGBLUP-H: " * "... " + bold(h) * ",  " * bold(h) %~% N(bold(0), sigma[h]^2 * bold(K)[S] * " ∘ " * bold(G)[T])), cex = .66)
txt(87, 37.3, bquote(c[k] %~% N(0, sigma[c]^2) * ",   " * bold(g) %~% N(bold(0), sigma[g]^2 * bold(G)[T])), cex = .7)
box(87, 15, 24, 10, "#9ECAE1"); txt(87, 17.5, "Prediction of testing lines", fontface = "bold")
txt(87, 13.5, bquote(hat(y)[i] == hat(mu) + sum(hat(c)[k] * s[ik], k == 1, K) + hat(g)[i] + hat(h)[i]))
arr(71.5, 32.5, 75, 37); arr(87, 35, 87, 20)
txt(50, 3, "Only the p-dimensional vectors of marker effects leave the source environments (summary-level, multi-source transfer)", cex = .72, fontface = "italic")
txt(50, 0.8, expression("In TLGBLUP-H, " * bold(K)[S] == tilde(bold(S)) * tilde(bold(S))^T / K * " and ∘ is the Hadamard product (interaction between " * bold(S) * bold(c) * " and " * bold(g) * ")"), cex = .66)
dev.off()

## ---- Table 1: data description ---------------------------------------------
TB1 <- list(); YY <- list()
for (dn in names(DS)) { D <- load_dataset(dn); y <- as.data.table(D$Y)
  for (tr in D$traits) { yt <- y[!is.na(get(tr))]
    TB1[[paste(dn, tr)]] <- data.table(dataset = dn, crop = DS[[dn]]$crop, trait = tr, environments = uniqueN(yt$Env), sources = uniqueN(yt$Env) - 1,
      target = D$target, lines = uniqueN(yt$Line), records = nrow(yt), target_lines = sum(yt$Env == D$target),
      target_seen_before = mean(yt[Env == D$target]$Line %in% yt[Env != D$target]$Line),
      markers_raw = D$n_markers_raw, markers_used = D$n_markers)
    YY[[paste(dn, tr)]] <- data.table(dataset = dn, trait = tr, Env = yt$Env, value = yt[[tr]], role = ifelse(yt$Env == D$target, "Target", "Source")) } }
t1 <- rbindlist(TB1); fwrite(t1, file.path(TD, "Table1_data.csv"))
EN <- rbindlist(lapply(names(DS), function(dn) { D <- load_dataset(dn); D$Y |> as.data.table() -> y
  y[, .(lines = .N), by = Env][, `:=`(dataset = dn, role = ifelse(Env == D$target, "target", "source"))] }))
fwrite(EN[, .(dataset, Env, lines, role)], file.path(TD, "TableS1_environments.csv"))
rm(D); .XCACHE <- new.env(); invisible(gc())
YY <- lab(rbindlist(YY))
if (nrow(YY[crop == "Rice"])) {
p2 <- ggplot(YY[crop == "Rice"], aes(Env, value, fill = role)) + geom_boxplot(outlier.size = .4, linewidth = .3) +
  facet_wrap(DSL[dataset] ~ trait, scales = "free", ncol = 4) + scale_fill_manual(values = c(Source = "#FDE9C9", Target = "#9ECAE1")) +
  labs(x = "Year (environment)", y = "Phenotype", fill = NULL) + th + theme(axis.text.x = element_text(angle = 0, hjust = .5))
save_fig(p2, "Figure2_phenotypes_rice", 10, 6.2) }
YM <- YY[crop == "Maize"]; YM[, EnvS := sub("_20[0-9]{2}$", "", Env)]
p3 <- ggplot(YM, aes(EnvS, value, fill = role)) + geom_boxplot(outlier.size = .3, linewidth = .25) +
  facet_wrap(~ case, scales = "free", ncol = 2) + scale_fill_manual(values = c(Source = "#FDE9C9", Target = "#9ECAE1")) +
  labs(x = if (any(grepl("^LC", YM$dataset))) "Environment (year for LC; location for G2F)" else "Environment (location)", y = if (any(grepl("^LC", YM$dataset))) "Phenotype" else "Grain moisture (%)", fill = NULL) + th + theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6.5))
save_fig(p3, "Figure3_phenotypes_maize", 10, 8.5)
YW <- YY[crop == "Wheat"]
p3b <- ggplot(YW, aes(Env, value, fill = role)) + geom_boxplot(outlier.size = .3, linewidth = .25) +
  facet_wrap(~ case, scales = "free", ncol = 4) + scale_fill_manual(values = c(Source = "#FDE9C9", Target = "#9ECAE1")) +
  labs(x = "Environment (management condition)", y = "Phenotype", fill = NULL) + th + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
save_fig(p3b, "Figure3b_phenotypes_wheat", 10, 9)


## ---- metrics recomputed from the predictions of every testing line ---------------------
## r, NRMSE = RMSE / SD of the observed values (scale- and location-free), NRMSE_mean = RMSE / mean (version 3)
## and the coincidence index: share of the observed best 20% of the testing lines that are also among the
## predicted best 20%.  "Best" = largest values for GY and PHR and smallest values for GC, PH and GM.
HIGH <- c(GY = TRUE, PHR = TRUE, GC = FALSE, PH = FALSE, GM = FALSE, DTHD = FALSE, DTMT = FALSE)   # wheat: earliness is favourable
coinc <- function(obs, pred, high) { k <- max(1, round(.2 * length(obs))); sg <- if (high) -1 else 1
  length(intersect(order(sg * obs)[1:k], order(sg * pred)[1:k])) / k }
PR <- fread(file.path(R, "predictions.csv.gz"))
M <- PR[, .(n_test = .N, r = cor(pred, obs), NRMSE = sqrt(mean((obs - pred)^2)) / sd(obs),
            NRMSE_mean = sqrt(mean((obs - pred)^2)) / mean(obs), CI20 = coinc(obs, pred, HIGH[trait[1]])),
        by = .(dataset, trait, scheme, rep, method)]
NTR <- fread(file.path(R, "metrics.csv"))[, .(n_train = n_train[1]), by = .(dataset, trait, scheme, rep, method)]
M <- merge(M, NTR, by = c("dataset", "trait", "scheme", "rep", "method"))
rm(PR); invisible(gc())
M <- lab(M); fwrite(M[, -c("case", "schemeL", "crop", "design")], file.path(TD, "TableS5_metrics_by_partition.csv"))
MAINM <- setdiff(METH, "TLGBLUP-F")

## ---- Table 2 (per data set x trait) and Table 3 (overall summary) -------------------------
t2 <- M[, .(r_mean = mean(r), r_sd = sd(r), NRMSE_mean = mean(NRMSE), NRMSE_sd = sd(NRMSE), CI20_mean = mean(CI20),
            n_train = n_train[1], n_test = n_test[1]), by = .(dataset, case, crop, scheme, design, trait, method)][order(case, scheme, method)]
fwrite(t2[, -"case"], file.path(TD, "Table2_performance.csv"))
ov <- rbind(t2[, .(r = mean(r_mean), NRMSE = mean(NRMSE_mean), CI20 = mean(CI20_mean), n_cases = .N), by = .(crop, scheme, method)],
            t2[dataset %in% c("Indica", "Japonica", "G2F_2014", "G2F_2015", "G2F_2016"),
               .(crop = "Rice+G2F", r = mean(r_mean), NRMSE = mean(NRMSE_mean), CI20 = mean(CI20_mean), n_cases = .N), by = .(scheme, method)],
            t2[!dataset %in% c("LC_Optimal", "LC_Drought"),          # all data sets with CV1 (rice, G2F and wheat EYT)
               .(crop = "Rice+G2F+Wheat", r = mean(r_mean), NRMSE = mean(NRMSE_mean), CI20 = mean(CI20_mean), n_cases = .N), by = .(scheme, method)],
            t2[crop %in% c("Rice", "Maize") & !grepl("^CV1", scheme),   # the 13 combinations of version 4
               .(crop = "Rice+Maize", r = mean(r_mean), NRMSE = mean(NRMSE_mean), CI20 = mean(CI20_mean), n_cases = .N), by = .(scheme, method)],
            t2[!grepl("^CV1", scheme), .(crop = "All", r = mean(r_mean), NRMSE = mean(NRMSE_mean), CI20 = mean(CI20_mean), n_cases = .N), by = .(scheme, method)])
fwrite(ov[order(scheme, crop, method)], file.path(TD, "Table3_overall_summary.csv"))

## ---- paired differences vs GBLUP: corrected resampled t test --------------------------------
## Var(mean difference) = (1/J + n_test/n_train) s^2 (Nadeau & Bengio, 2003), J = 10 partitions; Holm adjustment over
## the data set x trait combinations within each scheme and model (J = 10 partitions; J = 5 in the wheat EYT data sets).
nb <- function(d, ntr, nte) { J <- length(d); m <- mean(d); v <- (1 / J + nte / ntr) * var(d); se <- sqrt(v); q <- qt(.975, J - 1)
  c(mean = m, lo = m - q * se, hi = m + q * se, p = 2 * pt(-abs(m / se), J - 1), lo_naive = m - q * sd(d) / sqrt(J), hi_naive = m + q * sd(d) / sqrt(J)) }
B <- M[method == "GBLUP", .(dataset, trait, scheme, rep, r0 = r, n0 = NRMSE, c0 = CI20)]
DD <- merge(M[method != "GBLUP"], B, by = c("dataset", "trait", "scheme", "rep"))
DD[, `:=`(dr = r - r0, dn = 100 * (NRMSE - n0) / n0, dc = CI20 - c0)]
G5 <- DD[, { a <- nb(dr, n_train[1], n_test[1]); b <- nb(dn, n_train[1], n_test[1])
  .(dr = a[["mean"]], dr_lo = a[["lo"]], dr_hi = a[["hi"]], p_r = a[["p"]], dr_lo_naive = a[["lo_naive"]], dr_hi_naive = a[["hi_naive"]],
    dn = b[["mean"]], dn_lo = b[["lo"]], dn_hi = b[["hi"]], dc = mean(dc), win_r = mean(dr > 0), win_n = mean(dn < 0)) },
  by = .(dataset, case, crop, trait, scheme, schemeL, design, method)]
G5[, p_r_holm := p.adjust(p_r, "holm"), by = .(scheme, method)]
G5[, sig := fifelse(dr_lo > 0, "gain", fifelse(dr_hi < 0, "loss", "ns"))]
G5[, sig_naive := fifelse(dr_lo_naive > 0, "gain", fifelse(dr_hi_naive < 0, "loss", "ns"))]
fwrite(G5[order(case, scheme, method), -c("schemeL", "case")], file.path(TD, "Table4_differences_vs_GBLUP.csv"))
SIG <- G5[, .(cases = .N, mean_dr = mean(dr), gains_corrected = sum(sig == "gain"), losses_corrected = sum(sig == "loss"),
              gains_holm = sum(p_r_holm < .05 & dr > 0), losses_holm = sum(p_r_holm < .05 & dr < 0),
              gains_naive = sum(sig_naive == "gain"), losses_naive = sum(sig_naive == "loss"), wins_80pct = sum(win_r >= .8)),
          by = .(scheme, method)][order(scheme, method)]
fwrite(SIG, file.path(TD, "Table4b_significance_counts.csv"))

## ---- Figure 4: accuracy as a function of the training proportion (CV2) -------------------------
A4 <- t2[design == "CV2" & method %in% MAINM, .(r = r_mean, sd = r_sd), by = .(case, scheme, method)]; A4[, prop := PROP[scheme]]
A4[, method := factor(method, METH)]
p4 <- ggplot(A4, aes(prop, r, colour = method)) + geom_line(linewidth = .45) + geom_point(size = 1.2) +
  facet_wrap(~ case, ncol = 6, scales = "free_y") + scale_colour_manual(values = COL) + scale_x_continuous(breaks = c(10, 20, 50, 80)) +
  labs(x = "Target lines used for training (%)", y = "Mean Pearson correlation (10 or 5 partitions)", colour = "Model") + th +
  theme(axis.text.x = element_text(angle = 0, hjust = .5))
save_fig(p4, "Figure4_accuracy_vs_training", 12, 10.5)

## ---- Figure 5: gains over GBLUP with corrected confidence intervals (CV2) -----------------------
pd <- position_dodge(.7)
p5 <- ggplot(G5[design == "CV2" & method %in% MAINM], aes(case, dr, colour = method)) + geom_hline(yintercept = 0, linetype = 2) +
  geom_pointrange(aes(ymin = dr_lo, ymax = dr_hi), position = pd, size = .2, linewidth = .35) +
  facet_wrap(~ schemeL, ncol = 1, scales = "free_y") + scale_colour_manual(values = COL[-1]) +
  labs(x = NULL, y = "Difference in Pearson correlation vs GBLUP\n(mean and corrected 95% CI)", colour = "Model") + th
save_fig(p5, "Figure5_gain_vs_GBLUP", 14, 11)

## ---- Figure 6: new lines (CV1) vs lines evaluated in the sources (CV2) ----------------------------
C6 <- G5[scheme %in% c("S20", "S50", "CV1_S20", "CV1_S50") & method %in% MAINM & !dataset %in% c("LC_Optimal", "LC_Drought")]
if (nrow(C6[design == "CV1"])) {
  C6[, train := paste0(PROP[scheme], "% training")]
  p6 <- ggplot(C6, aes(case, dr, colour = method, shape = design)) + geom_hline(yintercept = 0, linetype = 2) +
    geom_point(position = position_dodge(.75), size = 1.8) + facet_wrap(~ train, ncol = 1) +
    scale_colour_manual(values = COL[-1]) + scale_shape_manual(values = c(CV1 = 17, CV2 = 16)) +
    labs(x = NULL, y = "Mean difference in Pearson correlation vs GBLUP", colour = "Model",
         shape = "Testing lines in the sources") + th
  save_fig(p6, "Figure6_CV1_vs_CV2", 13, 7.5) }

## ---- Figure 7: TLGBLUP vs angleTL (partition by partition) -------------------------------------------
TA <- dcast(M[method %in% c("TLGBLUP", "angleTL", "TLGBLUP-H", "MEGBLUP", "MEGBLUP-M", "GBLUP")],
            dataset + case + crop + trait + scheme + schemeL + design + rep ~ method, value.var = c("r", "NRMSE"))
setnames(TA, gsub("-", "_", names(TA)))
if ("r_angleTL" %in% names(TA)) {
  p7 <- ggplot(TA[!is.na(r_angleTL) & !is.na(r_TLGBLUP)], aes(r_angleTL, r_TLGBLUP, colour = crop)) + geom_abline(linetype = 2) +
    geom_point(size = .9, alpha = .6) + facet_wrap(~ schemeL, ncol = 3) + coord_equal() +
    scale_colour_manual(values = CROPCOL) +
    labs(x = "Pearson correlation, angleTL (cross-validated)", y = "Pearson correlation, TLGBLUP", colour = NULL) + th +
    theme(axis.text.x = element_text(angle = 0))
  save_fig(p7, "Figure7_TLGBLUP_vs_angleTL", 10, 7.2) }
fwrite(TA[, .(TL_minus_angleTL = mean(r_TLGBLUP - r_angleTL), win_TL_angleTL = mean(r_TLGBLUP > r_angleTL),
              TL_minus_ME = mean(r_TLGBLUP - r_MEGBLUP), win_TL_ME = mean(r_TLGBLUP > r_MEGBLUP),
              TLH_minus_ME = mean(r_TLGBLUP_H - r_MEGBLUP), win_TLH_ME = mean(r_TLGBLUP_H > r_MEGBLUP),
              TL_minus_MEM = mean(r_TLGBLUP - r_MEGBLUP_M), win_TL_MEM = mean(r_TLGBLUP > r_MEGBLUP_M),
              H_minus_TL = mean(r_TLGBLUP_H - r_TLGBLUP), win_H_TL = mean(r_TLGBLUP_H > r_TLGBLUP)),
           by = .(dataset, case, trait, scheme)][order(case, scheme), -"case"], file.path(TD, "Table5_pairwise_comparisons.csv"))
AT <- fread(file.path(R, "angletl_tuning.csv")); fwrite(AT[, .(c_mean = mean(c), share_mean_source = mean(source == "mean"),
  delta_median = median(delta)), by = .(dataset, trait, scheme)], file.path(TD, "TableS6_angleTL_tuning.csv"))

## ---- Figures 8-9: transfer weights and transfer share ------------------------------------------------
W <- lab(fread(file.path(R, "transfer_weights.csv"))); SH <- lab(fread(file.path(R, "shares.csv")))
W[, srcL := sub("_20[0-9]{2}$", "", source_year)]
p8 <- ggplot(W[method == "TLGBLUP" & scheme %in% c("S20", "S50", "S80")], aes(srcL, c, fill = schemeL)) + geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(outlier.size = .25, linewidth = .2) + facet_wrap(~ case, scales = "free", ncol = 5) +
  scale_fill_manual(values = c("#DEEBF7", "#9ECAE1", "#2171B5")) + labs(x = "Source environment", y = expression("Transfer weight " * hat(c)[k]), fill = NULL) + th +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6))
save_fig(p8, "Figure8_transfer_weights", 12, 13)
fwrite(W[, .(c_mean = mean(c), c_sd = sd(c)), by = .(dataset, trait, scheme, method, source_year)], file.path(TD, "TableS3_transfer_weights.csv"))
TSH <- SH[method %in% c("TLGBLUP", "TLGBLUP-H"), .(share_mean = mean(share), share_sd = sd(share), shareH_mean = mean(share_H)), by = .(dataset, case, trait, scheme, method)]
fwrite(TSH[, -"case"], file.path(TD, "Table6_transfer_share.csv"))
p9 <- ggplot(SH[method == "TLGBLUP"], aes(case, share, fill = schemeL)) + geom_boxplot(outlier.size = .25, linewidth = .2, position = position_dodge(.85)) +
  scale_fill_manual(values = c("#EFF3FF", "#BDD7E7", "#6BAED6", "#2171B5", "#FDD0A2", "#E6550D")) + coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "Share of the target genetic variance\nexplained by the transferred scores (TLGBLUP)", fill = NULL) + th
save_fig(p9, "Figure9_transfer_share", 14, 5.5)

## ---- Figure 10: G2F, every location as the target -----------------------------------------------------
rf <- file.path(R, "rotation_metrics.csv")
if (file.exists(rf)) {
  RO <- fread(rf); RO[, method := factor(method, METH)]
  RA <- RO[, .(r = mean(r), NRMSE = mean(NRMSE_SD), share = mean(share)), by = .(dataset, target, method)]
  RB <- merge(RA[method != "GBLUP"], RA[method == "GBLUP", .(dataset, target, r0 = r, n0 = NRMSE)], by = c("dataset", "target"))
  RB[, dr := r - r0]
  RS <- RB[, .(targets = .N, mean_dr = mean(dr), median_dr = median(dr), sd_dr = sd(dr), p_t = t.test(dr)$p.value,
               p_wilcox = suppressWarnings(wilcox.test(dr)$p.value), targets_better = mean(dr > 0)), by = .(method)][order(method)]
  RS2 <- RB[, .(targets = .N, mean_dr = mean(dr), targets_better = mean(dr > 0)), by = .(dataset, method)]
  RG <- RA[, .(mean_r = mean(r)), by = .(method)]
  fwrite(merge(RS, RG, by = "method"), file.path(TD, "Table7_rotation_summary.csv")); fwrite(RS2, file.path(TD, "TableS7_rotation_by_year.csv"))
  fwrite(RA, file.path(TD, "TableS8_rotation_by_target.csv"))
  RB[, yearL := sub("_", "-", dataset)]
  p10 <- ggplot(RB[method %in% MAINM], aes(method, dr, fill = method)) + geom_hline(yintercept = 0, linetype = 2) +
    geom_boxplot(outlier.size = .5, linewidth = .3) + geom_jitter(width = .12, size = .5, alpha = .5) + facet_wrap(~ yearL, nrow = 1) +
    scale_fill_manual(values = COL, guide = "none") +
    labs(x = NULL, y = "Difference in Pearson correlation vs GBLUP\n(mean of 5 partitions per target location)") + th
  save_fig(p10, "Figure10_G2F_rotation", 11, 4.8) }

## ---- Supplementary figures ------------------------------------------------------------------------
lim <- M[, quantile(NRMSE, .995)]
pS1 <- ggplot(M[design == "CV2"], aes(case, NRMSE, fill = method)) + geom_boxplot(outlier.size = .3, linewidth = .2, position = position_dodge(.85), width = .8) +
  facet_wrap(~ schemeL, ncol = 1, scales = "free_y") + scale_fill_manual(values = COL) + coord_cartesian(ylim = c(NA, lim)) +
  labs(x = NULL, y = "NRMSE = RMSE / SD (testing set)", fill = "Model") + th
save_fig(pS1, "FigureS1_nrmse", 14, 11)
pS2 <- ggplot(M[design == "CV2"], aes(case, r, fill = method)) + geom_boxplot(outlier.size = .3, linewidth = .2, position = position_dodge(.85), width = .8) +
  facet_wrap(~ schemeL, ncol = 1, scales = "free_y") + scale_fill_manual(values = COL) +
  labs(x = NULL, y = "Pearson correlation (testing set)", fill = "Model") + th
save_fig(pS2, "FigureS2_pearson_boxplots", 14, 11)
pS3 <- ggplot(TA[!is.na(r_TLGBLUP) & !is.na(r_MEGBLUP)], aes(r_MEGBLUP, r_TLGBLUP, colour = crop)) + geom_abline(linetype = 2) + geom_point(size = .9, alpha = .6) +
  facet_wrap(~ schemeL, ncol = 3) + coord_equal() + scale_colour_manual(values = CROPCOL) +
  labs(x = "Pearson correlation, MEGBLUP", y = "Pearson correlation, TLGBLUP", colour = NULL) + th + theme(axis.text.x = element_text(angle = 0))
save_fig(pS3, "FigureS3_TLGBLUP_vs_MEGBLUP", 10, 7.2)
pS4 <- ggplot(TA[!is.na(r_TLGBLUP_H) & !is.na(r_TLGBLUP)], aes(r_TLGBLUP, r_TLGBLUP_H, colour = crop)) + geom_abline(linetype = 2) + geom_point(size = .9, alpha = .6) +
  facet_wrap(~ schemeL, ncol = 3) + coord_equal() + scale_colour_manual(values = CROPCOL) +
  labs(x = "Pearson correlation, TLGBLUP", y = "Pearson correlation, TLGBLUP-H", colour = NULL) + th + theme(axis.text.x = element_text(angle = 0))
save_fig(pS4, "FigureS4_TLGBLUPH_vs_TLGBLUP", 10, 7.2)

## ---- sensitivity analyses (05_sensitivity.R) ---------------------------------------------------------
sf <- file.path(R, "sens_scores.csv")
if (file.exists(sf)) { SE <- fread(sf); SE[, NRMSE := NRMSE_SD]
  REF <- M[method %in% c("GBLUP", "TLGBLUP") & scheme %in% c("S20", "S50"), .(dataset, trait, scheme, rep, method, r)]
  SE2 <- rbind(REF, SE[, .(dataset, trait, scheme, rep, method, r)])
  fwrite(dcast(SE2[, .(r = mean(r)), by = .(dataset, trait, scheme, method)], dataset + trait + scheme ~ method, value.var = "r"),
         file.path(TD, "TableS9_sensitivity_scores_standardization.csv")) }
ef <- file.path(R, "sens_exact.csv")
if (file.exists(ef)) { EX <- fread(ef); EX <- merge(EX[, .(dataset, trait, scheme, rep, r_exact = r)],
    M[method == "MEGBLUP", .(dataset, trait, scheme, rep, r_trunc = r)], by = c("dataset", "trait", "scheme", "rep"))
  fwrite(EX[, .(partitions = .N, r_exact = mean(r_exact), r_truncated = mean(r_trunc), mean_diff = mean(r_exact - r_trunc),
                max_abs_diff = max(abs(r_exact - r_trunc))), by = .(dataset, trait, scheme)], file.path(TD, "TableS10_MEGBLUP_exact_kernels.csv")) }
mf <- file.path(R, "sens_mcmc.csv")
if (file.exists(mf)) { MC <- fread(mf); AG <- fread(file.path(R, "sens_mcmc_agreement.csv"))
  fwrite(MC[, .(minESS = min(ESS), r = r[1], seconds = seconds[1]), by = .(dataset, trait, method, chain)], file.path(TD, "TableS11_mcmc_ess.csv"))
  fwrite(AG[, .(maxPSRF = max(PSRF), cor_pred_std_long = cor_pred_std_long[1], cor_pred_long_long = cor_pred_long_long[1],
                max_abs_diff_pred_sd = max_abs_diff_pred_sd[1]), by = .(dataset, trait, method)], file.path(TD, "TableS12_mcmc_agreement.csv")) }
fwrite(S, file.path(TD, "TableS2_source_models.csv"))
cv1s <- file.path(R, "sources_cv1.csv"); if (file.exists(cv1s)) fwrite(fread(cv1s)[, .(n_mean = mean(n), h2_mean = mean(h2)), by = .(dataset, trait, scheme, source_year)],
                                                                       file.path(TD, "TableS2b_source_models_CV1.csv"))
TM <- fread(file.path(R, "timing.csv"))
fwrite(TM[, .(seconds_mean = mean(seconds), seconds_max = max(seconds)), by = .(dataset, method)][order(dataset, method)], file.path(TD, "TableS4_timing.csv"))
message("figures and tables written")
