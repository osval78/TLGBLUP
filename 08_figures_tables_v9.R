###############################################################################
# 08_figures_tables_v9.R -- figures (PNG 600 dpi, TIFF 600 dpi LZW and PDF) and tables (CSV) of manuscript
# version 9 (G2F maize and EYT wheat data sets; eight models including MEGBLUP-FA fitted by REML).
#   Main text : Figures 1-5, Tables 2-5 (Table 1 = data description, Table 6 = practical guide, both static).
#   Supporting information: Figures S1-S10, Tables S1-S14.
# Data sets with several traits (EYT) are summarized across traits in the main tables and figures; results by
# trait are in the supplementary tables and figures.
# Inputs: ../outputs/results/{metrics,predictions,transfer_weights,shares,sources,timing,angletl_tuning,
#         rotation_metrics,sens_*}.csv(.gz)  and  fa_{metrics,predictions,details,rcor,timing}.csv(.gz)
# Usage:  LC_ALL=C.UTF-8 Rscript 08_figures_tables_v9.R
###############################################################################
suppressPackageStartupMessages({ library(ggplot2); library(data.table); library(grid) })
source("00_functions.R")
KEEP <- c("G2F_2014", "G2F_2015", "G2F_2016", "EYT_4", "EYT_5", "EYT_6", "EYT_7"); DS <- DS[KEEP]
R <- file.path("..", "outputs", "results"); FD <- file.path("..", "outputs", "figures_v9"); TD <- file.path("..", "outputs", "tables_v9")
dir.create(FD, showWarnings = FALSE); dir.create(TD, showWarnings = FALSE)
rd <- function(f) { d <- fread(file.path(R, f)); if ("dataset" %in% names(d)) d <- d[dataset %in% KEEP]; d }
METH <- c("GBLUP", "MEGBLUP-M", "MEGBLUP", "MEGBLUP-FA", "angleTL", "TLGBLUP-F", "TLGBLUP", "TLGBLUP-H")
## intense, clearly distinguishable colours (and one shape per model for grey-scale printing)
COL <- c(GBLUP = "#000000", `MEGBLUP-M` = "#FF8C00", MEGBLUP = "#E8000B", `MEGBLUP-FA` = "#8B00D4",
         angleTL = "#00B140", `TLGBLUP-F` = "#00C8FF", TLGBLUP = "#0028D9", `TLGBLUP-H` = "#FF00B4")
SHP <- c(GBLUP = 16, `MEGBLUP-M` = 15, MEGBLUP = 17, `MEGBLUP-FA` = 18, angleTL = 4, `TLGBLUP-F` = 1, TLGBLUP = 19, `TLGBLUP-H` = 8)
SCH <- c(S10 = "CV2, 10% training", S20 = "CV2, 20% training", S50 = "CV2, 50% training", S80 = "CV2, 80% training",
         CV1_S20 = "CV1, 20% training", CV1_S50 = "CV1, 50% training")
PROP <- c(S10 = 10, S20 = 20, S50 = 50, S80 = 80, CV1_S20 = 20, CV1_S50 = 50)
DSL <- c(G2F_2014 = "G2F-2014", G2F_2015 = "G2F-2015", G2F_2016 = "G2F-2016", EYT_4 = "EYT-4", EYT_5 = "EYT-5", EYT_6 = "EYT-6", EYT_7 = "EYT-7")
CASES <- unlist(lapply(names(DS), function(d) paste(DSL[d], names(DS[[d]]$traits))))
CROPCOL <- c(Maize = "#E8000B", Wheat = "#0028D9")
lab <- function(x) { x[, case := factor(paste(DSL[dataset], trait), CASES)]; x[, dsl := factor(DSL[dataset], DSL)]
  x[, crop := ifelse(grepl("^G2F", dataset), "Maize", "Wheat")]
  if ("scheme" %in% names(x)) { x[, schemeL := factor(SCH[scheme], SCH)]; x[, design := ifelse(grepl("^CV1", scheme), "CV1", "CV2")] }
  if ("method" %in% names(x)) x[, method := factor(method, METH)]; x }
th <- theme_bw(base_size = 11) + theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        strip.background = element_rect(fill = "grey92"), strip.text = element_text(face = "bold"),
        legend.key.width = unit(1.2, "lines"))
save_fig <- function(p, f, w, h) {
  ggsave(file.path(FD, paste0(f, ".png")), p, width = w, height = h, dpi = 600, bg = "white")
  ggsave(file.path(FD, paste0(f, ".tiff")), p, width = w, height = h, dpi = 600, compression = "lzw", bg = "white")
  ggsave(file.path(FD, paste0(f, ".pdf")), p, width = w, height = h, bg = "white") }
SHP["MEGBLUP-FA"] <- 23
gmodel <- list(scale_colour_manual(values = COL, drop = FALSE), scale_shape_manual(values = SHP, drop = FALSE),
               scale_fill_manual(values = COL, drop = FALSE, guide = "none"),
               guides(colour = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)))

## ---- metrics recomputed from the predictions of every testing line ---------------------------------
HIGH <- c(GY = TRUE, PH = FALSE, GM = FALSE, DTHD = FALSE, DTMT = FALSE)
coinc <- function(obs, pred, high) { k <- max(1, round(.2 * length(obs))); sg <- if (high) -1 else 1
  length(intersect(order(sg * obs)[1:k], order(sg * pred)[1:k])) / k }
PR <- rbind(rd("predictions.csv.gz"), rd("fa_predictions.csv.gz"))
M <- PR[, .(n_test = .N, r = cor(pred, obs), NRMSE = sqrt(mean((obs - pred)^2)) / sd(obs),
            NRMSE_mean = sqrt(mean((obs - pred)^2)) / mean(obs), CI20 = coinc(obs, pred, HIGH[trait[1]])),
        by = .(dataset, trait, scheme, rep, method)]
NTR <- rbind(rd("metrics.csv"), rd("fa_metrics.csv"), fill = TRUE)[, .(n_train = n_train[1]), by = .(dataset, trait, scheme, rep, method)]
M <- merge(M, NTR, by = c("dataset", "trait", "scheme", "rep", "method")); rm(PR); invisible(gc())
M <- lab(M[method %in% METH])
fwrite(M[order(dataset, trait, scheme, rep, method), .(dataset, trait, scheme, rep, method, n_train, n_test, r, NRMSE, NRMSE_mean, CI20)],
       file.path(TD, "TableS_all_metrics_by_partition.csv"))

## ---- means by data set x trait (Table S3) and by data set across traits (Tables 2-3) ----------------
t2 <- M[, .(r_mean = mean(r), r_sd = sd(r), NRMSE_mean = mean(NRMSE), CI20_mean = mean(CI20), parts = .N),
        by = .(dataset, dsl, case, crop, trait, scheme, design, method)][order(case, scheme, method)]
fwrite(t2[, -c("case", "dsl")], file.path(TD, "TableS3_r_by_dataset_trait.csv"))
tD <- t2[, .(r = mean(r_mean), NRMSE = mean(NRMSE_mean), CI20 = mean(CI20_mean), traits = .N), by = .(dataset, dsl, crop, scheme, design, method)]
tA <- rbind(t2[, .(dataset = "Maize (3)", r = mean(r_mean), NRMSE = mean(NRMSE_mean), CI20 = mean(CI20_mean), traits = .N), by = .(crop, scheme, design, method)][crop == "Maize"],
            t2[, .(dataset = "Wheat (16)", r = mean(r_mean), NRMSE = mean(NRMSE_mean), CI20 = mean(CI20_mean), traits = .N), by = .(crop, scheme, design, method)][crop == "Wheat"],
            t2[, .(crop = "All", dataset = "All (19)", r = mean(r_mean), NRMSE = mean(NRMSE_mean), CI20 = mean(CI20_mean), traits = .N), by = .(scheme, design, method)])
fwrite(dcast(tD, dataset + scheme ~ method, value.var = "r"), file.path(TD, "Table2_3_r_by_dataset_across_traits.csv"))
fwrite(dcast(tA, dataset + scheme ~ method, value.var = "r"), file.path(TD, "Table2_3_r_overall.csv"))
fwrite(rbind(tD[, .(dataset = as.character(dsl), crop, scheme, method, r, NRMSE, CI20)], tA[, .(dataset, crop, scheme, method, r, NRMSE, CI20)]),
       file.path(TD, "TableS4_r_NRMSE_CI20_summary.csv"))

## ---- paired differences vs GBLUP: corrected resampled t test (Nadeau & Bengio, 2003) ----------------
nb <- function(d, ntr, nte) { J <- length(d); m <- mean(d); v <- (1 / J + nte / ntr) * var(d); se <- sqrt(v); q <- qt(.975, J - 1)
  c(mean = m, lo = m - q * se, hi = m + q * se, p = 2 * pt(-abs(m / se), J - 1), lo_naive = m - q * sd(d) / sqrt(J), hi_naive = m + q * sd(d) / sqrt(J)) }
B <- M[method == "GBLUP", .(dataset, trait, scheme, rep, r0 = r)]
DD <- merge(M[method != "GBLUP"], B, by = c("dataset", "trait", "scheme", "rep")); DD[, dr := r - r0]
G5 <- DD[, { a <- nb(dr, n_train[1], n_test[1]); .(dr = a[["mean"]], dr_lo = a[["lo"]], dr_hi = a[["hi"]], p_r = a[["p"]],
  dr_lo_naive = a[["lo_naive"]], dr_hi_naive = a[["hi_naive"]], win_r = mean(dr > 0)) },
  by = .(dataset, dsl, case, crop, trait, scheme, schemeL, design, method)]
G5[, p_r_holm := p.adjust(p_r, "holm"), by = .(scheme, method)]
G5[, sig := fifelse(dr_lo > 0, "gain", fifelse(dr_hi < 0, "loss", "ns"))]
G5[, sig_naive := fifelse(dr_lo_naive > 0, "gain", fifelse(dr_hi_naive < 0, "loss", "ns"))]
fwrite(G5[order(case, scheme, method), -c("schemeL", "case", "dsl")], file.path(TD, "TableS5_differences_vs_GBLUP_by_trait.csv"))
SIG <- G5[, .(cases = .N, mean_dr = mean(dr), wins_80pct = sum(win_r >= .8), gains_naive = sum(sig_naive == "gain"), losses_naive = sum(sig_naive == "loss"),
              gains_corrected = sum(sig == "gain"), losses_corrected = sum(sig == "loss"),
              gains_holm = sum(p_r_holm < .05 & dr > 0), losses_holm = sum(p_r_holm < .05 & dr < 0)), by = .(scheme, method)][order(scheme, method)]
fwrite(SIG, file.path(TD, "Table4_significance_counts.csv"))
## data-set level (EYT: mean over the four traits within each partition; the traits share the partitions)
DDs <- DD[, .(dr = mean(dr), n_train = n_train[1], n_test = n_test[1]), by = .(dataset, dsl, crop, scheme, schemeL, design, method, rep)]
GD <- DDs[, { a <- nb(dr, n_train[1], n_test[1]); .(dr = a[["mean"]], dr_lo = a[["lo"]], dr_hi = a[["hi"]], p_r = a[["p"]], win_r = mean(dr > 0)) },
          by = .(dataset, dsl, crop, scheme, schemeL, design, method)]
fwrite(GD[order(dataset, scheme, method), -"schemeL"], file.path(TD, "TableS6_differences_vs_GBLUP_by_dataset.csv"))
## TLGBLUP (and TLGBLUP-H) versus MEGBLUP-FA, paired by partition
pair_test <- function(m1, m2) {
  X <- merge(M[method == m1, .(dataset, trait, scheme, rep, r1 = r, n_train, n_test)], M[method == m2, .(dataset, trait, scheme, rep, r2 = r)],
             by = c("dataset", "trait", "scheme", "rep"))
  X[, { a <- nb(r1 - r2, n_train[1], n_test[1]); .(d = a[["mean"]], lo = a[["lo"]], hi = a[["hi"]], win = mean(r1 > r2)) }, by = .(dataset, trait, scheme)][,
    .(model = m1, versus = m2, cases = .N, mean_diff = mean(d), better_sig = sum(lo > 0), worse_sig = sum(hi < 0), wins_80pct = sum(win >= .8), losses_80pct = sum(win <= .2)),
    by = .(scheme)] }
PT <- rbindlist(list(pair_test("TLGBLUP", "MEGBLUP-FA"), pair_test("TLGBLUP-H", "MEGBLUP-FA"), pair_test("MEGBLUP-FA", "MEGBLUP"), pair_test("TLGBLUP", "angleTL")))
PT2 <- rbindlist(lapply(list(c("TLGBLUP", "MEGBLUP-FA"), c("TLGBLUP-H", "MEGBLUP-FA"), c("MEGBLUP-FA", "MEGBLUP"), c("TLGBLUP", "angleTL")), function(m) {
  X <- merge(M[method == m[1], .(dataset, crop, trait, scheme, rep, r1 = r)], M[method == m[2], .(dataset, trait, scheme, rep, r2 = r)], by = c("dataset", "trait", "scheme", "rep"))
  X[, .(model = m[1], versus = m[2], mean_diff = mean(r1 - r2)), by = .(crop, scheme)] }))
fwrite(PT, file.path(TD, "TableS7_pairwise_TLGBLUP_MEGBLUPFA.csv")); fwrite(PT2, file.path(TD, "TableS7b_pairwise_by_crop.csv"))

## ---- Figure 1: workflow ---------------------------------------------------------------------------------
for (fmt in c("png", "tiff", "pdf")) {
  fn <- file.path(FD, paste0("Figure1_TLGBLUP_workflow.", fmt))
  if (fmt == "png") png(fn, width = 10, height = 5.4, units = "in", res = 600) else if (fmt == "tiff") tiff(fn, width = 10, height = 5.4, units = "in", res = 600, compression = "lzw") else pdf(fn, width = 10, height = 5.4)
  grid.newpage(); pushViewport(viewport(xscale = c(0, 100), yscale = c(0, 54)))
  box <- function(x, y, w, h, fill, col = "grey10") grid.roundrect(unit(x, "native"), unit(y, "native"), unit(w, "native"), unit(h, "native"),
    r = unit(2, "mm"), gp = gpar(fill = fill, col = col, lwd = 1.4))
  txt <- function(x, y, lab, cex = .78, col = "black", ...) grid.text(lab, unit(x, "native"), unit(y, "native"), gp = gpar(cex = cex, col = col, ...))
  arr <- function(x0, y0, x1, y1) grid.lines(unit(c(x0, x1), "native"), unit(c(y0, y1), "native"),
    arrow = arrow(length = unit(2.2, "mm"), type = "closed"), gp = gpar(fill = "grey10", col = "grey10", lwd = 1.6))
  ys <- c(47, 36, 25, 14)
  for (k in 1:4) {
    box(11, ys[k], 19, 8, "#FFB000"); txt(11, ys[k] + 1.6, if (k < 4) paste("Source environment", k) else "Source environment K", fontface = "bold")
    txt(11, ys[k] - 1.8, if (k < 4) bquote(bold(y)[.(k)] == mu[.(k)] * bold("1") + bold(g)[.(k)] + bold(e)[.(k)]) else bquote(bold(y)[K] == mu[K] * bold("1") + bold(g)[K] + bold(e)[K]))
    box(36, ys[k], 19, 8, "#FFE066"); txt(36, ys[k] + 1.6, "Marker effects (summary)", fontface = "bold")
    txt(36, ys[k] - 1.8, if (k < 4) bquote(hat(bold(w))[.(k)] == bold(Z)[.(k)]^T * bold(G)[.(k)]^"+" * hat(bold(g))[.(k)] / p) else bquote(hat(bold(w))[K] == bold(Z)[K]^T * bold(G)[K]^"+" * hat(bold(g))[K] / p))
    arr(20.5, ys[k], 26.5, ys[k]); arr(45.5, ys[k], 54.5, 30.5) }
  box(63, 30.5, 17, 12, "#5BC0FF"); txt(63, 33.5, "Source genomic scores", fontface = "bold")
  txt(63, 29, bquote(bold(s)[k] == bold(Z)[T] * hat(bold(w))[k]), cex = .9); txt(63, 25.8, "(all target lines)", cex = .7)
  box(87, 43, 24, 16, "#1E64FF"); txt(87, 48.5, "Target environment (training lines)", fontface = "bold", col = "white")
  txt(87, 44.5, bquote("TLGBLUP: " * bold(y) == mu * bold("1") + bold(S) * bold(c) + bold(g) + bold(e)), cex = .74, col = "white")
  txt(87, 40.8, bquote("TLGBLUP-H: " * "... " + bold(h) * ",  " * bold(h) %~% N(bold(0), sigma[h]^2 * bold(K)[S] * " \u2218 " * bold(G)[T])), cex = .66, col = "white")
  txt(87, 37.3, bquote(c[k] %~% N(0, sigma[c]^2) * ",   " * bold(g) %~% N(bold(0), sigma[g]^2 * bold(G)[T])), cex = .7, col = "white")
  box(87, 15, 24, 10, "#00B140"); txt(87, 17.5, "Prediction of testing lines", fontface = "bold", col = "white")
  txt(87, 13.5, bquote(hat(y)[i] == hat(mu) + sum(hat(c)[k] * s[ik], k == 1, K) + hat(g)[i] + hat(h)[i]), col = "white")
  arr(71.5, 32.5, 75, 37); arr(87, 35, 87, 20)
  txt(50, 3, "Only the p-dimensional vectors of marker effects leave the source environments (summary-level, multi-source transfer)", cex = .72, fontface = "italic")
  txt(50, 0.8, expression("In TLGBLUP-H, " * bold(K)[S] == tilde(bold(S)) * tilde(bold(S))^T / K * " and \u2218 is the Hadamard product (interaction between " * bold(S) * bold(c) * " and " * bold(g) * ")"), cex = .66)
  dev.off() }

## ---- Figure 2: prediction schemes -----------------------------------------------------------------------
sc <- CJ(line = 1:12, env = c("Source 1", "Source 2", "Source 3", "Target"), panel = c("CV2: testing lines have records in the sources", "CV1: testing lines removed from the sources"))
testl <- c(2, 5, 7, 9, 10, 12)
sc[, status := fifelse(env == "Target", fifelse(line %in% testl, "Target testing (predicted)", "Target training"),
                       fifelse(grepl("^CV1", panel) & line %in% testl, "Removed from sources", "Source record"))]
sc[, panel := factor(panel, c("CV2: testing lines have records in the sources", "CV1: testing lines removed from the sources"))]; sc[, status := factor(status, c("Source record", "Target training", "Target testing (predicted)", "Removed from sources"))]
pF2 <- ggplot(sc, aes(env, factor(line, 12:1), fill = status)) + geom_tile(colour = "white", linewidth = 1) + facet_wrap(~ panel) +
  scale_fill_manual(values = c("Source record" = "#0028D9", "Target training" = "#5BC0FF", "Target testing (predicted)" = "#FF8C00", "Removed from sources" = "#D9D9D9")) +
  scale_y_discrete(labels = paste("Line", 12:1)) + labs(x = "Environment", y = NULL, fill = NULL) + th +
  theme(panel.grid = element_blank(), panel.border = element_blank(), strip.background = element_blank(), strip.text = element_text(size = 11))
save_fig(pF2, "Figure2_prediction_schemes", 9, 5)

## ---- Figure 3: accuracy vs training proportion, by data set (across traits) --------------------------------
A3 <- rbind(tD[design == "CV2", .(panel = as.character(dsl), scheme, method, r)], tA[design == "CV2" & dataset == "All (19)", .(panel = "All 19 combinations", scheme, method, r)])
A3[, panel := factor(panel, c(DSL, "All 19 combinations"))]; A3[, prop := PROP[scheme]]
p3 <- ggplot(A3, aes(prop, r, colour = method, shape = method, fill = method)) + geom_line(linewidth = .8) + geom_point(size = 2.4, stroke = 1) +
  facet_wrap(~ panel, ncol = 4, scales = "free_y") + scale_x_continuous(breaks = c(10, 20, 50, 80)) + gmodel +
  labs(x = "Target lines used for training (%)", y = "Mean Pearson correlation (testing lines)", colour = "Model", shape = "Model") + th +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = 3, linewidth = 1)))
save_fig(p3, "Figure3_accuracy_vs_training", 12, 7.2)

## ---- Figure 4: difference with GBLUP, CV2 vs CV1, by data set (across traits; corrected 95% CIs) ------------
F4 <- GD[scheme %in% c("S20", "S50", "CV1_S20", "CV1_S50")]
F4[, panel := factor(ifelse(design == "CV2", "CV2: testing lines evaluated in the sources", "CV1: new lines (removed from the sources)"),
                     c("CV2: testing lines evaluated in the sources", "CV1: new lines (removed from the sources)"))]
F4[, train := factor(paste0(PROP[scheme], "% training"), c("20% training", "50% training"))]
p4 <- ggplot(F4, aes(dsl, dr, colour = method, shape = method, fill = method)) + geom_hline(yintercept = 0, linetype = 2) +
  geom_pointrange(aes(ymin = dr_lo, ymax = dr_hi), position = position_dodge(.8), size = .35, linewidth = .5, stroke = .9) +
  facet_grid(panel ~ train, scales = "free_y") + gmodel + labs(x = NULL, y = "Difference in Pearson correlation vs GBLUP\n(mean and corrected 95% CI)", colour = "Model", shape = "Model") + th +
  guides(colour = guide_legend(nrow = 1, override.aes = list(size = .7, linewidth = .7)))
p4 <- suppressMessages(p4 + scale_colour_manual(values = COL[-1]) + scale_shape_manual(values = SHP[-1]) + scale_fill_manual(values = COL[-1], guide = "none"))
save_fig(p4, "Figure4_gain_vs_GBLUP_CV2_CV1", 13, 8.5)

## ---- Figure 5 and Table 5: transfer weights and transfer share --------------------------------------------------
W <- lab(rd("transfer_weights.csv")); SH <- lab(rd("shares.csv")); W[, srcL := sub("_20[0-9]{2}$", "", source_year)]
p5 <- ggplot(W[method == "TLGBLUP" & scheme %in% c("S20", "S50", "S80")], aes(srcL, c, fill = schemeL)) + geom_hline(yintercept = 0, linetype = 2) +
  geom_boxplot(outlier.size = .25, linewidth = .25) + facet_wrap(~ case, scales = "free", ncol = 5) +
  scale_fill_manual(values = c("#FFE066", "#FF8C00", "#E8000B")) + labs(x = "Source environment", y = expression("Transfer weight " * hat(c)[k] * " (TLGBLUP)"), fill = NULL) + th +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6.5))
save_fig(p5, "Figure5_transfer_weights", 12, 13)
fwrite(W[, .(c_mean = mean(c), c_sd = sd(c)), by = .(dataset, trait, scheme, method, source_year)], file.path(TD, "TableS9_transfer_weights.csv"))
TSH <- SH[method %in% c("TLGBLUP", "TLGBLUP-H"), .(share = mean(share), shareH = mean(share_H)), by = .(dataset, dsl, case, trait, scheme, method)]
fwrite(TSH[, -c("case", "dsl")], file.path(TD, "TableS8_transfer_share_by_trait.csv"))
T5 <- rbind(TSH[method == "TLGBLUP" & scheme %in% c("S20", "S50", "CV1_S50"), .(v = share, col = paste("TLGBLUP", scheme)), by = .(dataset, dsl, trait)],
            TSH[method == "TLGBLUP-H" & scheme %in% c("S20", "S50"), .(v = shareH, col = paste("H", scheme)), by = .(dataset, dsl, trait)])
fwrite(T5[, .(mean = mean(v), min = min(v), max = max(v)), by = .(dataset, dsl, col)][, -"dsl"], file.path(TD, "Table5_transfer_share_by_dataset.csv"))

## ---- MEGBLUP-FA: convergence, number of factors and genetic correlations target vs sources ----------------------
FAd <- rd("fa_details.csv"); FAr <- rd("fa_rcor.csv")
fwrite(FAd[, .(k = k[1], J = J[1], fits = .N, converged = mean(converged), iter_mean = mean(iter), rg_mean = mean(rg_mean), rg_min = mean(rg_min), rg_max = mean(rg_max)),
           by = .(dataset, trait, design = ifelse(grepl("^CV1", scheme), "CV1", "CV2"))][order(dataset, trait, design)], file.path(TD, "TableS12_MEGBLUP_FA_details.csv"))
FR <- lab(FAr[!grepl("^CV1", scheme), .(rg = mean(rg)), by = .(dataset, trait, source, scheme)])
FR[, srcL := sub("_20[0-9]{2}$", "", source)]
pS10 <- ggplot(FR[scheme %in% c("S20", "S50", "S80")], aes(srcL, rg, colour = schemeL)) + geom_hline(yintercept = 0, linetype = 2) +
  geom_point(size = 1.8, position = position_dodge(.6)) + facet_wrap(~ case, scales = "free_x", ncol = 5) + coord_cartesian(ylim = c(-1, 1)) +
  scale_colour_manual(values = c("#FFB000", "#FF4500", "#8B00D4")) + labs(x = "Source environment", y = "Genetic correlation with the target (MEGBLUP-FA, REML)", colour = NULL) + th +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 6.5))
save_fig(pS10, "FigureS10_FA_genetic_correlations", 12, 12)

## ---- Supplementary figures S1-S9 ------------------------------------------------------------------------------------
YY <- list(); EN <- list()
for (dn in names(DS)) { D <- load_dataset(dn); y <- as.data.table(D$Y)
  for (tr in D$traits) { yt <- y[!is.na(get(tr))]
    YY[[paste(dn, tr)]] <- data.table(dataset = dn, trait = tr, Env = yt$Env, value = yt[[tr]], role = ifelse(yt$Env == D$target, "Target", "Source")) }
  EN[[dn]] <- y[, .(lines = .N), by = Env][, `:=`(dataset = dn, role = ifelse(Env == D$target, "target", "source"))] }
rm(D); invisible(gc()); YY <- lab(rbindlist(YY)); fwrite(rbindlist(EN)[, .(dataset, Env, lines, role)], file.path(TD, "TableS1_environments.csv"))
YM <- YY[crop == "Maize"]; YM[, EnvS := sub("_20[0-9]{2}$", "", Env)]
pS1 <- ggplot(YM, aes(EnvS, value, fill = role)) + geom_boxplot(outlier.size = .3, linewidth = .3) + facet_wrap(~ case, scales = "free", ncol = 1) +
  scale_fill_manual(values = c(Source = "#FFB000", Target = "#0028D9")) + labs(x = "Environment (location)", y = "Grain moisture (%)", fill = NULL) + th +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 7))
save_fig(pS1, "FigureS1_phenotypes_maize", 10, 9)
pS2 <- ggplot(YY[crop == "Wheat"], aes(Env, value, fill = role)) + geom_boxplot(outlier.size = .3, linewidth = .3) + facet_wrap(~ case, scales = "free", ncol = 4) +
  scale_fill_manual(values = c(Source = "#FFB000", Target = "#0028D9")) + labs(x = "Environment (management condition)", y = "Phenotype", fill = NULL) + th +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7.5))
save_fig(pS2, "FigureS2_phenotypes_wheat", 10, 9)
A4 <- t2[design == "CV2"]; A4[, prop := PROP[scheme]]
pS3 <- ggplot(A4, aes(prop, r_mean, colour = method, shape = method, fill = method)) + geom_line(linewidth = .6) + geom_point(size = 1.6, stroke = .8) +
  facet_wrap(~ case, ncol = 5, scales = "free_y") + scale_x_continuous(breaks = c(10, 20, 50, 80)) + gmodel +
  labs(x = "Target lines used for training (%)", y = "Mean Pearson correlation (10 or 5 partitions)", colour = "Model", shape = "Model") + th
save_fig(pS3, "FigureS3_accuracy_vs_training_by_trait", 13, 12)
pS4 <- ggplot(G5[design == "CV2"], aes(case, dr, colour = method, shape = method, fill = method)) + geom_hline(yintercept = 0, linetype = 2) +
  geom_pointrange(aes(ymin = dr_lo, ymax = dr_hi), position = position_dodge(.8), size = .22, linewidth = .4, stroke = .7) +
  facet_wrap(~ schemeL, ncol = 1, scales = "free_y") + gmodel + labs(x = NULL, y = "Difference in Pearson correlation vs GBLUP\n(mean and corrected 95% CI)", colour = "Model", shape = "Model") +
  th + theme(axis.text.x = element_text(angle = 40, hjust = 1)) + guides(colour = guide_legend(nrow = 1, override.aes = list(size = .6, linewidth = .7)))
pS4 <- suppressMessages(pS4 + scale_colour_manual(values = COL[-1]) + scale_shape_manual(values = SHP[-1]) + scale_fill_manual(values = COL[-1], guide = "none"))
save_fig(pS4, "FigureS4_gain_vs_GBLUP_by_trait", 14, 12)
C6 <- G5[scheme %in% c("S20", "S50", "CV1_S20", "CV1_S50")]; C6[, train := paste0(PROP[scheme], "% training")]
pS5 <- ggplot(C6, aes(case, dr, colour = method, shape = design)) + geom_hline(yintercept = 0, linetype = 2) + geom_point(position = position_dodge(.8), size = 2) +
  facet_wrap(~ train, ncol = 1) + scale_colour_manual(values = COL) + scale_shape_manual(values = c(CV1 = 17, CV2 = 16)) +
  labs(x = NULL, y = "Mean difference in Pearson correlation vs GBLUP", colour = "Model", shape = "Testing lines in the sources") + th +
  theme(axis.text.x = element_text(angle = 40, hjust = 1)) + guides(colour = guide_legend(nrow = 1, override.aes = list(size = 3)))
save_fig(pS5, "FigureS5_CV1_vs_CV2_by_trait", 14, 8)
TA <- dcast(M, dataset + crop + trait + scheme + schemeL + design + rep ~ method, value.var = "r"); setnames(TA, gsub("-", "_", names(TA)))
pS6 <- ggplot(TA, aes(angleTL, TLGBLUP, colour = crop)) + geom_abline(linetype = 2) + geom_point(size = 1.1, alpha = .7) + facet_wrap(~ schemeL, ncol = 3) +
  coord_equal() + scale_colour_manual(values = CROPCOL) + labs(x = "Pearson correlation, angleTL (cross-validated)", y = "Pearson correlation, TLGBLUP", colour = NULL) + th
save_fig(pS6, "FigureS6_TLGBLUP_vs_angleTL", 10, 7.5)
pS7 <- ggplot(SH[method == "TLGBLUP"], aes(case, share, fill = schemeL)) + geom_boxplot(outlier.size = .25, linewidth = .25, position = position_dodge(.85)) +
  scale_fill_manual(values = c("#FFE066", "#FFB000", "#FF4500", "#B30000", "#5BC0FF", "#0028D9")) + coord_cartesian(ylim = c(0, 1)) +
  labs(x = NULL, y = "Share of the target genetic variance\nexplained by the transferred scores (TLGBLUP)", fill = NULL) + th + theme(axis.text.x = element_text(angle = 40, hjust = 1))
save_fig(pS7, "FigureS7_transfer_share", 14, 6)
rf <- file.path(R, "rotation_metrics.csv")
RO <- fread(rf); RO[, method := factor(method, METH)]
RA <- RO[, .(r = mean(r), share = mean(share)), by = .(dataset, target, method)]
RB <- merge(RA[method != "GBLUP"], RA[method == "GBLUP", .(dataset, target, r0 = r)], by = c("dataset", "target")); RB[, dr := r - r0]
RS <- RB[, .(targets = .N, mean_dr = mean(dr), median_dr = median(dr), p_t = t.test(dr)$p.value, p_wilcox = suppressWarnings(wilcox.test(dr)$p.value),
             targets_better = mean(dr > 0)), by = .(method)][order(method)]
fwrite(merge(RS, RA[, .(mean_r = mean(r)), by = method], by = "method"), file.path(TD, "TableS10_G2F_rotation_summary.csv"))
RB[, yearL := sub("_", "-", dataset)]
pS8 <- ggplot(RB, aes(method, dr, fill = method)) + geom_hline(yintercept = 0, linetype = 2) + geom_boxplot(outlier.size = .5, linewidth = .3) +
  geom_jitter(width = .12, size = .6, alpha = .6) + facet_wrap(~ yearL, nrow = 1) + scale_fill_manual(values = COL, guide = "none") +
  labs(x = NULL, y = "Difference in Pearson correlation vs GBLUP\n(mean over partitions per target location)") + th + theme(axis.text.x = element_text(angle = 40, hjust = 1))
save_fig(pS8, "FigureS8_G2F_rotation", 11, 5)
pS9 <- ggplot(TA[!is.na(MEGBLUP_FA)], aes(MEGBLUP_FA, TLGBLUP, colour = crop)) + geom_abline(linetype = 2) + geom_point(size = 1.1, alpha = .7) + facet_wrap(~ schemeL, ncol = 3) +
  coord_equal() + scale_colour_manual(values = CROPCOL) + labs(x = "Pearson correlation, MEGBLUP-FA (REML)", y = "Pearson correlation, TLGBLUP", colour = NULL) + th
save_fig(pS9, "FigureS9_TLGBLUP_vs_MEGBLUP_FA", 10, 7.5)

## ---- other supplementary tables ------------------------------------------------------------------------------------
S <- rd("sources.csv"); fwrite(S[, .(environments = .N, n_min = min(n), n_max = max(n), h2_min = min(h2), h2_max = max(h2), h2_mean = mean(h2)), by = .(dataset, trait)],
                             file.path(TD, "TableS2_source_heritabilities.csv"))
AT <- rd("angletl_tuning.csv"); fwrite(AT[, .(c_mean = mean(c), share_mean_source = mean(source == "mean")), by = .(dataset, trait, scheme)], file.path(TD, "TableS13_angleTL_tuning.csv"))
TM <- rbind(rd("timing.csv"), rd("fa_timing.csv")); TM[, crop := ifelse(grepl("^G2F", dataset), "Maize", "Wheat")]
fwrite(TM[method %in% METH, .(seconds_mean = mean(seconds), seconds_min_dataset = NA_real_), by = .(crop, dataset, method)][order(crop, dataset, method)][, -"seconds_min_dataset"],
       file.path(TD, "TableS11_timing.csv"))
SE <- rd("sens_scores.csv"); REF <- M[method %in% c("GBLUP", "TLGBLUP") & scheme %in% c("S20", "S50") & crop == "Maize", .(dataset, trait, scheme, rep, method, r)]
SE2 <- rbind(REF, SE[, .(dataset, trait, scheme, rep, method, r)])
fwrite(dcast(SE2[, .(r = mean(r)), by = .(dataset, trait, scheme, method)], dataset + trait + scheme ~ method, value.var = "r"), file.path(TD, "TableS14_sensitivity.csv"))
EX <- merge(rd("sens_exact.csv")[, .(dataset, trait, scheme, rep, r_exact = r)], M[method == "MEGBLUP", .(dataset, trait, scheme, rep, r_trunc = r)], by = c("dataset", "trait", "scheme", "rep"))
fwrite(EX[, .(partitions = .N, r_exact = mean(r_exact), r_truncated = mean(r_trunc), mean_diff = mean(r_exact - r_trunc)), by = .(dataset, trait, scheme)], file.path(TD, "TableS14b_MEGBLUP_exact_kernels.csv"))
MC <- rd("sens_mcmc.csv"); AG <- rd("sens_mcmc_agreement.csv")
fwrite(MC[, .(minESS = min(ESS)), by = .(dataset, trait, method, chain)], file.path(TD, "TableS14c_mcmc_ess.csv"))
fwrite(AG[, .(maxPSRF = max(PSRF), cor_pred_std_long = cor_pred_std_long[1]), by = .(dataset, trait, method)], file.path(TD, "TableS14d_mcmc_agreement.csv"))
message("v9 figures and tables written")
