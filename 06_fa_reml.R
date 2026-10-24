###############################################################################
# 06_fa_reml.R -- MEGBLUP-FA: multi-environment genomic model with a factor-analytic
# (FA) genetic covariance between environments, fitted by restricted maximum likelihood.
#
#   y_ij = mu_j + g_ij + e_ij,   vec(g) ~ N(0, (L L' + Psi) (x) G),   e_ij ~ N(0, s2_j)
#
# L (J x k) are the loadings of the k common factors, Psi = diag(psi_j) the environment-specific
# genetic variances and s2_j heterogeneous residual variances (one per environment); mu_j are fixed
# environment means.  Identifiability: the upper triangle of L is set to zero (L[1, 2:k] = 0).
#
# Algorithm (exact REML; average-information (AI) algorithm):
#  * Common factors: f_k = F a_k, a_k ~ N(0, I_r), F = U D^{1/2} with G_L = U D U' (lines with records),
#    so F F' = G_L.  Environment-specific terms are integrated out: within environment j,
#    B_j = psi_j G_j + s2_j I = Q_j diag(psi_j e_j + s2_j) Q_j' (eigendecomposition of G_j computed once).
#  * The mixed-model system for theta = (mu, a_1, ..., a_k) has dimension J + k r; with flat priors on mu
#    its solution gives P y, log|V| + log|X'V^-1 X| and y'Py, i.e. the exact REML log-likelihood.
#  * Score: Fisher identity (expected complete-data score under the posterior of theta).
#  * AI matrix: 0.5 * W' P W with W = [dV/dtheta_a P y]; Newton steps with step halving.
#  * log-scale parameters for psi_j and s2_j (both kept positive).
# Environments whose records cover all lines of G_L use the eigenbasis of G_L itself (diagonal blocks).
###############################################################################

fa_reml <- function(dd, y, G, target, test_lines, k = 2, maxit = 100, tol = 1e-4, verbose = FALSE, init = NULL) {
  ## dd: data.frame with Env, Line (one row per record); y: phenotypes with NA for masked records
  ## G : genomic relationship matrix with rownames = lines;  target: name of the target environment
  ## test_lines: lines of the target environment to be predicted
  obs <- !is.na(y); dd <- dd[obs, ]; y <- y[obs]
  envs <- sort(unique(dd$Env)); J <- length(envs)
  if (!(target %in% envs)) stop("target environment without records")
  L <- rownames(G)[rownames(G) %in% unique(dd$Line)]; NL <- length(L)
  GL <- G[L, L]; eg <- eigen(GL, symmetric = TRUE); eg$values <- pmax(eg$values, 0)
  keep <- eg$values > 1e-8 * max(eg$values); r <- sum(keep)
  U <- eg$vectors[, keep, drop = FALSE]; d <- eg$values[keep]; Fm <- sweep(U, 2, sqrt(d), "*")   # F F' = G_L
  ## environment blocks (records ordered as in L)
  E <- lapply(envs, function(en) {
    i <- which(dd$Env == en); i <- i[order(match(dd$Line[i], L))]; li <- match(dd$Line[i], L)
    full <- length(li) == NL
    if (full) { Q <- eg$vectors; ev <- eg$values; P <- NULL } else {
      e <- eigen(GL[li, li, drop = FALSE], symmetric = TRUE); Q <- e$vectors; ev <- pmax(e$values, 0)
      P <- crossprod(Fm[li, , drop = FALSE], Q) }                      # r x n_j
    yy <- y[i]; list(env = en, li = li, n = length(li), full = full, Q = Q, ev = ev, P = P,
                     q1 = drop(crossprod(Q, rep(1, length(li)))), qy = drop(crossprod(Q, yy)), v = var(yy)) })
  ## P' x  and  P x  (shortcut for environments covering all lines: P = [diag(sqrt(d)) 0])
  Pt <- function(b, x) if (b$full) c(sqrt(d) * x, rep(0, b$n - r)) else drop(crossprod(b$P, x))
  Px <- function(b, x) if (b$full) sqrt(d) * x[seq_len(r)] else drop(b$P %*% x)
  PxM <- function(b, X) if (b$full) sqrt(d) * X[seq_len(r), , drop = FALSE] else b$P %*% X
  ## parameters
  lam_free <- matrix(TRUE, J, k); if (k > 1) for (q in 2:k) lam_free[seq_len(q - 1), q] <- FALSE
  v <- sapply(E, `[[`, "v")
  Lam <- matrix(0, J, k); Lam[, 1] <- sqrt(.45 * v)
  if (k > 1) for (q in 2:k) Lam[, q] <- .2 * sqrt(v) * rep(c(1, -1), length.out = J)
  Lam[!lam_free] <- 0; psi <- .10 * v; s2 <- .40 * v
  if (!is.null(init)) { Lam <- init$Lam; psi <- init$psi; s2 <- init$s2 }
  lo_psi <- log(1e-6 * v); lo_s2 <- log(1e-4 * v)
  pack <- function(Lam, psi, s2) c(Lam[lam_free], log(psi), log(s2))
  unpack <- function(th) { Lm <- matrix(0, J, k); Lm[lam_free] <- th[seq_len(sum(lam_free))]
    nl <- sum(lam_free); list(Lam = Lm, psi = exp(th[nl + 1:J]), s2 = exp(th[nl + J + 1:J])) }
  nl <- sum(lam_free); npar <- nl + 2 * J
  idx_mu <- seq_len(J); idx_a <- function(q) J + (q - 1) * r + seq_len(r); nth <- J + k * r

  core <- function(th, full = FALSE) {
    p <- unpack(th); Lm <- p$Lam
    Om <- matrix(0, nth, nth); b <- numeric(nth); ldB <- 0; syy <- 0
    for (j in seq_len(J)) { B <- E[[j]]; S <- 1 / (p$psi[j] * B$ev + p$s2[j]); ldB <- ldB - sum(log(S))
      Sq1 <- S * B$q1; Sqy <- S * B$qy; syy <- syy + sum(B$qy * Sqy)
      Om[j, j] <- sum(B$q1 * Sq1); b[j] <- sum(B$q1 * Sqy)
      u1 <- Px(B, Sq1); uy <- Px(B, Sqy)
      if (B$full) Mdiag <- d * S[seq_len(r)] else M <- tcrossprod(sweep(B$P, 2, sqrt(S), "*"))
      E[[j]]$S <<- S; E[[j]]$u1 <<- u1; E[[j]]$uy <<- uy
      if (B$full) { E[[j]]$Mdiag <<- Mdiag; E[[j]]$M <<- NULL } else E[[j]]$M <<- M
      for (q in seq_len(k)) { iq <- idx_a(q)
        Om[j, iq] <- Lm[j, q] * u1; b[iq] <- b[iq] + Lm[j, q] * uy
        for (s in q:k) { is <- idx_a(s); w <- Lm[j, q] * Lm[j, s]; if (w == 0) next
          if (B$full) { ii <- cbind(iq, is); Om[ii] <- Om[ii] + w * Mdiag } else Om[iq, is] <- Om[iq, is] + w * M } } }
    for (q in seq_len(k)) { iq <- idx_a(q); Om[cbind(iq, iq)] <- Om[cbind(iq, iq)] + 1
      if (q < k) for (s in (q + 1):k) Om[idx_a(s), iq] <- t(Om[iq, idx_a(s)]) }
    Om[idx_a(1)[1]:nth, idx_mu] <- t(Om[idx_mu, idx_a(1)[1]:nth])
    R <- tryCatch(chol(Om), error = function(e) NULL); if (is.null(R)) return(list(ll = -Inf))
    m <- backsolve(R, forwardsolve(t(R), b))
    ll <- -0.5 * (ldB + 2 * sum(log(diag(R))) + syy - sum(b * m))
    list(ll = ll, R = R, m = m, p = p)
  }

  score_ai <- function(cr) {
    p <- cr$p; Lm <- p$Lam; R <- cr$R; m <- cr$m
    C <- chol2inv(R); mmu <- m[idx_mu]; A <- matrix(m[-idx_mu], r, k)
    sc_lam <- matrix(0, J, k); sc_psi <- sc_s2 <- numeric(J); z <- vector("list", J)
    Caa <- function(q, s) C[idx_a(q), idx_a(s), drop = FALSE]
    for (j in seq_len(J)) { B <- E[[j]]; S <- B$S; Al <- drop(A %*% Lm[j, ])
      rt <- B$qy - mmu[j] * B$q1 - Pt(B, Al); z[[j]] <- S * rt
      cma <- numeric(r); Cl <- matrix(0, r, r)
      for (q in seq_len(k)) { if (Lm[j, q] == 0) next; cma <- cma + Lm[j, q] * C[j, idx_a(q)]
        for (s in seq_len(k)) if (Lm[j, s] != 0) Cl <- Cl + Lm[j, q] * Lm[j, s] * Caa(q, s) }
      dPCP <- if (B$full) c(d * diag(Cl), rep(0, B$n - r)) else colSums(B$P * (Cl %*% B$P))
      vr <- B$q1^2 * C[j, j] + 2 * B$q1 * Pt(B, cma) + dPCP; E2 <- rt^2 + vr
      sc_psi[j] <- -0.5 * sum(S * B$ev) + 0.5 * sum(S^2 * B$ev * E2)
      sc_s2[j] <- -0.5 * sum(S) + 0.5 * sum(S^2 * E2)
      MA <- if (B$full) B$Mdiag * A else B$M %*% A
      for (q in seq_len(k)) if (lam_free[j, q]) {
        tr <- 0; for (s in seq_len(k)) if (Lm[j, s] != 0) {
          tr <- tr + Lm[j, s] * (if (B$full) sum(B$Mdiag * diag(Caa(s, q))) else sum(B$M * Caa(s, q))) }
        sc_lam[j, q] <- sum(A[, q] * (B$uy - mmu[j] * B$u1 - drop(MA %*% Lm[j, ]))) - sum(C[j, idx_a(q)] * B$u1) - tr } }
    score <- c(sc_lam[lam_free], sc_psi * p$psi, sc_s2 * p$s2)          # log-scale for psi and s2
    ## AI matrix: columns of W = dV/dtheta P y, in the Q basis of every environment
    Wt <- matrix(0, NL, J); for (j in seq_len(J)) Wt[E[[j]]$li, j] <- drop(E[[j]]$Q %*% z[[j]])
    GW <- GL %*% Wt
    QW <- lapply(seq_len(J), function(j) matrix(0, E[[j]]$n, npar))
    QGW <- lapply(seq_len(J), function(j) crossprod(E[[j]]$Q, GW[E[[j]]$li, , drop = FALSE]))
    col <- 0
    for (q in seq_len(k)) { GWl <- lapply(seq_len(J), function(bb) drop(QGW[[bb]] %*% Lm[, q]))
     for (j in seq_len(J)) if (lam_free[j, q]) { col <- col + 1
      for (bb in seq_len(J)) QW[[bb]][, col] <- QGW[[bb]][, j] * Lm[bb, q] + (if (bb == j) GWl[[bb]] else 0) } }
    for (j in seq_len(J)) QW[[j]][, nl + j] <- p$psi[j] * E[[j]]$ev * z[[j]]           # psi (log scale)
    for (j in seq_len(J)) QW[[j]][, nl + J + j] <- p$s2[j] * z[[j]]                    # s2 (log scale)
    rhs <- matrix(0, nth, npar)
    for (j in seq_len(J)) { SQW <- E[[j]]$S * QW[[j]]; rhs[j, ] <- colSums(E[[j]]$q1 * SQW)
      PS <- PxM(E[[j]], SQW); for (q in seq_len(k)) if (Lm[j, q] != 0) rhs[idx_a(q), ] <- rhs[idx_a(q), ] + Lm[j, q] * PS }
    sol <- backsolve(R, forwardsolve(t(R), rhs)); AI <- matrix(0, npar, npar)
    for (j in seq_len(J)) { B <- E[[j]]; X <- matrix(0, r, npar)
      for (q in seq_len(k)) if (Lm[j, q] != 0) X <- X + Lm[j, q] * sol[idx_a(q), , drop = FALSE]
      PtX <- if (B$full) rbind(sqrt(d) * X, matrix(0, B$n - r, npar)) else crossprod(B$P, X)
      zW <- B$S * (QW[[j]] - outer(B$q1, sol[j, ]) - PtX); AI <- AI + crossprod(QW[[j]], zW) }
    list(score = score, AI = 0.5 * (AI + t(AI)) / 1, z = z, A = A, mmu = mmu)
  }

  th <- pack(Lam, psi, s2); cr <- core(th); it <- 0; conv <- FALSE; hist <- cr$ll
  while (it < maxit) { it <- it + 1
    sa <- score_ai(cr); H <- sa$AI; g <- sa$score
    ## parameters at the lower bound with a negative score are kept fixed
    lb <- c(rep(-Inf, nl), lo_psi, lo_s2); fix <- th <= lb + 1e-8 & g < 0; g[fix] <- 0
    H2 <- H; H2[fix, ] <- 0; H2[, fix] <- 0; diag(H2)[fix] <- 1
    step <- tryCatch(solve(H2 + diag(1e-8 * max(1, max(abs(diag(H2)))), npar), g), error = function(e) g / max(1, max(abs(diag(H2)))))
    mx <- max(abs(step)); if (mx > 2) step <- step * 2 / mx                       # trust region
    a <- 1; ok <- FALSE
    for (h in 1:12) { tn <- pmax(th + a * step, lb); cn <- core(tn)
      if (is.finite(cn$ll) && cn$ll >= cr$ll - 1e-10) { ok <- TRUE; break }; a <- a / 2 }
    if (!ok) { conv <- TRUE; break }
    dll <- cn$ll - cr$ll; th <- tn; cr <- cn; hist <- c(hist, cr$ll)
    if (verbose) message(sprintf("it %2d  ll = %.4f  dll = %.2e  step = %.3f", it, cr$ll, dll, a * mx))
    if (dll < tol) { conv <- TRUE; break } }
  cr <- core(th); sa <- score_ai(cr); p <- cr$p
  ## predictions of the testing lines in the target environment
  jT <- match(target, envs); BT <- E[[jT]]
  te <- intersect(test_lines, rownames(G))
  Fte <- (G[te, L, drop = FALSE] %*% U) %*% diag(1 / sqrt(d), r)                 # F for any line (in or out of L)
  common <- drop(Fte %*% (sa$A %*% p$Lam[jT, ]))
  oT <- L[BT$li]; spec <- p$psi[jT] * drop(G[te, oT, drop = FALSE] %*% (BT$Q %*% sa$z[[jT]]))
  pred <- setNames(sa$mmu[jT] + common + spec, te); if (verbose) { message("mu ", sa$mmu[jT]); message("common ", paste(round(head(common),4), collapse=" ")); message("spec ", paste(round(head(spec),4),collapse=" ")) }
  Sg <- tcrossprod(p$Lam) + diag(p$psi, J); dimnames(Sg) <- list(envs, envs)
  list(pred = pred, ll = cr$ll, iter = it, converged = conv, max_score = max(abs(sa$score)),
       Lambda = p$Lam, psi = setNames(p$psi, envs), s2 = setNames(p$s2, envs), Sigma_g = Sg,
       rcor = cov2cor(Sg), rank = r, n_lines = NL, k = k, hist = hist)
}

## number of factors: largest k <= 2 that is identified, (J - k)^2 >= J + k
fa_k <- function(J, kmax = 2) { for (k in kmax:1) if ((J - k)^2 >= J + k) return(k); 1 }
