#' PROFIT with IPCW weights (Route A)
#'
#' Weighted variant of profit_data_driven(). Source this AFTER Profit_V1.R;
#' it defines a separate function name so existing calls are unaffected.
#'
#' Weights enter in two places:
#'   (1) the marginal covariance, via row-scaling of the demeaned data by
#'       sqrt(w) - since sum_i w_i x_i x_i' = sum_i (sqrt(w_i) x_i)(...)',
#'       this is eq. 8 of Koner et al. (2024) with v_i proportional to w_ij;
#'   (2) the pre-whitened design in the pLRT step, via row-scaling of
#'       X.hat, Z.hat, Y.hat by sqrt(w), giving the weighted pseudo-
#'       likelihood. RLRsim simulates from these matrices, so the null
#'       reference distribution matches the realized weighted design.
#'
#' The bivariate mean (fbps) is left unweighted by design: it only demeans
#' the data prior to the covariance step and therefore affects the choice
#' of eigenbasis, not the test statistic, which is computed from raw Y.
#'
#' @param wts observation-level weights, length nrow(Y), in the SAME ROW
#'   ORDER as Y. Internally rescaled to mean 1. NULL reproduces the
#'   unweighted test exactly.
#' @param weight_vc if TRUE, also weight the variance-component lme via
#'   varFixed(~ 1/w). Off by default: it slows convergence and the
#'   residual covariance is a nuisance parameter here.

profit_data_driven_wtd <- function(Y, s, tt, m, visit.time, ind.tgrid,
                                   pve = 0.85, num_knots = 2,
                                   wts = NULL, weight_vc = FALSE) {
  
  nt            <- length(tt)
  R             <- length(s)
  n             <- length(m)
  id            <- unlist(purrr::map(1:n, ~ rep(.x, m[.x])))
  not.represent <- setdiff(1:nt, unique(unlist(ind.tgrid)))
  
  ## ---- weights ---------------------------------------------------
  N <- nrow(Y)
  if (is.null(wts)) wts <- rep(1, N)
  stopifnot(length(wts) == N, all(is.finite(wts)), all(wts > 0))
  w  <- wts / mean(wts)          # mean 1 preserves the effective sample size
  sw <- sqrt(w)
  
  ## ---- Part 1: mean, marginal covariance, eigenbasis --------------
  Ymat.fit      <- fbps(Y, covariates = list(x = visit.time, z = s))
  Ymat.demeaned <- Y - Ymat.fit$Yhat
  
  # WEIGHTED marginal covariance: row i scaled by sqrt(w_i)
  Ymat.fpca <- fpca.face(Ymat.demeaned * sw, pve = pve)
  
  Phi.hat   <- (Ymat.fpca$efunctions) * sqrt(R)
  npc       <- Ymat.fpca$npc
  Ymat.proj <- Y %*% Phi.hat / R
  
  ## ---- Part 2: testing in the projected model ---------------------
  p     <- 1
  K     <- num_knots
  knots <- unname(quantile(tt, probs = seq(0, 1, length = K + 2)[-c(1, K + 2)]))
  X     <- do.call(cbind, lapply(0:p, function(i) visit.time^i))
  Z     <- outer(visit.time, knots, FUN = "-")
  Z     <- {Z * (Z > 0)}^p
  
  TestStat <- matrix(NA, 2, npc)
  pValue   <- matrix(NA, 2, npc)
  
  for (irow in 1:npc) {
    
    y <- Ymat.proj[, irow]
    
    dat.fpca <- data.frame(id = c(rep(0, length(not.represent)), id),
                           tindex = c(not.represent, unlist(ind.tgrid)),
                           value = y) %>%
      tidyr::spread(tindex, value) %>%
      dplyr::filter(id > 0) %>% dplyr::select(-id) %>% as.matrix
    
    dat.imputed <- dat.fpca
    col_means   <- colMeans(dat.imputed, na.rm = TRUE)
    for (i in 1:nrow(dat.imputed)) {
      na_idx <- which(is.na(dat.imputed[i, ]))
      if (length(na_idx) > 0) dat.imputed[i, na_idx] <- col_means[na_idx]
    }
    
    pcaObj        <- prcomp(dat.imputed, center = TRUE, scale. = TRUE)
    explained_var <- cumsum(pcaObj$sdev^2) / sum(pcaObj$sdev^2)
    npc_sub       <- which(explained_var >= pve)[1]
    
    ef     <- pcaObj$rotation
    ef.all <- do.call("rbind", lapply(ind.tgrid, function(t) ef[t, , drop = FALSE]))
    
    dat_lme        <- cbind.data.frame(rep(1, length(id)), id, y, X, Z, ef.all)
    names(dat_lme) <- c("one", "id", "Y", paste0("X", 1:ncol(X)),
                        paste0("Z", 1:ncol(Z)), paste0("EF", 1:ncol(ef.all)))
    
    random.model <- as.formula(paste0("~ -1 +", paste0("Z", 1:ncol(Z), collapse = "+")))
    ef.model     <- as.formula(paste0("~ -1 +", paste0("EF", 1:ncol(ef.all), collapse = "+")))
    fixed.model  <- as.formula(paste0("Y ~ -1 + ", paste0("X", 1:ncol(X), collapse = " + ")))
    
    if (weight_vc) {
      dat_lme$.invw <- 1 / w
      lmeObj <- tryCatch(
        lme(fixed.model,
            random  = list(one = pdIdent(random.model), id = ef.model),
            weights = nlme::varFixed(~ .invw),
            data = dat_lme, control = lmeControl(opt = "optim")),
        error = function(e) {
          warning("weighted varFixed lme failed at component ", irow,
                  "; falling back to unweighted variance components")
          lme(fixed.model,
              random = list(one = pdIdent(random.model), id = ef.model),
              data = dat_lme, control = lmeControl(opt = "optim"))
        })
    } else {
      lmeObj <- lme(fixed.model,
                    random = list(one = pdIdent(random.model), id = ef.model),
                    data = dat_lme, control = lmeControl(opt = "optim"))
    }
    
    sig <- lmeObj$sigma^2
    vc  <- lapply(rev(lmeObj$modelStruct$reStruct), VarCorr)   # renamed from `m`
    ev  <- as.vector(vc$id[, "Variance"]) * sig
    
    vcov       <- tcrossprod(ef * rep(sqrt(ev), rep.int(nrow(ef), length(ev))))
    diag(vcov) <- diag(vcov) + ifelse(sig > 0, sig, 1e-5)
    sighat     <- Matrix::bdiag(purrr::map(ind.tgrid, ~ `[`(vcov, i = .x, j = .x)))
    
    chol_sig <- t(chol(sighat))
    X.hat    <- forwardsolve(chol_sig, X)
    Z.hat    <- forwardsolve(chol_sig, Z)
    Y.hat    <- as.vector(forwardsolve(chol_sig, y))
    
    # WEIGHTED pseudo-likelihood: scale the whitened design by sqrt(w)
    X.hat <- X.hat * sw
    Z.hat <- Z.hat * sw
    Y.hat <- Y.hat * sw
    
    ## Method 1: nlme + RLRsim
    modeldat        <- cbind.data.frame(Y.hat, X.hat, Z.hat)
    names(modeldat) <- c("Y", paste0("X", 1:ncol(X.hat)), paste0("Z", 1:ncol(Z.hat)))
    modeldat$one    <- rep(1, nrow(modeldat))
    random.model    <- paste0("~ -1 +", paste0("Z", 1:ncol(Z.hat), collapse = "+"))
    
    mA_model <- rlang::call2("lme", !!!list(
      fixed  = reformulate(termlabels = paste0("X", 1:ncol(X.hat)),
                           response = "Y", intercept = FALSE),
      random = list(one = pdIdent(as.formula(random.model))),
      data   = quote(modeldat), method = "ML"))
    mA   <- eval(mA_model)
    m0   <- lm(Y ~ X1 - 1, data = modeldat)
    eLRT <- RLRsim::exactLRT(m = mA, m0 = m0, nsim = 1e5)
    
    ## Method 2: direct optimization
    XhtXh     <- crossprod(X.hat)
    ZhtZh     <- crossprod(Z.hat)
    ZhtXh     <- crossprod(Z.hat, X.hat)
    XhtYh     <- as.vector(crossprod(X.hat, Y.hat))
    ZhtYh     <- as.vector(crossprod(Z.hat, Y.hat))
    YhtYh     <- sum(Y.hat^2)
    eig.ZhtZh <- svd(ZhtZh)$d
    
    neg_profile_log_lik <- function(loglam, XtX, ZtX, XtY, ZtZ, ZtY, YtY, eig.ZtZ) {
      lam          <- exp(loglam)
      H            <- lam * ZtZ
      diag(H)      <- diag(H) + 1
      chol.H       <- chol(H)
      choltinvZtX  <- forwardsolve(t(chol.H), ZtX)
      XtHinvX      <- XtX - lam * crossprod(choltinvZtX)
      choltinvZtY  <- as.vector(forwardsolve(t(chol.H), ZtY))
      XtHinvY      <- XtY - lam * as.vector(crossprod(choltinvZtX, choltinvZtY))
      chol.XtHinvX <- chol(XtHinvX)
      p2.lik       <- as.vector(forwardsolve(t(chol.XtHinvX), XtHinvY))
      p1.lik2      <- YtY - lam * sum(choltinvZtY^2)
      p1.lik2 - sum(p2.lik^2) + sum(log(1 + lam * eig.ZtZ))
    }
    
    neglikfit <- optimize(f = neg_profile_log_lik, interval = c(-40, 40),
                          XtX = XhtXh, ZtX = ZhtXh, XtY = XhtYh, ZtZ = ZhtZh,
                          ZtY = ZhtYh, YtY = YhtYh, eig.ZtZ = eig.ZhtZh)
    
    negLik.H0uHA <- neglikfit$objective
    negLik.H0    <- YhtYh - (1 / XhtXh[1, 1]) * ((XhtYh[1])^2)
    pLRT         <- negLik.H0 - negLik.H0uHA
    sim.LRT      <- RLRsim::LRTSim(X = X.hat, Z = Z.hat, q = ncol(X.hat) - 1,
                                   sqrt.Sigma = diag(ncol(Z.hat)))
    pval.LRT     <- mean(sim.LRT > pLRT)
    
    pValue[,  irow] <- c(eLRT$p.value,   pval.LRT)
    TestStat[, irow] <- c(eLRT$statistic, pLRT)
  }
  
  list(max_TestStat    = apply(round(TestStat, 3), 1, max),
       pValue          = round(pValue, 6),
       min_pValue      = apply(round(pValue, 3), 1, min),
       K               = npc,
       K_with_min_pVal = apply(round(pValue, 3), 1, which.min),
       phi             = Phi.hat,
       pc_score        = Ymat.fpca$scores,
       eigen_val       = Ymat.fpca$evalues,
       weights_used    = !is.null(wts) && var(w) > 0)
}


# ============================================================
# Usage
# ============================================================
# source("Profit_V1.R")
# source("Profit_V1_weighted.R")
#
# ## rebuild inputs in a guaranteed-consistent row order
# dw          <- dw %>% arrange(spid, round)          # from nhats_ipcw.R
# func_cols   <- grep("^dmeanmin", names(dw), value = TRUE)
# Yw          <- as.matrix(dw[, func_cols])
# ids         <- unique(dw$spid)
# mw          <- as.integer(table(factor(dw$spid, levels = ids)))
# vtl         <- split(dw$round, factor(dw$spid, levels = ids))
# visit.time  <- unlist(vtl, use.names = FALSE)
# tt          <- sort(unique(dw$round))
# ind.tgrid   <- lapply(vtl, function(x) match(x, tt))
# wts         <- dw$w                                  # SAME row order as Yw
#
# stopifnot(length(wts) == nrow(Yw),
#           identical(as.character(dw$spid), as.character(rep(ids, mw))))
#
# ## unweighted reference (should reproduce your published numbers)
# fit_unw <- profit_data_driven_wtd(Y = Yw, s = seq_along(func_cols), tt = tt,
#                                   m = mw, visit.time = visit.time,
#                                   ind.tgrid = ind.tgrid,
#                                   pve = 0.9, num_knots = 3, wts = NULL)
#
# ## IPCW-weighted
# fit_ipcw <- profit_data_driven_wtd(Y = Yw, s = seq_along(func_cols), tt = tt,
#                                    m = mw, visit.time = visit.time,
#                                    ind.tgrid = ind.tgrid,
#                                    pve = 0.9, num_knots = 3, wts = wts)
#
# summarise_fit <- function(f) {
#   p <- as.numeric(f$pValue[1, ])
#   data.frame(K = f$K, min_p = min(p), cutoff = 0.05 / f$K,
#              reject = min(p) < 0.05 / f$K)
# }
# rbind(unweighted = summarise_fit(fit_unw),
#       ipcw       = summarise_fit(fit_ipcw))