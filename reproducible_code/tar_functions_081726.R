########################################################################
# tar_functions_081726.R
#
# Core estimation and data-generation functions for the TAR simulation
# study. Defines four PTE-estimation approaches (full surrogate set,
# robust regularization, TAR weights-only, TAR), the eight simulation
# settings, population-truth calculators, and the plotting/table-writing
# functions used to build the manuscript's figures and tables.
#
# Sourced by tar_run_setting_081726.R (one setting at a time) and by
# tar_combine_results_081726.R (plots, tables, and the large-sample R_Q
# diagnostic, after all eight settings have been run).
########################################################################

library(glmnet)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(Rsurrogate)
library(MASS)
library(future.apply)

###################################################################
## Post-lasso debiasing helpers
###################################################################

debias_ols <- function(original_coefs, data_train) {
  nonzero_coefs <- original_coefs[original_coefs != 0]
  selected_cols <- c(names(nonzero_coefs))
  X_selected <- as.matrix(data_train[, selected_cols])
  X_selected <- cbind(Intercept = 1, X_selected)
  # If the number of selected columns is not less than the number of rows,
  # the refit is not identified; return an all-NA coefficient vector of the
  # right shape rather than let qr.solve() silently return partial NAs.
  if (ncol(X_selected) >= nrow(X_selected)) {
    debias_coef <- matrix(NA_real_, nrow = ncol(X_selected), ncol = 1)
    rownames(debias_coef) <- c("Intercept", selected_cols)
    return(debias_coef)
  }
  # qr.solve() gives the OLS solution via QR decomposition, which is more
  # numerically stable than solve(t(X) %*% X) when selected surrogates are
  # correlated.
  debias_coef <- as.matrix(qr.solve(X_selected, data_train$y))
  rownames(debias_coef)[-1] = selected_cols
  return(debias_coef)
}

debias_pred = function(debias_coef, data_val) {
  preds = as.matrix(data_val[ , match(rownames(debias_coef), names(data_val)) ])%*%as.matrix(debias_coef)
  return(preds)
}

debias_p = function(debias_coef, penalty) {
  noint = debias_coef[-1]
  return(sum(noint != 0 ))
}

debias_deltas = function(debias_coef, data_val, control_data) {
  debias_coef_subset = debias_coef[-1, , drop = FALSE]
  debias_coef_subset_vec = as.vector(debias_coef_subset)
  names(debias_coef_subset_vec) = rownames(debias_coef_subset)
  # If this fold x lambda's refit was unidentified (see debias_ols()), do
  # not pass NA coefficients into make.delta.s(); return NA so the caller
  # can average with na.rm = TRUE.
  if (any(is.na(debias_coef_subset_vec))) return(NA_real_)
  return(make.delta.s(debias_coef_subset_vec, control = control_data, treat = data_val))
}

make.delta.s <- function(coefs_vec, treat, control) {
  if (length(coefs_vec) == 0) {
    return(mean(treat$y) - mean(control$y))
  } else {
    selected_surrogates <- names(coefs_vec)
    control.surrogates <- as.matrix(control[, selected_surrogates, drop = FALSE])
    treat.surrogates <- as.matrix(treat[, selected_surrogates, drop = FALSE])
    y.treat <- treat$y

    s.tilde.0 <- control.surrogates %*% coefs_vec
    s.tilde.1 <- treat.surrogates %*% coefs_vec

    mu.1.s0 <- vapply(s.tilde.0, pred.smooth, FUN.VALUE = numeric(1), zz = s.tilde.1, y1 = y.treat)
    if (sum(is.na(mu.1.s0)) > 0) {
      # Nearest-neighbor fallback for any point the kernel smoother could
      # not evaluate.
      c.mat = cbind(s.tilde.0, mu.1.s0)
      for (o in 1:length(mu.1.s0)) {
        if (is.na(mu.1.s0[o])) {
          distance = abs(s.tilde.0 - s.tilde.0[o])
          c.temp = cbind(c.mat, distance)
          c.temp = c.temp[!is.na(c.temp[, 2]), ]
          new.est = c.temp[c.temp[, 3] == min(c.temp[, 3]), 2]
          mu.1.s0[o] = new.est[1]
        }
      }
    }
    return(mean(mu.1.s0) - mean(control$y))
  }
}

te.calculate = function(s, treat, control) {
  mean_diff <- mean(treat[[s]]) - mean(control[[s]])
  pooled_sd <- sqrt((var(treat[[s]]) + var(control[[s]])) / 2)
  return(abs(mean_diff / pooled_sd))
}

###################################################################
## .bootstrap_var(): resamples both arms with replacement and reruns
## fit_once() (the full estimation algorithm, including re-selection of
## lambda*) boot.num times, returning bootstrap SE + 95% percentile CI
## for delta, delta.s, and R.s (the PTE). Resamples run in parallel via
## future_lapply() using whatever backend future::plan() is set to in the
## calling script.
###################################################################

.bootstrap_var <- function(data, fit_once, boot.num) {
  control <- data[data$treat == 0, ]
  treat   <- data[data$treat == 1, ]

  boot_results <- future_lapply(seq_len(boot.num), function(b) {
    control_b <- control[sample(nrow(control), nrow(control), replace = TRUE), ]
    treat_b   <- treat[sample(nrow(treat), nrow(treat), replace = TRUE), ]
    dat_b <- rbind(treat_b, control_b)
    tryCatch(fit_once(dat_b),
             error = function(e) list(delta = NA_real_, delta.s = NA_real_, pte = NA_real_))
  }, future.seed = TRUE)

  boot_delta   <- vapply(boot_results, function(r) r$delta,   numeric(1))
  boot_delta.s <- vapply(boot_results, function(r) r$delta.s, numeric(1))
  boot_R.s     <- vapply(boot_results, function(r) r$pte,     numeric(1))

  list(
    delta_se   = sd(boot_delta,   na.rm = TRUE),
    delta_ci   = unname(quantile(boot_delta,   c(.025, .975), na.rm = TRUE)),
    delta.s_se = sd(boot_delta.s, na.rm = TRUE),
    delta.s_ci = unname(quantile(boot_delta.s, c(.025, .975), na.rm = TRUE)),
    R.s_se     = sd(boot_R.s,     na.rm = TRUE),
    R.s_ci     = unname(quantile(boot_R.s,     c(.025, .975), na.rm = TRUE)),
    boot.num   = boot.num
  )
}

###################################################################
## full.surrogate.approach(): full-candidate-set PTE via
## Rsurrogate::R.s.estimate(), no selection.
###################################################################

full.surrogate.approach <- function(data, surrogate_names, var = FALSE, boot.num = 200) {
  fit_once <- function(dat) {
    r <- tryCatch({
      R.s.estimate(yone = dat$y[dat$treat == 1], yzero = dat$y[dat$treat == 0],
                   sone = as.matrix(dat[dat$treat == 1, surrogate_names]),
                   szero = as.matrix(dat[dat$treat == 0, surrogate_names]),
                   number = "multiple", extrapolate = TRUE)
    }, error = function(e) e)
    if (inherits(r, "error")) {
      return(list(delta = NA_real_, delta.s = NA_real_, pte = NA_real_))
    }
    list(delta = r$delta, delta.s = r$delta.s, pte = r$R.s)
  }

  point <- fit_once(data)
  point$var <- if (var) .bootstrap_var(data, fit_once, boot.num) else NULL
  point
}

###################################################################
## smoothed.regularized(): robust regularization -- ordinary (unweighted)
## CV lasso, selected via lambda.1se, followed by a post-lasso OLS refit.
###################################################################

smoothed.regularized = function(data, surrogate_names, nlambda = 20, lambda_rule = "lambda.1se",
                                 var = FALSE, boot.num = 200) {
  fit_once <- function(dat) {
    y.treat <- dat$y[dat$treat == 1]
    X.treat <- dat[dat$treat == 1, setdiff(names(dat), "y")]
    X.treat <- as.matrix(X.treat[, setdiff(names(X.treat), "treat")])

    mean1 <- mean(dat$y[dat$treat == 1])
    mean0 <- mean(dat$y[dat$treat == 0])
    delta <- mean1 - mean0

    cvfit <- cv.glmnet(X.treat, y.treat, alpha = 1, nlambda = nlambda)
    coefs <- coef(cvfit, s = lambda_rule)
    coefs_vec <- as.vector(coefs)
    names(coefs_vec) <- rownames(coefs)

    nonzero_coefs <- coefs_vec[coefs_vec != 0]
    nonzero_coefs <- nonzero_coefs[names(nonzero_coefs) != "(Intercept)"]

    selected_cols <- c(names(nonzero_coefs))
    data.treat <- dat[dat$treat == 1, ]
    X_selected <- as.matrix(data.treat[, selected_cols])
    X_selected <- cbind(Intercept = 1, X_selected)
    # See debias_ols() for the rank-deficiency guard. An empty selection
    # falls through fine: intercept-only refit, and make.delta.s()'s
    # length-0 branch gives delta.s = delta (PTE = 0).
    if (ncol(X_selected) >= nrow(X_selected)) {
      return(list(selection = nonzero_coefs, delta = delta, delta.s = NA_real_, pte = NA_real_))
    }
    debias_coef <- as.matrix(qr.solve(X_selected, data.treat$y))
    debias_coef <- debias_coef[-1]
    names(debias_coef) <- selected_cols
    delta.s <- make.delta.s(coefs_vec = debias_coef, treat = dat[dat$treat == 1, ], control = dat[dat$treat == 0, ])

    list(selection = nonzero_coefs, delta = delta, delta.s = delta.s, pte = 1 - delta.s / delta)
  }

  point <- fit_once(data)
  point$var <- if (var) .bootstrap_var(data, fit_once, boot.num) else NULL
  point
}

###################################################################
## tar.approach(): treatment-aware regularization -- weighted lasso with
## multi-objective selection of lambda* (CV prediction error, subject to
## a treatment-effect-preservation constraint within tolerance epsilon).
###################################################################

tar.approach = function(data, folds = 5, epsilon = 0.05, nlambda = 20, surrogate_names, c.constant = 1, small.n = 1e-6,
                         var = FALSE, boot.num = 200) {

  fit_once <- function(dat) {
    control <- dat[dat$treat == 0, ]
    treat <- dat[dat$treat == 1, ]

    treatment_effect_magnitude <- sapply(surrogate_names, te.calculate, treat = treat, control = control)
    penalty_factors <- (1 / (treatment_effect_magnitude + small.n))^c.constant
    penalty_factors <- penalty_factors / mean(penalty_factors)
    names(penalty_factors) <- surrogate_names
    penalty_factors_fulldata <- penalty_factors

    y.treat <- treat$y
    X.treat <- treat[setdiff(names(treat), c("y", "treat"))]
    cvfit_fulldata <- glmnet(as.matrix(X.treat), y.treat, alpha = 1, penalty.factor = penalty_factors_fulldata[colnames(X.treat)], nlambda = nlambda)
    lambda.grid <- cvfit_fulldata$lambda

    optim.matrix <- matrix(nrow = length(lambda.grid), ncol = 4)
    optim.matrix[, 1] <- lambda.grid
    mse.mat <- matrix(nrow = length(lambda.grid), ncol = folds)
    sparse.mat <- matrix(nrow = length(lambda.grid), ncol = folds)
    residuals.mat <- matrix(nrow = length(lambda.grid), ncol = folds)

    # k-fold CV within the treated arm: every treated subject assigned to
    # exactly one of `folds` disjoint groups.
    n1_total <- sum(dat$treat == 1)
    fold.id <- sample(rep(1:folds, length.out = n1_total))

    for (uu in 1:folds) {
      ind.cv <- as.integer(fold.id != uu)   # 1 = train, 0 = validate
      sub.treat <- treat[ind.cv == 1, ]
      treatment_effect_magnitude <- sapply(surrogate_names, te.calculate, treat = sub.treat, control = control)

      penalty_factors <- (1 / (treatment_effect_magnitude + small.n))^c.constant
      penalty_factors <- penalty_factors / mean(penalty_factors)
      names(penalty_factors) <- surrogate_names

      y.treat <- sub.treat$y
      X.treat <- sub.treat[setdiff(names(sub.treat), c("y", "treat"))]

      cvfit <- glmnet(as.matrix(X.treat), y.treat, alpha = 1, penalty.factor = penalty_factors[colnames(X.treat)], lambda = lambda.grid)

      sub.treat.val <- treat[ind.cv == 0, ]
      y.treat.val <- sub.treat.val$y
      X.treat.val <- sub.treat.val[setdiff(names(sub.treat.val), c("y", "treat"))]

      coefs_all <- coef(cvfit, cvfit$lambda)[-1, ]
      result_list <- lapply(seq_len(ncol(coefs_all)), function(j) {
        column <- coefs_all[, j]
        debias_ols(column, data_train = sub.treat)
      })

      debias_predictions <- sapply(result_list, debias_pred, data_val = cbind(Intercept = 1, X.treat.val))
      mse.mat[, uu] <- length(y.treat.val)^(-1) * apply(as.matrix((debias_predictions - y.treat.val)^2), 2, sum)
      sparse.mat[, uu] <- sapply(result_list, debias_p)
      residuals.mat[, uu] <- vapply(result_list, debias_deltas, FUN.VALUE = numeric(1), data_val = sub.treat.val, control_data = control)
    }

    # Average across folds with na.rm = TRUE, so one fold's unrecoverable
    # kernel-smoothing evaluation does not poison the whole estimate for
    # that lambda.
    optim.matrix[, 2] <- apply(mse.mat, 1, mean, na.rm = TRUE)
    optim.matrix[, 3] <- apply(sparse.mat, 1, mean, na.rm = TRUE)
    optim.matrix[, 4] <- apply(residuals.mat, 1, mean, na.rm = TRUE)

    n_na_cells <- sum(is.na(residuals.mat))
    n_total_cells <- length(residuals.mat)
    if (n_na_cells > 0) {
      pct_na <- 100 * n_na_cells / n_total_cells
      if (pct_na > 5) {
        warning(sprintf("tar.approach: %.1f%% of (fold x lambda) residual evaluations were NA (%d of %d cells).",
                         pct_na, n_na_cells, n_total_cells))
      }
    }

    # Feasibility threshold: delta.s + epsilon*delta from the full
    # (unpenalized, all-candidate) working model. Falls back to the
    # least-penalized identifiable point on TAR's own lasso path if the
    # full-candidate R.s.estimate() call fails (e.g. k close to or
    # exceeding n1).
    rr_attempt <- tryCatch({
      R.s.estimate(yone = dat$y[dat$treat == 1], yzero = dat$y[dat$treat == 0],
                   sone = as.matrix(dat[dat$treat == 1, names(dat) %in% surrogate_names]),
                   szero = as.matrix(dat[dat$treat == 0, names(dat) %in% surrogate_names]),
                   number = "multiple", extrapolate = TRUE)
    }, error = function(e) e)

    if (!inherits(rr_attempt, "error")) {
      threshold <- rr_attempt$delta.s + epsilon * rr_attempt$delta
    } else {
      warning(sprintf("tar.approach: full-candidate-set R.s.estimate() failed (%s) -- falling back to the least-penalized identifiable point on TAR's own lasso path to define the threshold.",
                       conditionMessage(rr_attempt)))
      delta_full <- mean(treat$y) - mean(control$y)
      coefs_full_path <- coef(cvfit_fulldata, cvfit_fulldata$lambda)[-1, , drop = FALSE]
      delta_s_proxy <- NA_real_
      for (jcol in rev(seq_len(ncol(coefs_full_path)))) {
        dc <- debias_ols(coefs_full_path[, jcol], data_train = treat)
        dc_subset <- dc[-1, , drop = FALSE]
        dc_vec <- as.vector(dc_subset)
        names(dc_vec) <- rownames(dc_subset)
        if (!any(is.na(dc_vec))) {
          delta_s_proxy <- make.delta.s(dc_vec, treat = treat, control = control)
          break
        }
      }
      if (is.na(delta_s_proxy)) delta_s_proxy <- delta_full
      threshold <- delta_s_proxy + epsilon * delta_full
    }

    feasible <- !is.na(optim.matrix[, 4]) & (optim.matrix[, 4] <= threshold)
    n_feasible <- sum(feasible)

    if (n_feasible == 0) {
      optimal.lambda <- 0
    }
    if (n_feasible == 1) {
      optim.temp <- optim.matrix[feasible, , drop = FALSE]
      optimal.lambda <- optim.temp[1, 1]
    }
    if (n_feasible > 1) {
      optim.matrix.f <- optim.matrix[feasible, , drop = FALSE]
      optimal.lambda <- optim.matrix.f[which.min(optim.matrix.f[, 3]), 1]
    }

    y.treat <- treat$y
    X.treat <- treat[setdiff(names(treat), c("y", "treat"))]
    cvfit <- glmnet(as.matrix(X.treat), y.treat, alpha = 1, penalty.factor = penalty_factors_fulldata[colnames(X.treat)], lambda = optimal.lambda)
    coefs <- coef(cvfit, optimal.lambda)[-1, ]
    coefs_vec <- as.vector(coefs)
    names(coefs_vec) <- names(coefs)
    nonzero_coefs <- coefs_vec[coefs_vec != 0]
    selected_cols <- c(names(nonzero_coefs))
    X_selected <- as.matrix(treat[, selected_cols])
    X_selected <- cbind(Intercept = 1, X_selected)
    delta <- mean(treat$y) - mean(control$y)
    if (ncol(X_selected) >= nrow(X_selected)) {
      return(list(selection = nonzero_coefs, delta = delta, delta.s = NA_real_, pte = NA_real_))
    }
    debias_coef <- as.matrix(qr.solve(X_selected, treat$y))
    rownames(debias_coef)[-1] <- selected_cols
    debias_coef_subset <- debias_coef[-1, , drop = FALSE]
    debias_coef_subset_vec <- as.vector(debias_coef_subset)
    names(debias_coef_subset_vec) <- rownames(debias_coef_subset)
    delta.s <- make.delta.s(debias_coef_subset_vec, control = control, treat = treat)
    list(selection = nonzero_coefs, delta = delta, delta.s = delta.s, pte = 1 - delta.s / delta)
  }

  point <- fit_once(data)
  point$var <- if (var) .bootstrap_var(data, fit_once, boot.num) else NULL
  point
}

###################################################################
## tar.weightsonly(): TAR penalty weights within a single ordinary-CV
## lasso fit (no multi-objective lambda* selection).
###################################################################

tar.weightsonly = function(data, nlambda = 20, surrogate_names, c.constant = 1, small.n = 1e-6,
                            lambda_rule = "lambda.1se", var = FALSE, boot.num = 200) {

  fit_once <- function(dat) {
    control <- dat[dat$treat == 0, ]
    treat <- dat[dat$treat == 1, ]

    treatment_effect_magnitude <- sapply(surrogate_names, te.calculate, treat = treat, control = control)
    penalty_factors <- (1 / (treatment_effect_magnitude + small.n))^c.constant
    penalty_factors <- penalty_factors / mean(penalty_factors)
    names(penalty_factors) <- surrogate_names
    penalty_factors_fulldata <- penalty_factors

    y.treat <- dat$y[dat$treat == 1]
    X.treat <- dat[dat$treat == 1, setdiff(names(dat), "y")]
    X.treat <- as.matrix(X.treat[, setdiff(names(X.treat), "treat")])

    mean1 <- mean(dat$y[dat$treat == 1])
    mean0 <- mean(dat$y[dat$treat == 0])
    delta <- mean1 - mean0

    cvfit <- cv.glmnet(X.treat, y.treat, alpha = 1, penalty.factor = penalty_factors_fulldata[colnames(X.treat)], nlambda = nlambda)
    coefs <- coef(cvfit, s = lambda_rule)
    coefs_vec <- as.vector(coefs)
    names(coefs_vec) <- rownames(coefs)

    nonzero_coefs <- coefs_vec[coefs_vec != 0]
    nonzero_coefs <- nonzero_coefs[names(nonzero_coefs) != "(Intercept)"]

    selected_cols <- c(names(nonzero_coefs))
    data.treat <- dat[dat$treat == 1, ]
    X_selected <- as.matrix(data.treat[, selected_cols])
    X_selected <- cbind(Intercept = 1, X_selected)
    if (ncol(X_selected) >= nrow(X_selected)) {
      return(list(selection = nonzero_coefs, delta = delta, delta.s = NA_real_, pte = NA_real_))
    }
    debias_coef <- as.matrix(qr.solve(X_selected, data.treat$y))
    debias_coef <- debias_coef[-1]
    names(debias_coef) <- selected_cols
    delta.s <- make.delta.s(coefs_vec = debias_coef, treat = dat[dat$treat == 1, ], control = dat[dat$treat == 0, ])

    list(selection = nonzero_coefs, delta = delta, delta.s = delta.s, pte = 1 - delta.s / delta)
  }

  point <- fit_once(data)
  point$var <- if (var) .bootstrap_var(data, fit_once, boot.num) else NULL
  point
}

###################################################################
## Simulation settings 1-4 (linear, Normal candidates)
###################################################################

.gen_lin_base <- function(n1, n0, rho, mu1, mu0, sd2_1, sd2_0) {
  R <- matrix(rho, 9, 9); diag(R) <- 1
  sd_t <- sqrt(sd2_1); sd_c <- sqrt(sd2_0)
  Sigma_t <- (sd_t %*% t(sd_t)) * R
  Sigma_c <- (sd_c %*% t(sd_c)) * R
  s.treat   <- mvrnorm(n1, mu1, Sigma_t)
  s.control <- mvrnorm(n0, mu0, Sigma_c)
  list(s.treat = s.treat, s.control = s.control)
}

.assemble <- function(y1, y0, s.treat, s.control, k) {
  dat <- as.data.frame(rbind(cbind(y1, s.treat), cbind(y0, s.control)))
  names(dat) <- c("y", paste0("s", 1:k))
  dat$treat <- c(rep(1, length(y1)), rep(0, length(y0)))
  dat
}

###################################################################
## Simulation settings 5-7 (nonlinear, heterogeneous candidates)
###################################################################

.gen_mixed_corr <- function(n, rho, arm) {
  R <- matrix(rho, 9, 9); diag(R) <- 1
  Z <- mvrnorm(n, mu = rep(0, 9), Sigma = R)
  U <- pnorm(Z)
  if (arm == "treat") {
    s1 <- qnorm(U[,1], mean = 6, sd = 2)
    s2 <- qbinom(U[,2], size = 1, prob = 0.6)
    s3 <- qexp(U[,3], rate = 1/2)
  } else {
    s1 <- qnorm(U[,1], mean = 5, sd = 1)
    s2 <- qbinom(U[,2], size = 1, prob = 0.1)
    s3 <- qexp(U[,3], rate = 1/0.5)
  }
  s4 <- qnorm(U[,4],  mean = 0,  sd = 2)
  s5 <- qnorm(U[,5],  mean = 1,  sd = 2)
  s6 <- qgamma(U[,6], shape = 2, scale = 2)
  s7 <- qnorm(U[,7],  mean = 0,  sd = 2)
  s8 <- qnorm(U[,8],  mean = -1, sd = 0.5)
  s9 <- qgamma(U[,9], shape = 2, scale = 2)
  cbind(s1, s2, s3, s4, s5, s6, s7, s8, s9)
}

y1.nl <- function(s) 5   + 0.7*s[,1]^2*s[,2] + 0.5*s[,1]*s[,2]*s[,3] + exp(0.2*s[,1]*s[,2]) + exp(rnorm(nrow(s), 0, 0.09))
y0.nl <- function(s) 0.5 + 0.5*s[,1]^2*s[,2] + 0.3*s[,1]*s[,2]*s[,3] + exp(0.1*s[,1]*s[,2]) + exp(rnorm(nrow(s), 0, 0.09))

y1.nl.rogue <- function(s, coef = c(1, 0.5, 0.5))
  y1.nl(s) + coef[1]*s[,4] + coef[2]*s[,5] + coef[3]*s[,6]
y0.nl.rogue <- function(s, coef = c(1, 0.5, 0.5))
  y0.nl(s) + coef[1]*s[,4] + coef[2]*s[,5] + coef[3]*s[,6]

.gen_redundant_nl <- function(n, arm) {
  # S1, S2 are identically distributed in both arms, so they carry no
  # treatment-effect information on their own even though they have
  # nonzero coefficients in the outcome model below. Only S3 differs by
  # arm. True minimal set S0 = {s3}.
  s1 <- rnorm(n, 6, 1.5)
  s2 <- rbinom(n, size = 1, prob = 0.5)
  if (arm == "treat") {
    s3 <- rexp(n, rate = 1/3)
  } else {
    s3 <- rexp(n, rate = 1/0.3)
  }
  s4 <- rnorm(n, 0, 2); s5 <- rnorm(n, 1, 2); s6 <- rgamma(n, shape = 2, scale = 2)
  s7 <- rnorm(n, 0, 2); s8 <- rnorm(n, -1, 0.5); s9 <- rgamma(n, shape = 2, scale = 2)
  cbind(s1, s2, s3, s4, s5, s6, s7, s8, s9)
}

y1.redundant.nl <- function(s) 7.0 + 0.15*s[,1]^2*s[,2] + 2.2*s[,1]*s[,2]*s[,3] + exp(0.05*s[,1]*s[,2]) + rnorm(nrow(s), 0, 1.5)
y0.redundant.nl <- function(s) 2.0 + 0.15*s[,1]^2*s[,2] + 2.2*s[,1]*s[,2]*s[,3] + exp(0.05*s[,1]*s[,2]) + rnorm(nrow(s), 0, 1.5)

###################################################################
## Setting 8 (many-candidate pool)
###################################################################

.gen_setting8 <- function(n1, n0, k_noise = 77, block = 10, rho_block = 0.2) {
  R3 <- matrix(0.1, 3, 3); diag(R3) <- 1
  sd_t3 <- sqrt(c(4, 4, 8)); sd_c3 <- sqrt(c(1, 1, 2))
  Sigma_t3 <- (sd_t3 %*% t(sd_t3)) * R3
  Sigma_c3 <- (sd_c3 %*% t(sd_c3)) * R3
  s_true_t <- mvrnorm(n1, c(6, 6, 4), Sigma_t3)
  s_true_c <- mvrnorm(n0, c(5, 5, 2), Sigma_c3)

  n_blocks <- ceiling(k_noise / block)
  Rb <- matrix(rho_block, block, block); diag(Rb) <- 1
  gen_noise <- function(n) {
    cols <- vector("list", n_blocks)
    for (bIdx in 1:n_blocks) {
      cols[[bIdx]] <- mvrnorm(n, rep(0, block), Rb)
    }
    do.call(cbind, cols)[, 1:k_noise, drop = FALSE]
  }
  noise_t <- gen_noise(n1)
  noise_c <- gen_noise(n0)

  y1 <- 5.35 + 2*s_true_t[,1] + 4*s_true_t[,2] + 2.5*s_true_t[,3] + rnorm(n1, 0, 3)
  y0 <- 2    + 1.92*s_true_c[,1] + 3.5*s_true_c[,2] + 2.3*s_true_c[,3] + rnorm(n0, 0, 3)

  dat <- as.data.frame(rbind(
    cbind(y1, s_true_t, noise_t),
    cbind(y0, s_true_c, noise_c)
  ))
  names(dat) <- c("y", paste0("s", 1:(3 + k_noise)))
  dat$treat <- c(rep(1, n1), rep(0, n0))
  dat
}

###################################################################
## Master dispatcher: draws one dataset (n1 treated, n0 control) for the
## given setting (1-8).
###################################################################

simulate_setting <- function(n1, n0, setting) {
  if (setting == 1) {
    g <- .gen_lin_base(n1, n0, rho = 0,
                        mu1 = c(6,6,4,1,1,4,1,1,4), mu0 = c(5,5,2,1,1,4,1,1,4),
                        sd2_1 = c(4,4,8,4,4,8,4,4,8), sd2_0 = c(1,1,2,4,4,8,4,4,8))
    y1 <- 5.35 + 2*g$s.treat[,1] + 4*g$s.treat[,2] + 2.5*g$s.treat[,3] + rnorm(n1, 0, 3)
    y0 <- 2    + 1.92*g$s.control[,1] + 3.5*g$s.control[,2] + 2.3*g$s.control[,3] + rnorm(n0, 0, 3)
    return(.assemble(y1, y0, g$s.treat, g$s.control, 9))
  }

  if (setting == 2) {
    g <- .gen_lin_base(n1, n0, rho = 0.5,
                        mu1 = c(6,6,4,1,1,4,1,1,4), mu0 = c(5,5,2,1,1,4,1,1,4),
                        sd2_1 = c(4,4,8,4,4,8,4,4,8), sd2_0 = c(1,1,2,4,4,8,4,4,8))
    y1 <- 5.35 + 2*g$s.treat[,1] + 4*g$s.treat[,2] + 2.5*g$s.treat[,3] + rnorm(n1, 0, 3)
    y0 <- 2    + 1.92*g$s.control[,1] + 3.5*g$s.control[,2] + 2.3*g$s.control[,3] + rnorm(n0, 0, 3)
    return(.assemble(y1, y0, g$s.treat, g$s.control, 9))
  }

  if (setting == 3) {
    # Rogue variables (identical distribution both arms, so no treatment
    # effect, despite predicting Y). True minimal set S0 = {s1,s2,s3}.
    g <- .gen_lin_base(n1, n0, rho = 0.1,
                        mu1 = c(6,6,4,1,1,4,1,1,4), mu0 = c(5,5,2,1,1,4,1,1,4),
                        sd2_1 = c(4,4,8,4,4,8,4,4,8), sd2_0 = c(1,1,2,4,4,8,4,4,8))
    y1 <- 5.35 + 2*g$s.treat[,1]   + 4*g$s.treat[,2]   + 2.5*g$s.treat[,3]   +
                 0.6*g$s.treat[,4]   + 0.3*g$s.treat[,5] + 0.3*g$s.treat[,6]   + rnorm(n1, 0, 3)
    y0 <- 2    + 1.92*g$s.control[,1] + 3.5*g$s.control[,2] + 2.3*g$s.control[,3] +
                 0.6*g$s.control[,4]   + 0.3*g$s.control[,5] + 0.3*g$s.control[,6]   + rnorm(n0, 0, 3)
    return(.assemble(y1, y0, g$s.treat, g$s.control, 9))
  }

  if (setting == 4) {
    # Redundant surrogates / minimal set. True minimal set S0 = {s3}.
    d <- 0.13
    mu1 <- c(6 + d, 6 + d, 4, 1, 1, 4, 1, 1, 4)
    mu0 <- c(6,     6,     2, 1, 1, 4, 1, 1, 4)
    g <- .gen_lin_base(n1, n0, rho = 0.08, mu1 = mu1, mu0 = mu0,
                        sd2_1 = c(4,4,8,4,4,8,4,4,8), sd2_0 = c(4,4,2,4,4,8,4,4,8))
    y1 <- 2 + 4*g$s.treat[,1]   + 4*g$s.treat[,2]   + 2.5*g$s.treat[,3]   + rnorm(n1, 0, 3)
    y0 <- 2 + 4*g$s.control[,1] + 4*g$s.control[,2] + 1.6*g$s.control[,3] + rnorm(n0, 0, 3)
    return(.assemble(y1, y0, g$s.treat, g$s.control, 9))
  }

  if (setting == 5) {
    s.treat   <- .gen_mixed_corr(n1, rho = 0.3, arm = "treat")
    s.control <- .gen_mixed_corr(n0, rho = 0.3, arm = "control")
    y1 <- y1.nl(s.treat)
    y0 <- y0.nl(s.control)
    return(.assemble(y1, y0, s.treat, s.control, 9))
  }

  if (setting == 6) {
    s.treat   <- .gen_mixed_corr(n1, rho = 0.1, arm = "treat")
    s.control <- .gen_mixed_corr(n0, rho = 0.1, arm = "control")
    y1 <- y1.nl.rogue(s.treat)
    y0 <- y0.nl.rogue(s.control)
    return(.assemble(y1, y0, s.treat, s.control, 9))
  }

  if (setting == 7) {
    # Nonlinear/heterogeneous redundant surrogates / minimal set. True
    # minimal set S0 = {s3}.
    s.treat   <- .gen_redundant_nl(n1, arm = "treat")
    s.control <- .gen_redundant_nl(n0, arm = "control")
    y1 <- y1.redundant.nl(s.treat)
    y0 <- y0.redundant.nl(s.control)
    return(.assemble(y1, y0, s.treat, s.control, 9))
  }

  if (setting == 8) {
    return(.gen_setting8(n1, n0, k_noise = 77, block = 10, rho_block = 0.2))
  }

  stop("Unknown setting: ", setting)
}

true_minimal_set <- list(
  `1` = c("s1","s2","s3"),
  `2` = c("s1","s2","s3"),
  `3` = c("s1","s2","s3"),
  `4` = c("s3"),
  `5` = c("s1","s2","s3"),
  `6` = c("s1","s2","s3"),
  `7` = c("s3"),
  `8` = c("s1","s2","s3")
)

###################################################################
## Exact analytic truth for the 5 linear-Normal settings (1,2,3,4,8): the
## multivariate-Normal conditional-mean formula gives E(Y(1) | S_subset =
## s) exactly, so Delta_subset reduces to plugging the control arm's
## subset mean into that linear function.
###################################################################

analytic_truth_linear <- function(mu1, Sigma1, mu0, b1_0, b1, b0_0, b0, subset_idx) {
  k <- length(mu1)
  rest_idx <- setdiff(seq_len(k), subset_idx)
  EY1 <- b1_0 + sum(b1 * mu1)
  EY0 <- b0_0 + sum(b0 * mu0)
  Delta <- EY1 - EY0
  s_eval <- mu0[subset_idx]
  if (length(rest_idx) == 0) {
    EY1_given_s <- b1_0 + sum(b1[subset_idx] * s_eval)
  } else {
    Sig_rs <- Sigma1[rest_idx, subset_idx, drop = FALSE]
    Sig_ss <- Sigma1[subset_idx, subset_idx, drop = FALSE]
    cond_rest <- mu1[rest_idx] + Sig_rs %*% solve(Sig_ss, s_eval - mu1[subset_idx])
    EY1_given_s <- b1_0 + sum(b1[subset_idx] * s_eval) + sum(b1[rest_idx] * cond_rest)
  }
  Delta_subset <- EY1_given_s - EY0
  R <- 1 - Delta_subset / Delta
  list(R = R, Delta = Delta, Delta_subset = Delta_subset)
}

## Builds (mu1, Sigma1, mu0, b1_0, b1, b0_0, b0) for each of the 5
## linear-Normal settings, matching simulate_setting()'s parameters
## exactly. Setting 8 is represented by just its 3 true surrogates (the 77
## noise candidates have zero outcome coefficient and are independent of
## the true 3, so R_full == R_S0 exactly).
.linear_setting_params <- function(setting) {
  mk_R <- function(rho, k = 9) { R <- matrix(rho, k, k); diag(R) <- 1; R }
  if (setting %in% c(1, 2, 3)) {
    rho <- switch(as.character(setting), `1` = 0, `2` = 0.5, `3` = 0.1)
    mu1 <- c(6,6,4,1,1,4,1,1,4); mu0 <- c(5,5,2,1,1,4,1,1,4)
    sd1 <- sqrt(c(4,4,8,4,4,8,4,4,8))
    Sigma1 <- outer(sd1, sd1) * mk_R(rho)
    if (setting == 3) {
      b1 <- c(2,4,2.5,0.6,0.3,0.3,0,0,0); b1_0 <- 5.35
      b0 <- c(1.92,3.5,2.3,0.6,0.3,0.3,0,0,0); b0_0 <- 2
    } else {
      b1 <- c(2,4,2.5,0,0,0,0,0,0); b1_0 <- 5.35
      b0 <- c(1.92,3.5,2.3,0,0,0,0,0,0); b0_0 <- 2
    }
    return(list(mu1 = mu1, Sigma1 = Sigma1, mu0 = mu0, b1_0 = b1_0, b1 = b1, b0_0 = b0_0, b0 = b0))
  }
  if (setting == 4) {
    d <- 0.13
    mu1 <- c(6+d, 6+d, 4, 1,1,4,1,1,4); mu0 <- c(6, 6, 2, 1,1,4,1,1,4)
    sd1 <- sqrt(c(4,4,8,4,4,8,4,4,8))
    Sigma1 <- outer(sd1, sd1) * mk_R(0.08)
    b1 <- c(4,4,2.5,0,0,0,0,0,0); b1_0 <- 2
    b0 <- c(4,4,1.6,0,0,0,0,0,0); b0_0 <- 2
    return(list(mu1 = mu1, Sigma1 = Sigma1, mu0 = mu0, b1_0 = b1_0, b1 = b1, b0_0 = b0_0, b0 = b0))
  }
  if (setting == 8) {
    mu1 <- c(6,6,4); mu0 <- c(5,5,2)
    sd1 <- sqrt(c(4,4,8))
    Sigma1 <- outer(sd1, sd1) * mk_R(0.1, k = 3)
    b1 <- c(2,4,2.5); b1_0 <- 5.35
    b0 <- c(1.92,3.5,2.3); b0_0 <- 2
    return(list(mu1 = mu1, Sigma1 = Sigma1, mu0 = mu0, b1_0 = b1_0, b1 = b1, b0_0 = b0_0, b0 = b0))
  }
  stop("setting ", setting, " is not one of the linear-Normal settings (1,2,3,4,8)")
}

## Returns list(R_S0 = ..., R_full = ...) for one of the 5 linear-Normal
## settings, using the exact analytic formula above.
linear_truth <- function(setting) {
  s0_idx_map <- list(`1` = 1:3, `2` = 1:3, `3` = 1:3, `4` = 3, `8` = 1:3)
  p <- .linear_setting_params(setting)
  s0_idx <- s0_idx_map[[as.character(setting)]]
  full_idx <- seq_len(length(p$mu1))
  res_S0   <- analytic_truth_linear(p$mu1, p$Sigma1, p$mu0, p$b1_0, p$b1, p$b0_0, p$b0, s0_idx)
  res_full <- analytic_truth_linear(p$mu1, p$Sigma1, p$mu0, p$b1_0, p$b1, p$b0_0, p$b0, full_idx)
  list(R_S0 = res_S0$R, R_full = res_full$R)
}

###################################################################
## Plug-in Monte Carlo truth for the 3 nonlinear/heterogeneous settings
## (5,6,7): draw a large independent sample from each arm; to evaluate
## Delta for a given surrogate subset, take the large treated-arm draw,
## replace its subset columns with an independent large control-arm draw
## of those columns, and evaluate the true outcome-generating function
## directly on the resulting hybrid covariates.
###################################################################

plugin_truth_nonlinear <- function(setting, subset_cols, n_big = 1e6, seed = 1) {
  set.seed(seed)
  full_cols <- paste0("s", 1:9)

  treat_big   <- simulate_setting(n_big, n_big, setting)
  treat_big   <- treat_big[treat_big$treat == 1, ]
  control_big <- simulate_setting(n_big, n_big, setting)
  control_big <- control_big[control_big$treat == 0, ]

  Delta <- mean(treat_big$y) - mean(control_big$y)

  hyp <- treat_big
  hyp[, subset_cols] <- control_big[, subset_cols]

  y_hyp <- switch(as.character(setting),
    `5` = y1.nl(as.matrix(hyp[, full_cols])),
    `6` = y1.nl.rogue(as.matrix(hyp[, full_cols])),
    `7` = y1.redundant.nl(as.matrix(hyp[, full_cols])),
    stop("setting ", setting, " is not one of the nonlinear settings (5,6,7)")
  )
  Delta_subset <- mean(y_hyp) - mean(control_big$y)
  R <- 1 - Delta_subset / Delta
  list(R = R, Delta = Delta, Delta_subset = Delta_subset)
}

###################################################################
## .fast_nw(): Nadaraya-Watson kernel smoother (Gaussian kernel), chunked
## over both the evaluation points and the reference points so that no
## single intermediate matrix exceeds chunk_eval x chunk_x cells,
## regardless of how large eval_pts/x are. Falls back to nearest-neighbor
## for any evaluation point with zero total weight.
###################################################################

.fast_nw <- function(eval_pts, x, y, h, chunk_eval = 500, chunk_x = 20000) {
  n_eval <- length(eval_pts)
  n_x <- length(x)
  out <- numeric(n_eval)
  i <- 1
  while (i <= n_eval) {
    idx <- i:min(i + chunk_eval - 1, n_eval)
    q <- eval_pts[idx]
    wsum <- numeric(length(idx))
    wy   <- numeric(length(idx))
    j <- 1
    while (j <= n_x) {
      jdx <- j:min(j + chunk_x - 1, n_x)
      xc <- x[jdx]; yc <- y[jdx]
      w <- exp(-0.5 * (outer(q, xc, "-") / h)^2)
      wsum <- wsum + rowSums(w)
      wy   <- wy + as.vector(w %*% yc)
      j <- j + chunk_x
    }
    safe <- wsum > 0
    res <- numeric(length(idx))
    if (any(safe))  res[safe] <- wy[safe] / wsum[safe]
    if (any(!safe)) {
      for (k in which(!safe)) {
        nn <- which.min(abs(x - q[k]))
        res[k] <- y[nn]
      }
    }
    out[idx] <- res
    i <- i + chunk_eval
  }
  out
}

###################################################################
## large_sample_RQ(): large-sample approximation of R_Q, the population
## target that the linear working-model estimators (Full, RR, TARw, and
## TAR when it selects a fixed set) converge to under model
## misspecification, as opposed to the true nonlinear PTE
## (plugin_truth_nonlinear()). Mirrors the real estimator's index +
## kernel-smoothing structure (debias_ols()/make.delta.s()):
##   1. Fit OLS of Y(1) on an intercept + subset_cols, treated arm only,
##      on a large sample (n_beta) -- the population working-model
##      projection.
##   2. Form the scalar index Q = beta0 + S %*% beta in both arms, on an
##      independent sample (n_eval).
##   3. Kernel-smooth E(Y(1)|Q(1)=q) using the treated arm's (Q1,Y1),
##      evaluated at each control point's Q0 (.fast_nw()).
##   4. Delta_Q = mean(smoothed prediction) - mean(Y0); R_Q = 1 -
##      Delta_Q/Delta.
## Also reports Q-overlap diagnostics: the population-level fraction of
## Q0 falling outside Q1's observed range / 1st-99th percentile range,
## and the same fraction at a realistic n=500-per-arm finite sample
## (averaged over n_reps_finite replicates).
##
## n_eval = 5e4 (not larger): step 3 is a pairwise kernel smooth costing
## O(n_eval^2), so n_eval = 2e6 would be computationally infeasible even
## with chunked memory. 5e4 is two orders of magnitude larger than the
## real estimator's n=500, so Monte Carlo/smoothing noise at this scale
## is negligible for the precision these diagnostics need.
###################################################################

large_sample_RQ <- function(setting, subset_cols, n_beta = 2e6, n_eval = 5e4,
                             n_finite = 500, n_reps_finite = 300, seed = 1) {
  full_cols <- paste0("s", 1:9)
  set.seed(seed)

  beta_big <- simulate_setting(n_beta, n_beta, setting)
  beta_big <- beta_big[beta_big$treat == 1, ]
  fit <- lm(as.formula(paste("y ~", paste(subset_cols, collapse = " + "))), data = beta_big)
  beta <- coef(fit)

  dat <- simulate_setting(n_eval, n_eval, setting)
  treat <- dat[dat$treat == 1, ]; control <- dat[dat$treat == 0, ]
  Q1 <- beta[1] + as.matrix(treat[, subset_cols, drop = FALSE])   %*% beta[-1]
  Q0 <- beta[1] + as.matrix(control[, subset_cols, drop = FALSE]) %*% beta[-1]
  h <- bw.nrd(Q1)
  mu1_at_Q0 <- .fast_nw(as.vector(Q0), as.vector(Q1), treat$y, h)
  Delta <- mean(treat$y) - mean(control$y)
  Delta_Q <- mean(mu1_at_Q0) - mean(control$y)
  R_Q <- 1 - Delta_Q / Delta

  Q1v <- as.vector(Q1); Q0v <- as.vector(Q0)
  ov_pop <- list(
    frac_outside_range = mean(Q0v < min(Q1v) | Q0v > max(Q1v)),
    frac_outside_1_99   = mean(Q0v < quantile(Q1v, 0.01) | Q0v > quantile(Q1v, 0.99))
  )

  fracs_finite <- numeric(n_reps_finite)
  for (r in seq_len(n_reps_finite)) {
    d_fin <- simulate_setting(n_finite, n_finite, setting)
    t_fin <- d_fin[d_fin$treat == 1, ]; c_fin <- d_fin[d_fin$treat == 0, ]
    Q1_fin <- beta[1] + as.matrix(t_fin[, subset_cols, drop = FALSE]) %*% beta[-1]
    Q0_fin <- beta[1] + as.matrix(c_fin[, subset_cols, drop = FALSE]) %*% beta[-1]
    fracs_finite[r] <- mean(Q0_fin < min(Q1_fin) | Q0_fin > max(Q1_fin))
  }

  list(Delta = Delta, Delta_Q = Delta_Q, R_Q = R_Q,
       overlap_population = ov_pop,
       overlap_finite_n500 = mean(fracs_finite),
       beta = beta)
}

###################################################################
## Plotting functions
###################################################################

plot_selection_grid <- function(selection_props, setting_label, s0_cols = NULL, compact = FALSE) {
  title_text <- if (compact) setting_label else paste0("Surrogate selection probability -- ", setting_label)
  base_sz <- if (compact) 8 else 12
  txt_sz  <- if (compact) 2.2 else 3

  p <- ggplot(selection_props, aes(x = Surrogate, y = Method, fill = Selection_Prob)) +
    geom_tile(color = "white", linewidth = 0.6) +
    geom_text(aes(label = sprintf("%.2f", Selection_Prob)), size = txt_sz, color = "#D55E00", fontface = "bold") +
    scale_fill_gradientn(colors = c("#f7fbff", "#6baed6", "#08306b"), limits = c(0, 1), name = "P(selected)") +
    labs(title = title_text, x = NULL, y = NULL) +
    theme_minimal(base_size = base_sz) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 0),
          plot.title = element_text(face = "bold", size = base_sz))
  if (!is.null(s0_cols)) {
    p <- p + annotate("rect",
                       xmin = match(s0_cols, levels(factor(selection_props$Surrogate))) - 0.5,
                       xmax = match(s0_cols, levels(factor(selection_props$Surrogate))) + 0.5,
                       ymin = 0.5, ymax = length(unique(selection_props$Method)) + 0.5,
                       fill = NA, color = "#d62728", linewidth = 0.9)
  }
  if (compact) p <- p + theme(legend.position = "none")
  p
}

plot_selection_setting8 <- function(sel_true_df, sel_noise_long, compact = FALSE) {
  # sel_true_df: Method, Surrogate (s1/s2/s3), Selection_Prob
  # sel_noise_long: Method, Selection_Prob (one row per noise candidate)
  base_sz <- if (compact) 8 else 12
  txt_sz  <- if (compact) 2.0 else 2.8

  p1 <- ggplot(sel_true_df, aes(x = Surrogate, y = Selection_Prob, fill = Method)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = sprintf("%.2f", Selection_Prob)),
              position = position_dodge(width = 0.8), vjust = -0.5, size = txt_sz) +
    ylim(0, 1.12) +
    labs(title = "True surrogates (s1-s3)", x = NULL, y = "P(selected)") +
    theme_minimal(base_size = base_sz) +
    theme(legend.position = "bottom", legend.title = element_blank(),
          legend.text = element_text(size = base_sz - 1))

  n_noise <- sum(sel_noise_long$Method == sel_noise_long$Method[1])
  p2 <- ggplot(sel_noise_long, aes(x = Method, y = Selection_Prob, fill = Method)) +
    geom_boxplot(outlier.size = 0.6, width = 0.5) +
    labs(title = sprintf("%d noise candidates", n_noise), x = NULL, y = "P(selected)") +
    theme_minimal(base_size = base_sz) +
    guides(fill = "none")

  (p1 + p2) +
    plot_layout(widths = c(1, 1.3), guides = "collect") &
    theme(
      legend.position = "bottom",
      legend.direction = "horizontal"
    )
}

combined_selection_figure_1to4 <- function(panels_1to4) {
  stopifnot(length(panels_1to4) == 4)
  wrap_plots(panels_1to4, ncol = 2, nrow = 2) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 10, face = "bold"))
}

combined_selection_figure_5to8 <- function(panels_5to7, panel8) {
  stopifnot(length(panels_5to7) == 3)
  all_panels <- c(panels_5to7, list(panel8))
  wrap_plots(all_panels, ncol = 2, nrow = 2) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(size = 10, face = "bold"))
}

###################################################################
## summarize_results(): builds one table row (population truths, Est,
## ASE, ESE, P(S0^C selected), mean # selected) for one method within one
## setting. ase_vec is a numeric vector of per-replicate bootstrap SEs
## (all-NA if run.bootstrap was FALSE, giving a blank ASE cell).
###################################################################

summarize_results <- function(R_est, selection_matrix, R_full, R_S0, s0_cols, method_name, ase_vec = NULL) {
  ese <- sd(R_est, na.rm = TRUE)

  if (!is.null(selection_matrix)) {
    non_s0_cols <- setdiff(colnames(selection_matrix), s0_cols)
    if (length(non_s0_cols) > 0) {
      p_s0c_selected <- mean(rowSums(selection_matrix[, non_s0_cols, drop = FALSE]) > 0)
    } else {
      p_s0c_selected <- 0
    }
    mean_num_selected <- mean(rowSums(selection_matrix))
  } else {
    p_s0c_selected <- NA
    mean_num_selected <- NA
  }

  if (is.null(ase_vec) || all(is.na(ase_vec))) {
    ase <- NA_real_
  } else {
    ase <- mean(ase_vec, na.rm = TRUE)
  }

  data.frame(Method = method_name,
             R_full = round(R_full, 3),
             R_S0 = round(R_S0, 3),
             Est = round(mean(R_est, na.rm = TRUE), 3),
             ASE = round(ase, 3),
             ESE = round(ese, 3),
             P_S0C_selected = round(p_s0c_selected, 3),
             Mean_num_selected = round(mean_num_selected, 3))
}

###################################################################
## write_grouped_table(): writes a .tex tabular fragment (R_full, R_S0,
## Est, ASE, ESE, P(S0^C sel.), mean # sel.) for the given settings. The
## Full surrogate row shows "--" for the two selection columns, since it
## does not perform selection.
###################################################################

write_grouped_table <- function(df, settings, setting_labels, file,
                                 method_order = c("Full surrogate", "Robust reg.", "TAR (weights only)", "TAR"),
                                 method_display = c("Full", "RR", "TARw", "TAR")) {
  fmt <- function(x, digits = 3) {
    if (length(x) == 0 || is.na(x)) return("")
    formatC(x, format = "f", digits = digits)
  }

  lines <- c("\\begin{tabular}{lrrrrrrr}",
             "\\hline",
             "Method & $R_{\\bS}$ & $R_{\\bS_0}$ & $\\widehat R$ & ASE & ESE & $P(\\bS_0^C \\text{ sel.})$ & Mean \\# sel. \\\\",
             "\\hline")

  for (s in settings) {
    lab <- setting_labels[[as.character(s)]]
    lines <- c(lines, sprintf("\\multicolumn{8}{l}{\\textit{%s}} \\\\", lab), "\\hline")
    block <- df[df$Setting == s, ]
    for (m_i in seq_along(method_order)) {
      row <- block[block$Method == method_order[m_i], ]
      if (nrow(row) == 0) next
      is_full <- method_order[m_i] == "Full surrogate"
      p_str <- if (is_full) "--" else fmt(row$P_S0C_selected)
      n_str <- if (is_full) "--" else fmt(row$Mean_num_selected, digits = 2)
      lines <- c(lines, sprintf("%s & %s & %s & %s & %s & %s & %s & %s \\\\",
                                 method_display[m_i],
                                 fmt(row$R_full), fmt(row$R_S0), fmt(row$Est), fmt(row$ASE), fmt(row$ESE),
                                 p_str, n_str))
    }
    lines <- c(lines, "\\hline")
  }
  lines <- c(lines, "\\end{tabular}")
  writeLines(lines, file)
  invisible(lines)
}

###################################################################
## write_setting8fp_table(): writes the Setting 8 mean-noise-candidates-
## selected .tex fragment.
###################################################################

write_setting8fp_table <- function(sel_smoothed, sel_taronly, sel_tar, noise_cols, file) {
  fp <- data.frame(
    Method = c("Robust reg.", "TAR (weights only)", "TAR"),
    Mean_num_noise_selected_per_rep = c(
      mean(rowSums(sel_smoothed[, noise_cols, drop = FALSE])),
      mean(rowSums(sel_taronly[, noise_cols, drop = FALSE])),
      mean(rowSums(sel_tar[, noise_cols, drop = FALSE]))
    )
  )
  lines <- c("\\begin{tabular}{lc}",
             "\\hline",
             "Method & Mean number of noise candidates selected \\\\",
             "\\hline",
             sprintf("%s & %s \\\\", fp$Method, formatC(fp$Mean_num_noise_selected_per_rep, format = "f", digits = 3)),
             "\\hline",
             "\\end{tabular}")
  writeLines(lines, file)
  invisible(fp)
}
