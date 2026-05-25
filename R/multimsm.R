#' Fit a Multi-regime Marginal Structural Cox Model for Three-way Treatment Switching
#'
#' @description
#' `multimsm()` fits a marginal structural Cox model with categorical treatment
#' for randomized clinical trials or trial-like longitudinal data with multi-way
#' treatment switching, assuming participants underwent regular longitudinal
#' follow-up visits. This method assumes a three-way treatment switch scenario:
#' control participants may crossover to experimental group and both control
#' and experimental participants may switch to a subsequent therapy.
#' `multimsm()` takes a longitudinal survival dataset with 3 treatment process
#' variables as primary inputs:
#'
#' * `rand`: __baseline__ randomized treatment arm. `1` denotes randomized
#'   experimental treatment and `0` denotes randomized control/SOC.
#' * `cross`: __time-varying__ control-to-experimental crossover indicator (1/0).
#' * `subseq`: __time-varying__ subsequent therapy initiation indicator (1/0).
#'
#' Internally, `multimsm()` converts this triplet of treatment trajectory into
#' the five-level _time-varying treatment regime_:
#'
#' \deqn{G(t) \in \{C, E, CE, CS, ES\}}
#'
#' where `C` and `E` denote __sustained control__ and __sustained experimental__
#' treatment regime at t, `CE` denotes __control-to-experimental crossover__, `CS`
#' denotes on __control-to-subsequent switch__ regime status at t, and `ES` denotes
#' on __experimental-to-subsequent switch__ regime status at t.
#'
#' The downstream estimator is a multinomial-regime stabilized-weight marginal
#' structural Cox model with the time-varying regime as predictor:
#'
#' \deqn{\lambda_{T^{\overline g}}(t)=\lambda_0(t)\exp\{\beta_e\mathbb{I}({G(t)=E}) +
#' \beta_{ce} \mathbb{I}({G(t)=CE}) + \beta_{cs} \mathbb{I}({G(t)=CS}) +
#' \beta_{es} \mathbb{I}({G(t)=ES})\}}
#'
#' The stabilized treatment-regime weight (S-IPTW) at each
#' post-baseline interval is the probability of the observed current regime under
#' a stabilizer model divided by the corresponding probability under a denominator
#' model. The cumulative product of these row-specific ratios is used as the
#' inverse probability weight. Optional ordinary censoring weights may also be
#' multiplied in, if informative censoring is present.
#'
#' @details
#' The treatment-switching process is assumed to follow the package's
#' five-regime structure G(t):
#'
#' \deqn{
#' G_i(t)=
#' \begin{cases}
#' C,  & R_i=0,\ Crs_i(t)=0,\ Subs_i(t)=0,\\
#' E,  & R_i=1,\ Crs_i(t)=0,\ Subs_i(t)=0,\\
#' CE, & R_i=0,\ Crs_i(t)=1,\ Subs_i(t)=0,\\
#' CS, & R_i=0,\ Crs_i(t)=0,\ Subs_i(t)=1,\\
#' ES, & R_i=1,\ Crs_i(t)=0,\ Subs_i(t)=1.
#' \end{cases}
#' }
#'
#' `R(t)` is `rand` randomized arm, `Crs(t)` is `cross` the absorbing crossover
#' status, and `Subs(t)` is `subseq` the absorbing subsequent initiation status.
#' Users may supply `cross` and `subseq` with either sustained absorbing
#' status indicator such as (`0, 0, 1, 1, 1`), or as switch-initiation action
#' indicators such as (`0, 0, 1, 0, 0`) - both acceptable.
#'
#' The current method allows _only one switch type_ per subject's follow up. Thus,
#' control-arm subjects may remain in `C`, transition to `CE`, or transition to `CS`;
#' experimental-arm subjects may remain in `E` or transition to `ES`.
#'
#' The outcome model uses generated `regime` as the exposure. By default
#' __`C` is the reference regime__, so the coefficient for `regimeE` estimates the
#' marginal sustained experimental-versus-control contrast aimed by the package's
#' primary no-switch estimand. The remaining coefficients describe the distinct
#' switched-regime contrasts averaging switch times relative to sustained control.
#'
#' @param dat_long A `data.frame` in long-format start-stop format with one row
#'   per subject follow-up interval.
#' @param id Character. Subject identifier variable name.
#' @param tstart Character. Interval start-time variable name.
#' @param tstop Character. Interval stop-time variable name.
#' @param event Character. Event indicator variable name (`1`: event, `0`: censor).
#' @param cens Optional character Ordinary right-censoring indicator used
#'   only if `ipcw_mod` is supplied. This is not the treatment-switching process.
#' @param rand Character. Baseline randomization variable name (`1`: randomized
#'   experimental treatment, `0`: randomized control). The value must be constant
#'   within subject over time.
#' @param cross Character. _Time-varying_ control-to-experimental crossover
#'   indicator. It may be one-time pulse at crossover time (...,0,1,0,0,...) or
#'   absorbing (...,0,1,1,1,...). It must be structurally 0 always for subjects
#'   randomized to experimental treatment (assumed no EC switch).
#' @param subseq Character. _Time-varying_ subsequent therapy initiation indicator.
#'   It may be pulse or absorbing.
#' @param base_cov Optional character vector of baseline covariates to include in
#'   the default numerator model when `iptw_num` is `NULL`.
#' @param iptw_num Optional _right-hand-side_ formula or character string for
#'   the stabilized numerator multinomial regime model. If `NULL`, the default is
#'   `~ regime_lag + factor(visit) + < base_cov >`.
#' @param iptw_den Required _right-hand-side formula_ or character string for
#'   the denominator multinomial regime model. This should typically include
#'   `regime_lag`, visit or time terms, baseline covariates, and time-varying
#'   confounders affected by prior treatment and predictive of subsequent
#'   switching.
#' @param ipcw_mod Optional right-hand-side formula or character string for an
#'   ordinary censoring model. If supplied, `cens` must also be supplied.
#' @param wt_trunc Optional numeric in `(0.5, 1)`. If supplied, the
#'   final combined inverse probability weight is truncated at the
#'   `(1 - wt_trunc)` and `wt_trunc` empirical quantiles before the Cox model is
#'   fit. For example, `wt_trunc = 0.95` caps the lower and upper tails at the
#'   5th and 95th percentiles. The default `NULL` performs no quantile
#'   truncation.
#' @param prob_bounds Numeric length-2 vector of lower and upper probability
#'   bounds used to stabilize extreme predicted probabilities away from 0 and 1.
#' @param normalize_weights Logical. If `TRUE`, normalize the final truncated
#'   weights to have mean 1.
#' @param robust Logical. If `TRUE`, fit the Cox model with `cluster(id)` to use
#'   the standard robust sandwich variance estimator.
#' @param check_inputs Logical. If `TRUE`, enforce structural checks for the
#'   five-regime treatment-switching system.
#' @param trace_multinom Logical. Passed to [nnet::multinom()].
#' @param maxit Integer maximum number of iterations for the multinomial models.
#'
#' @return An object of class `"multimsm"`; a list with components:
#'
#' \itemize{
#'   \item `coef_table` A `data.frame` of full estimated coefficient table.
#'   \item `coef` Log-hazard ratio estimates for treatment regimes against control.
#'   \item `hr` Hazard ratio estimates for treatment regimes against control.
#'   \item `hr_ci` Wald 95% confidence intervals for the hazard ratios.
#'   \item `fit` The fitted weighted Cox model from [survival::coxph()].
#'   \item `p.val` P values of regime effect estimates.
#'   \item `dat_long` Augmented long-format data containing the derived variables,
#'     including switch indicators, and final estimated weights
#'     (`.w_use`/`.w_use_trunc`).
#'   \item `regime_models` A list containing the fitted numerator and denominator
#'     mulinomial regime models.
#'   \item `censor_model` The fitted censoring model when `ipcw_mod` is supplied;
#'     othrwise `NULL`.
#'   \item `diagnostics` A list containing untruncated final-weight quantiles,
#'     final switch-regime counts/proportions, cohort size, and truncation bounds.
#'   \item `call` The matched function call.
#' }
#'
#' @references
#' Robins JM, Hernan MA, Brumback B. Marginal structural models and causal
#' inference in epidemiology. *Epidemiology*. 2000;11(5):550-560.
#'
#' Suarez D, Haro JM, Novick D, et al. Marginal structural models for multiple
#' treatment comparisons: an application to antipsychotic treatment for
#' schizophrenia. *Journal of Clinical Epidemiology*. 2008;61(6):525-530.
#'
#' @examples
#' set.seed(123)
#'
#' sim_obj <- simswc(
#'   n = 300,
#'   n_visit = 8,
#'   base_cov = c("L1", "L3"),
#'   trt_prob = 0.5,
#'   param_tdcov = c(-1.5, 0.3, 0.3, -0.2, 0.5),
#'   param_tdcov2 = c(0.05, 0.2, 0.2, -0.2, 0.7),
#'   param_sw = c(-3, 0.4, 0.4, -0.3, 0.3, 0.8),
#'   param_haz = c(0.1, log(0.8), log(1.1), log(1.1),
#'                 log(1.4), log(1.4), log(1.3), log(1), log(0.9)),
#'   param_cens = NULL,
#'   param_select = c(0, 0.5, 0.25),
#'   lagk = TRUE,
#'   true_hr = FALSE
#' )
#'
#' # Convert to the required treatment switch indicator triplet
#' dat <- sim_obj$dat_long
#' dat$rand <- stats::ave(dat$A, dat$id, FUN = function(x) rep(x[1], length(x)))
#' dat$cross <- dat$S_CE
#' dat$subseq <- pmax(dat$S_ES, dat$S_CS)
#'
#' fit <- multimsm(
#'   dat_long = dat, id = "id", tstart = "t.start", tstop = "t.stop",
#'   event = "event", rand = "rand", cross = "cross", subseq = "subseq",
#'   base_cov = c("L1", "L3"),
#'   iptw_num = ~ regime_lag + factor(visit) + L1 + L3,
#'   iptw_den = ~ regime_lag + factor(visit) + L1 + L3 + X + U + Alag1,
#'   wt_trunc = 0.95
#' )
#'
#' fit
#'
#' @export
multimsm <- function(dat_long, id = "id", tstart = "t.start", tstop = "t.stop",
                     event = "event", cens = NULL, rand = "rand", cross = "cross",
                     subseq = "subseq", base_cov = NULL, iptw_num = NULL, iptw_den,
                     ipcw_mod = NULL, wt_trunc = NULL, prob_bounds = c(1e-6, 1 - 1e-6),
                     normalize_weights = TRUE, robust = TRUE, check_inputs = TRUE,
                     trace_multinom = FALSE, maxit = 200) {

  cl <- match.call()
  regime_levels <- c("C", "E", "CE", "CS", "ES")

  if (missing(iptw_den) || is.null(iptw_den)) {
    stop("`iptw_den` is required and must be an RHS-only formula or string.",
         call. = FALSE)
  }
  .check_required_namespace("nnet")
  .check_required_namespace("survival")
  .check_scalar_character(id, "id")
  .check_scalar_character(tstart, "tstart")
  .check_scalar_character(tstop, "tstop")
  .check_scalar_character(event, "event")
  .check_scalar_character(rand, "rand")
  .check_scalar_character(cross, "cross")
  .check_scalar_character(subseq, "subseq")
  if (!is.null(cens)) .check_scalar_character(cens, "cens")
  .check_wt_trunc(wt_trunc)
  .check_probability_bounds(prob_bounds)
  if (!is.logical(normalize_weights) || length(normalize_weights) != 1L ||
      is.na(normalize_weights)) {
    stop("`normalize_weights` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(robust) || length(robust) != 1L || is.na(robust)) {
    stop("`robust` must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.logical(check_inputs) || length(check_inputs) != 1L ||
      is.na(check_inputs)) {
    stop("`check_inputs` must be TRUE or FALSE.", call. = FALSE)
  }

  df <- as.data.frame(dat_long)
  .require_columns(df, c(id, tstart, tstop, event, rand, cross, subseq),
                   context = "`dat_long`")
  if (!is.null(cens)) .require_columns(df, cens, context = "`dat_long`")

  df$.id <- df[[id]]
  df$.tstart <- df[[tstart]]
  df$.tstop <- df[[tstop]]
  df$.event <- df[[event]]
  df$.rand_raw <- df[[rand]]
  df$.cross_raw <- df[[cross]]
  df$.subseq_raw <- df[[subseq]]
  df$.cens <- if (!is.null(cens)) df[[cens]] else 0L

  df <- df[order(df$.id, df$.tstart, df$.tstop), , drop = FALSE]

  .validate_survival_columns(df)
  df$.rand <- .as_binary01(df$.rand_raw, "rand")
  df$.cross_raw <- .as_binary01(df$.cross_raw, "cross")
  df$.subseq_raw <- .as_binary01(df$.subseq_raw, "subseq")
  df$.event <- .as_binary01(df$.event, "event")
  if (anyNA(df$.rand)) {
    stop("`rand` cannot contain missing values.", call. = FALSE)
  }
  if (anyNA(df$.cross_raw)) {
    stop("`cross` cannot contain missing values; code no crossover as 0.",
         call. = FALSE)
  }
  if (anyNA(df$.subseq_raw)) {
    stop("`subseq` cannot contain missing values; code no subsequent therapy as 0.",
         call. = FALSE)
  }
  if (anyNA(df$.event)) {
    stop("`event` cannot contain missing values.", call. = FALSE)
  }
  if (!is.null(cens)) {
    df$.cens <- .as_binary01(df$.cens, "cens")
    if (anyNA(df$.cens)) {
      stop("`cens` cannot contain missing values when `ipcw_mod` is used.",
           call. = FALSE)
    }
  }

  df$.row_in_id <- unsplit(
    lapply(split(df$.id, df$.id, drop = TRUE), seq_along),
    df$.id
  )
  df$.visit <- df$.tstart
  df$visit <- df$.visit

  rand_by_id <- split(df$.rand, df$.id, drop = TRUE)
  nonconstant_rand <- names(rand_by_id)[vapply(
    rand_by_id,
    function(x) length(unique(stats::na.omit(x))) > 1L,
    logical(1L)
  )]
  if (length(nonconstant_rand) > 0L) {
    stop("`rand` must be constant within subject. Non-constant IDs include: ",
         paste(utils::head(nonconstant_rand, 5L), collapse = ", "),
         if (length(nonconstant_rand) > 5L) ", ..." else "",
         call. = FALSE)
  }
  df$.rand <- unsplit(
    lapply(rand_by_id, function(x) rep(x[which(!is.na(x))[1L]], length(x))),
    df$.id
  )

  df$.cross_abs <- .cummax_by_id(df$.cross_raw, df$.id)
  df$.subseq_abs <- .cummax_by_id(df$.subseq_raw, df$.id)

  if (isTRUE(check_inputs)) {
    .check_triplet_structure(df)
  }

  df <- .derive_regime_from_triplet(df, regime_levels = regime_levels)

  ## User-facing columns for formulas and diagnostics.
  df$regime <- df$.regime
  df$regime_lag <- df$.regime_lag
  df$cross_abs <- df$.cross_abs
  df$subseq_abs <- df$.subseq_abs
  df$S_CE_derived <- df$.S_CE
  df$S_CS_derived <- df$.S_CS
  df$S_ES_derived <- df$.S_ES

  if (is.null(iptw_num)) {
    iptw_num <- .default_numerator_formula(base_cov, df)
  } else if (!is.null(base_cov)) {
    base_cov <- as.character(base_cov)
    .require_columns(df, base_cov, context = "`dat_long`")
  }

  idx_trt <- which(df$.row_in_id > 1L)
  if (length(idx_trt) == 0L) {
    stop("No post-baseline rows found; cannot fit a longitudinal regime model.",
         call. = FALSE)
  }

  trt_df <- df[idx_trt, , drop = FALSE]
  trt_df$.regime <- droplevels(trt_df$.regime)
  trt_df$.regime_lag <- droplevels(trt_df$.regime_lag)
  trt_df$regime <- trt_df$.regime
  trt_df$regime_lag <- trt_df$.regime_lag

  num_fm <- .build_formula_from_rhs(iptw_num, lhs = ".regime")
  den_fm <- .build_formula_from_rhs(iptw_den, lhs = ".regime")

  num_mod <- nnet::multinom(
    formula = num_fm,
    data = trt_df,
    trace = trace_multinom,
    maxit = maxit
  )
  den_mod <- nnet::multinom(
    formula = den_fm,
    data = trt_df,
    trace = trace_multinom,
    maxit = maxit
  )

  p_num <- rep(1, nrow(df))
  p_den <- rep(1, nrow(df))
  p_num[idx_trt] <- .observed_multinom_prob(
    model = num_mod,
    newdata = trt_df,
    observed = trt_df$.regime,
    prob_bounds = prob_bounds
  )
  p_den[idx_trt] <- .observed_multinom_prob(
    model = den_mod,
    newdata = trt_df,
    observed = trt_df$.regime,
    prob_bounds = prob_bounds
  )

  df$.p_num <- p_num
  df$.p_den <- p_den
  df$.w_row_regime <- p_num / p_den
  df$.w_row_regime[df$.row_in_id == 1L] <- 1
  df$.w_row_regime <- pmin(
    pmax(df$.w_row_regime, prob_bounds[1L]),
    1 / prob_bounds[1L]
  )
  df$.wt_regime <- .cumprod_by_id(df$.w_row_regime, df$.id)

  cens_mod <- NULL
  df$.wt_cens <- 1
  if (!is.null(ipcw_mod)) {
    if (is.null(cens)) {
      stop("`cens` must be provided when `ipcw_mod` is supplied.",
           call. = FALSE)
    }
    cens_fm <- .build_formula_from_rhs(ipcw_mod, lhs = ".cens")
    cens_mod <- stats::glm(cens_fm, family = stats::binomial(), data = df)
    p_cens <- stats::predict(cens_mod, newdata = df, type = "response")
    p_cens <- pmin(pmax(p_cens, prob_bounds[1L]), prob_bounds[2L])
    p_observed <- pmax(1 - p_cens, prob_bounds[1L])
    df$.wt_cens <- 1 / pmax(.lagged_cumprod_by_id(p_observed, df$.id),
                            prob_bounds[1L])
  }

  df$.w_comb <- df$.wt_regime * df$.wt_cens
  df$.w_comb[!is.finite(df$.w_comb)] <- NA_real_
  if (all(is.na(df$.w_comb))) {
    stop("All final weights are missing or non-finite.", call. = FALSE)
  }
  med_w <- stats::median(df$.w_comb[is.finite(df$.w_comb)], na.rm = TRUE)
  df$.w_comb[!is.finite(df$.w_comb) | is.na(df$.w_comb)] <- med_w

  df$.w_use_untruncated <- df$.w_comb
  if (!is.null(wt_trunc)) {
    qlo <- as.numeric(stats::quantile(
      df$.w_use_untruncated,
      probs = 1 - wt_trunc,
      na.rm = TRUE,
      names = FALSE
    ))
    qhi <- as.numeric(stats::quantile(
      df$.w_use_untruncated,
      probs = wt_trunc,
      na.rm = TRUE,
      names = FALSE
    ))
    df$.w_use_trunc <- pmin(pmax(df$.w_use_untruncated, qlo), qhi)
  } else {
    qlo <- NA_real_
    qhi <- NA_real_
    df$.w_use_trunc <- df$.w_use_untruncated
  }

  if (isTRUE(normalize_weights)) {
    df$.w_use <- df$.w_use_trunc * nrow(df) / sum(df$.w_use_trunc, na.rm = TRUE)
  } else {
    df$.w_use <- df$.w_use_trunc
  }
  df$.wt_use <- df$.w_use

  switch_stats <- .patient_final_regime_summary(df, regime_levels = regime_levels)
  if (switch_stats$count["C"] < 20L || switch_stats$count["E"] < 20L) {
    warning("!Attention! Possible weak support: too few always E/ always C subject counts (< 20).",
            call. = FALSE)
  }

  df$.regime <- stats::relevel(df$.regime, ref = regime_levels[1L])
  df$regime <- stats::relevel(df$regime, ref = regime_levels[1L])

  if (isTRUE(robust)) {
    cox_fm <- stats::as.formula(
      "survival::Surv(.tstart, .tstop, .event) ~ regime + cluster(.id)"
    )
  } else {
    cox_fm <- stats::as.formula(
      "survival::Surv(.tstart, .tstop, .event) ~ regime"
    )
  }

  fit <- survival::coxph(cox_fm, data = df, weights = df$.w_use)

  beta_hat <- stats::coef(fit)
  vc <- stats::vcov(fit)
  se <- sqrt(diag(vc))
  z_val <- beta_hat / se
  p_val <- 2 * stats::pnorm(abs(z_val), lower.tail = FALSE)
  hr <- exp(beta_hat)
  ci <- cbind(
    lower = exp(beta_hat - 1.96 * se),
    upper = exp(beta_hat + 1.96 * se)
  )
  coef_table <- data.frame(
    logHR = as.numeric(beta_hat),
    HR = as.numeric(hr),
    lower = as.numeric(ci[, "lower"]),
    upper = as.numeric(ci[, "upper"]),
    p.val = as.numeric(p_val),
    .pstar = .p_stars(p_val),
    check.names = FALSE,
    row.names = .display_regime_coef_names(names(beta_hat))
  )

  wq_probs <- c(0, 0.05, 0.5, 0.95, 1)
  wq <- stats::quantile(df$.w_use, probs = wq_probs, na.rm = TRUE)

  out <- list(
    fit = fit,
    coef = beta_hat,
    hr = hr,
    hr_ci = ci,
    coef_table = coef_table,
    p.val = p_val,
    dat_long = df,
    regime_models = list(numerator = num_mod, denominator = den_mod),
    censor_model = cens_mod,
    diagnostics = list(
      weight_quantiles = wq,
      switch_count = switch_stats$count,
      switch_prop = switch_stats$prop,
      n_ids = switch_stats$n_ids,
      wt_trunc = wt_trunc,
      wt_trunc_bounds = c(lower = qlo, upper = qhi)
    ),
    call = cl
  )
  class(out) <- "multimsm"

  out
}

#' Print a fitted multi-regime marginal structural Cox model
#'
#' @param x An object returned by [multimsm()].
#' @param digits Number of decimal places for coefficient and weight summaries.
#' @param ... Currently unused.
#'
#' @return Invisibly returns `x`.
#' @export
print.multimsm <- function(x, digits = 3, ...) {
  if (!inherits(x, "multimsm")) {
    return(NextMethod())
  }

  title <- "Multi-regime marginal structural Cox model"
  cat(strrep("=", nchar(title)), "\n", sep = "")
  cat(title, "\n", sep = "")
  cat(strrep("=", nchar(title)), "\n", sep = "")

  switch_count <- x$diagnostics$switch_count
  switch_prop <- x$diagnostics$switch_prop
  if (!is.null(switch_count) && !is.null(switch_prop)) {
    n_ids <- x$diagnostics$n_ids
    if (is.null(n_ids) || length(n_ids) != 1L || is.na(n_ids)) {
      n_ids <- sum(as.integer(switch_count), na.rm = TRUE)
    }

    cat(sprintf("Switch summary (N = %d):\n", as.integer(n_ids)))
    labels <- c(
      C = "# of always C",
      E = "# of always E",
      CE = "# of CE switch",
      CS = "# of CS switch",
      ES = "# of ES switch"
    )
    for (nm in names(labels)) {
      cnt <- if (nm %in% names(switch_count)) as.integer(switch_count[[nm]]) else 0L
      pr <- if (nm %in% names(switch_prop)) as.numeric(switch_prop[[nm]]) else 0
      cat(sprintf("   *  %s: %d (%s%%)\n",
                  labels[[nm]], cnt,
                  formatC(100 * pr, format = "f", digits = 1L)))
    }
    cat(strrep("~", nchar(title)), "\n", sep = "")
  }

  cat("Hazard ratios (HR) by regime\n")
  if (!is.null(x$coef_table)) {
    tab <- x$coef_table
  } else if (length(x$hr) > 0L) {
    se <- sqrt(diag(stats::vcov(x$fit)))
    z_val <- x$coef / se
    p_val <- 2 * stats::pnorm(abs(z_val), lower.tail = FALSE)
    tab <- data.frame(
      logHR = as.numeric(x$coef),
      HR = as.numeric(x$hr),
      lower = as.numeric(x$hr_ci[, "lower"]),
      upper = as.numeric(x$hr_ci[, "upper"]),
      p.val = as.numeric(p_val),
      check.names = FALSE,
      row.names = .display_regime_coef_names(names(x$coef))
    )
    tab$.pstar <- .p_stars(p_val)
  } else {
    tab <- NULL
  }

  if (is.null(tab) || nrow(tab) == 0L) {
    cat("No non-reference regime coefficients were estimated.\n")
  } else {
    tab_print <- tab

    pstar <- rep("", nrow(tab_print))
    if (".pstar" %in% names(tab_print)) {
      pstar <- as.character(tab_print$.pstar)
      tab_print$.pstar <- NULL
    }

    wanted <- c("logHR", "HR", "lower", "upper", "p.val")
    keep <- intersect(wanted, names(tab_print))
    tab_print <- tab_print[, keep, drop = FALSE]

    for (cc in intersect(wanted, names(tab_print))) {
      tab_print[[cc]] <- formatC(as.numeric(tab_print[[cc]]),
                                 format = "f", digits = digits)
    }

    mat <- as.matrix(tab_print)
    mat <- cbind(mat, pstar)
    colnames(mat)[ncol(mat)] <- ""
    print(noquote(mat), quote = FALSE, right = TRUE)
  }

  if (!is.null(x$diagnostics$weight_quantiles)) {
    cat(strrep("~", nchar(title)), "\n", sep = "")
    cat("Weight quantiles:\n")

    wq <- x$diagnostics$weight_quantiles
    wq <- as.numeric(wq)
    q_names <- names(x$diagnostics$weight_quantiles)
    if (!is.null(q_names)) {
      names(wq) <- q_names
    }

    keep_q <- c("0%", "5%", "50%", "95%", "100%")
    if (!is.null(names(wq)) && all(keep_q %in% names(wq))) {
      wq <- wq[keep_q]
    } else {
      probs <- c(0, 0.05, 0.50, 0.95, 1.00)
      if (!is.null(x$dat_long$.w_use)) {
        wq <- stats::quantile(x$dat_long$.w_use, probs = probs, na.rm = TRUE)
      } else {
        wq <- stats::quantile(wq, probs = probs, na.rm = TRUE)
      }
      names(wq) <- keep_q
    }

    wq_print <- matrix(
      formatC(as.numeric(wq), format = "f", digits = digits),
      nrow = 1L
    )
    colnames(wq_print) <- keep_q
    rownames(wq_print) <- "S.IPTW"
    print(noquote(wq_print), quote = FALSE, right = TRUE)
  }

  cat(strrep("-", nchar(title)), "\n", sep = "")
  cat("A `multimsm` object\n")
  invisible(x)
}

.display_regime_coef_names <- function(x) {
  out <- as.character(x)
  out <- sub("^regime", "regime.", out)
  out
}

.p_stars <- function(p) {
  out <- rep("", length(p))
  out[!is.na(p) & p < 0.05] <- "*"
  out[!is.na(p) & p < 0.01] <- "**"
  out[!is.na(p) & p < 0.001] <- "***"
  out
}

.patient_final_regime_summary <- function(df, regime_levels) {
  ord <- order(df$.id, df$.tstart, df$.tstop)
  d <- df[ord, , drop = FALSE]
  last_idx <- !duplicated(d$.id, fromLast = TRUE)
  last_regime <- as.character(d$.regime[last_idx])
  n_ids <- length(last_regime)

  count <- stats::setNames(rep(0L, length(regime_levels)), regime_levels)
  obs_count <- table(factor(last_regime, levels = regime_levels))
  count[names(obs_count)] <- as.integer(obs_count)
  prop <- count / n_ids

  list(count = count, prop = prop, n_ids = n_ids)
}

.check_required_namespace <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop("Package '", pkg, "' is required.", call. = FALSE)
  }
  invisible(TRUE)
}

.check_scalar_character <- function(x, name) {
  if (!is.character(x) || length(x) != 1L || is.na(x) || !nzchar(x)) {
    stop("`", name, "` must be a non-missing character scalar.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.require_columns <- function(data, cols, context = "data") {
  missing_cols <- setdiff(cols, names(data))
  if (length(missing_cols) > 0L) {
    stop("Missing required column(s) in ", context, ": ",
         paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  invisible(TRUE)
}

.check_wt_trunc <- function(wt_trunc) {
  if (is.null(wt_trunc)) {
    return(invisible(TRUE))
  }
  if (!is.numeric(wt_trunc) || length(wt_trunc) != 1L || anyNA(wt_trunc) ||
      !is.finite(wt_trunc) || wt_trunc <= 0.5 || wt_trunc >= 1) {
    stop("`wt_trunc` must be NULL or a single numeric value with ",
         "0.5 < wt_trunc < 1. For example, use `wt_trunc = 0.95` ",
         "to truncate at the 5th and 95th percentiles.", call. = FALSE)
  }
  invisible(TRUE)
}

.check_probability_bounds <- function(prob_bounds) {
  if (!is.numeric(prob_bounds) || length(prob_bounds) != 2L ||
      anyNA(prob_bounds) || prob_bounds[1L] <= 0 || prob_bounds[2L] >= 1 ||
      prob_bounds[1L] >= prob_bounds[2L]) {
    stop("`prob_bounds` must be a numeric length-2 vector with ",
         "0 < prob_bounds[1] < prob_bounds[2] < 1.", call. = FALSE)
  }
  invisible(TRUE)
}

.validate_survival_columns <- function(df) {
  if (anyNA(df$.id)) {
    stop("`id` contains missing values.", call. = FALSE)
  }
  if (!is.numeric(df$.tstart) || !is.numeric(df$.tstop)) {
    stop("`tstart` and `tstop` must be numeric.", call. = FALSE)
  }
  if (anyNA(df$.tstart) || anyNA(df$.tstop)) {
    stop("`tstart` and `tstop` cannot contain missing values.", call. = FALSE)
  }
  if (any(df$.tstop <= df$.tstart)) {
    stop("Every interval must satisfy `tstop > tstart`.", call. = FALSE)
  }
  invisible(TRUE)
}

.as_binary01 <- function(x, name) {
  if (is.factor(x)) x <- as.character(x)
  if (is.logical(x)) x <- as.integer(x)
  if (is.character(x)) {
    x_trim <- trimws(tolower(x))
    x_trim[x_trim %in% c("false", "f", "no", "n")] <- "0"
    x_trim[x_trim %in% c("true", "t", "yes", "y")] <- "1"
    suppressWarnings(x_num <- as.numeric(x_trim))
  } else {
    suppressWarnings(x_num <- as.numeric(x))
  }
  bad <- !(is.na(x_num) | x_num %in% c(0, 1))
  if (any(bad)) {
    stop("`", name, "` must be coded as 0/1, TRUE/FALSE, or yes/no.",
         call. = FALSE)
  }
  as.integer(x_num)
}

.cummax_by_id <- function(x, id) {
  unsplit(lapply(split(x, id, drop = TRUE), function(z) {
    z[is.na(z)] <- 0L
    cummax(as.integer(z))
  }), id)
}

.cumprod_by_id <- function(x, id) {
  unsplit(lapply(split(x, id, drop = TRUE), cumprod), id)
}

.lagged_cumprod_by_id <- function(x, id) {
  unsplit(lapply(split(x, id, drop = TRUE), function(z) {
    cumprod(c(1, utils::head(z, -1L)))
  }), id)
}

.check_triplet_structure <- function(df) {
  bad_exp_cross <- unique(df$.id[df$.rand == 1L & df$.cross_abs == 1L])
  if (length(bad_exp_cross) > 0L) {
    stop("`cross` is a control-to-experimental crossover indicator and must be ",
         "0 for subjects randomized to experimental treatment. Violating IDs ",
         "include: ",
         paste(utils::head(bad_exp_cross, 5L), collapse = ", "),
         if (length(bad_exp_cross) > 5L) ", ..." else "",
         call. = FALSE)
  }

  bad_both <- unique(df$.id[df$.cross_abs == 1L & df$.subseq_abs == 1L])
  if (length(bad_both) > 0L) {
    stop("The current five-regime model allows only one switch type. ",
         "Subjects cannot have both `cross` and `subseq` after absorption. ",
         "Such paths would require an additional regime such as CES. ",
         "Violating IDs include: ",
         paste(utils::head(bad_both, 5L), collapse = ", "),
         if (length(bad_both) > 5L) ", ..." else "",
         call. = FALSE)
  }
  invisible(TRUE)
}

.derive_regime_from_triplet <- function(df, regime_levels) {
  df$.S_CE <- as.integer(df$.rand == 0L & df$.cross_abs == 1L &
                           df$.subseq_abs == 0L)
  df$.S_CS <- as.integer(df$.rand == 0L & df$.cross_abs == 0L &
                           df$.subseq_abs == 1L)
  df$.S_ES <- as.integer(df$.rand == 1L & df$.cross_abs == 0L &
                           df$.subseq_abs == 1L)

  regime <- ifelse(df$.rand == 1L, "E", "C")
  regime[df$.S_CE == 1L] <- "CE"
  regime[df$.S_CS == 1L] <- "CS"
  regime[df$.S_ES == 1L] <- "ES"

  bad_regime <- setdiff(unique(stats::na.omit(regime)), regime_levels)
  if (length(bad_regime) > 0L) {
    stop("Derived regimes not included in `regime_levels`: ",
         paste(bad_regime, collapse = ", "), call. = FALSE)
  }

  df$.regime <- factor(regime, levels = regime_levels)
  if (anyNA(df$.regime)) {
    stop("Unable to derive a valid treatment regime for every row.",
         call. = FALSE)
  }

  lag_regime_by_id <- function(x) {
    x_chr <- as.character(x)
    if (length(x_chr) == 0L) return(character(0L))
    c(x_chr[1L], utils::head(x_chr, -1L))
  }
  df$.regime_lag <- unsplit(
    lapply(split(df$.regime, df$.id, drop = TRUE), lag_regime_by_id),
    df$.id
  )
  df$.regime_lag <- factor(df$.regime_lag, levels = regime_levels)

  df
}

.default_numerator_formula <- function(base_cov, df) {
  if (!is.null(base_cov) && length(base_cov) > 0L) {
    base_cov <- as.character(base_cov)
    .require_columns(df, base_cov, context = "`dat_long`")
    rhs <- paste(c("regime_lag", "factor(visit)", base_cov), collapse = " + ")
  } else {
    rhs <- "regime_lag + factor(visit)"
  }
  stats::as.formula(paste("~", rhs))
}

.build_formula_from_rhs <- function(rhs, lhs) {
  if (inherits(rhs, "formula")) {
    rhs_chr <- paste(deparse(rhs), collapse = "")
  } else {
    rhs_chr <- paste(rhs, collapse = "")
  }
  rhs_chr <- trimws(rhs_chr)
  if (!nzchar(rhs_chr)) {
    stop("Model formula cannot be empty.", call. = FALSE)
  }
  if (grepl("^~", rhs_chr)) {
    out <- paste0(lhs, " ", rhs_chr)
  } else if (grepl("~", rhs_chr)) {
    out <- rhs_chr
  } else {
    out <- paste0(lhs, " ~ ", rhs_chr)
  }
  stats::as.formula(out)
}

.observed_multinom_prob <- function(model, newdata, observed, prob_bounds) {
  pr <- stats::predict(model, newdata = newdata, type = "probs")
  lev <- model$lev

  if (is.null(dim(pr))) {
    if (length(lev) == 2L) {
      pr <- cbind(1 - pr, pr)
      colnames(pr) <- lev
    } else {
      pr <- matrix(pr, ncol = length(lev))
      colnames(pr) <- lev
    }
  }
  pr <- as.matrix(pr)
  if (is.null(colnames(pr))) colnames(pr) <- lev

  obs_chr <- as.character(observed)
  out <- rep(NA_real_, length(obs_chr))
  for (j in seq_along(obs_chr)) {
    sj <- obs_chr[j]
    if (!is.na(sj) && sj %in% colnames(pr)) {
      out[j] <- pr[j, sj]
    }
  }

  out <- pmin(pmax(out, prob_bounds[1L]), prob_bounds[2L])
  out
}
