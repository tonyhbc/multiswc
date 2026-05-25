#' Simulate longitudinal survival data with three-way treatment switching
#'
#' @description
#' `simswc()` simulates a two-arm randomized clinical-trial-like longitudinal
#' survival dataset with a three-way post-randomization treatment swich
#' process. It is intended to generate data for evaluating and demonstrating
#' multi-regime marginal structural Cox model [multimsm()].
#'
#' Follow-up is represented on a discrete visit grid
#' \eqn{k = 0, 1, \ldots, K - 1}, where `K = n_visit`. Each subject is
#' randomized at baseline to experimental treatment \eqn{A_0 = 1} or control
#' treatment \eqn{A_0 = 0}. During follow-up, subjects may experience at most
#' one absorbing treatment switch:
#'
#' \itemize{
#'   \item `ES`: randomized experimental \eqn{\rightarrow} subsequent therapy;
#'   \item `CE`: randomized control \eqn{\rightarrow} experimental crossover;
#'   \item `CS`: randomized control \eqn{\rightarrow} subsequent therapy.
#' }
#'
#' The simulator generates baseline covariates, two time-varying confounders
#' (`X` and `U`), a generic switch process `V`, switch-type indicators
#' `S_ES`, `S_CE`, and `S_CS`, optional random censoring, and a survival outcome
#' from a piecewise exponential hazard model. The output includes both
#' subject-level wide-format data and counting-process long-format data suitable
#' for Cox modeling and for the package's multi-regime MSM workflow.
#'
#' @details
#' Let \eqn{L_i} denote the vector of selected baseline covariates named in
#' `base_cov`. The simulator first generates six possible baseline covariates:
#'
#' \itemize{
#'   \item `L1`: Binary age group indicator,
#'     \eqn{L_1 \sim \mathrm{Bernoulli}(0.3)}.
#'   \item `L2`: Three-level risk group taking values `0`, `1`, and `2` with
#'     probabilities `0.60`, `0.35`, and `0.05`, respectively.
#'   \item `L3`: Binary baseline metastasis/progression-like indicator,
#'     \eqn{L_3 \sim \mathrm{Bernoulli}(0.1)}.
#'   \item `L4`: Continuous tumor-size-like marker,
#'     \eqn{L_4 \sim N(60, 35^2)}.
#'   \item `L5`: Binary biomarker indicator,
#'     \eqn{L_5 \sim \mathrm{Bernoulli}(0.8)}.
#'   \item `L6`: Binary prior-treatment indicator,
#'     \eqn{L_6 \sim \mathrm{Bernoulli}(0.6)}.
#' }
#'
#' Only the variables included in `base_cov` are used in the treatment,
#' confounder, switch, censoring, and hazard linear predictors, and only these
#' selected baseline covariates are returned in the simulated datasets.
#'
#' Baseline treatment is assigned by stratified randomization within the
#' interaction strata formed by `base_cov`; within each stratum, approximately
#' `trt_prob` of subjects are assigned to \eqn{A_0 = 1}. Because `base_cov`
#' defines both the model covariate vector and the randomization strata, using a
#' continuous variable such as `L4` may create many small strata. In practical
#' simulation studies, it is often preferable to stratify on categorical or binary
#' baseline covariates, for example `base_cov = c("L1", "L3")`.
#'
#' The time-varying binary confounder `X` is initialized as \eqn{X_0 = 0} and is
#' made absorbing by construction. For \eqn{k = 1, \ldots, K - 1},
#'
#' \deqn{
#' \Pr(\widetilde X_{ik}=1 \mid H_{i,k-1}) =
#' \operatorname{expit}\{\eta_0 + \eta_L^\top L_i +
#' \eta_A A_{i,k-1} + \eta_X X_{i,k-1}\},
#' }
#'
#' and the stored value is
#'
#' \deqn{X_{ik} = \max(X_{i,k-1}, \widetilde X_{ik}).}
#'
#' The continuous time-varying confounder `U` is initialized as
#' \eqn{U_0 \sim N(0,1)}. For \eqn{k = 1, \ldots, K - 1},
#'
#' \deqn{
#' U_{ik} \sim N(\mu_{ik}, 1), \quad
#' \mu_{ik} =
#' \alpha_0 + \alpha_L^\top L_i +
#' \alpha_A A_{i,k-1} + \alpha_U U_{i,k-1}.
#' }
#'
#' The generic switch indicator `V` is initialized as \eqn{V_0 = 0}. Among
#' subjects who have not yet switched, the probability of a new switch at visit
#' \eqn{k} is
#'
#' \deqn{
#' \Pr(V_{ik}=1 \mid V_{i,k-1}=0, H_{ik}) =
#' \operatorname{expit}\{\delta_0 + \delta_L^\top L_i +
#' \delta_A A_{i,k-1} + \delta_X X_{ik} + \delta_U U_{ik}\}.
#' }
#'
#' Once a subject switches, `V` remains equal to 1. The switch destination is then
#' assigned as follows:
#'
#' \itemize{
#'   \item if \eqn{A_{i0}=1}, the subject switches to `ES`; `S_ES` is set to 1
#'     from that visit onward and `A` is set to 0 thereafter to represent being
#'     off the original experimental-treatment indicator and on subsequent
#'     therapy;
#'   \item if \eqn{A_{i0}=0}, the subject switches either to `CE` or to `CS`;
#'     under `CE`, `S_CE` is set to 1 and `A` is set to 1 thereafter; under `CS`,
#'     `S_CS` is set to 1 and `A` remains 0 thereafter.
#' }
#'
#' If `param_select = NULL`, the control-arm switch destination is chosen with
#' fixed probability `p_pick_CE` for `CE` and `1 - p_pick_CE` for `CS`. If
#' `param_select` is supplied, the probability of `CE` among newly switching
#' control-arm subjects is
#'
#' \deqn{
#' \Pr(\mathrm{CE}_{ik}=1 \mid \mathrm{new\ switch}, A_{i0}=0) =
#' \operatorname{expit}\{\xi_0 + \xi_X X_{ik} + \xi_U U_{ik}\}.
#' }
#'
#' The event time is generated from a piecewise exponential model. For interval
#' \eqn{[k, k+1)}, where \eqn{k = 0,\ldots,K-1}, the internal column index is
#' \eqn{m = k + 1}, and the interval-specific hazard is
#'
#' \deqn{
#' \lambda_{ik} =
#' \{\lambda_0 \exp(0.1m)\}
#' \exp\{\beta_A A_{ik} + \beta_L^\top L_i +
#' \beta_X X_{ik} + \beta_U U_{ik} +
#' \beta_{ES} S^{ES}_{ik} +
#' \beta_{CE} S^{CE}_{ik} +
#' \beta_{CS} S^{CS}_{ik}\}.
#' }
#'
#' Conditional on the interval covariates, an exponential waiting time with rate
#' \eqn{\lambda_{ik}} is drawn. If the waiting time is less than 1, the event is
#' placed inside that interval at exact continuous time \eqn{k + T^\ast};
#' otherwise simulation proceeds to the next interval.
#'
#' If `param_cens` is supplied, random right censoring is generated in discrete
#' time. At the end of interval \eqn{k}, among subjects not yet randomly
#' censored,
#'
#' \deqn{
#' \Pr(C_{i,k+1}=1 \mid C_{ik}=0, H_{ik}) =
#' \operatorname{expit}\{\gamma_0 + \gamma_L^\top L_i +
#' \gamma_A A_{ik} + \gamma_X X_{ik} + \gamma_U U_{ik}\}.
#' }
#'
#' If `param_cens = NULL`, no additional random censoring is generated. In all
#' cases, the simulator also applies an administrative censoring time
#' `C.a = study_end - T_e`, where `T_e` is sampled uniformly from the integers
#' `0:floor(n_visit / 3)`. The observed time is the minimum of event time,
#' random censoring time, and administrative censoring time.
#'
#' The long-format output uses start-stop counting-process intervals. The
#' terminal event time is kept as a continuous time when the event occurs inside
#' an interval. If `lagk = TRUE`, the long-format data also include
#' `Alag1`-`Alag3`, `Xlag1`-`Xlag3`, `Ulag1`-`Ulag3`, and one-step leads
#' `Xnext1` and `Unext1`.
#'
#' For the revised triplet-input [multimsm()] interface, the simulation output
#' can be converted by defining `rand` as the subject-level baseline treatment
#' `A.0`, `cross = S_CE`, and `subseq = pmax(S_ES, S_CS)`.
#'
#' @param n Integer sample size. The current implementation requires an integer
#'   greater than 49.
#' @param n_visit Integer number of planned discrete visits, \eqn{K}. The visit
#'   indices in the wide data are `0, ..., n_visit - 1`, and the long-format
#'   intervals are initially `[0,1), [1,2), ...`.
#' @param base_cov Character vector specifying which baseline covariates are
#'   included in the simulator's linear predictors and used to define strata for
#'   baseline randomization. Must be a non-empty subset of
#'   `c("L1", "L2", "L3", "L4", "L5", "L6")`. Coefficients in all parameter
#'   vectors must follow this exact order.
#' @param trt_prob Numeric scalar in `(0, 1)`. Target probability of assignment
#'   to baseline experimental treatment \eqn{A_0 = 1} within each stratum formed
#'   by `base_cov`.
#' @param param_tdcov Numeric vector of length `length(base_cov) + 3` giving the
#'   coefficients \eqn{\eta} for the binary time-varying confounder model for
#'   `X`. If `p = length(base_cov)`, the required order is:
#'   \itemize{
#'     \item `[1]`: intercept \eqn{\eta_0};
#'     \item `[2:(p+1)]`: baseline-covariate coefficients
#'       \eqn{\eta_L}, in the same order as `base_cov`;
#'     \item `[p+2]`: coefficient for prior treatment \eqn{A_{k-1}}
#'       \eqn{\eta_A};
#'     \item `[p+3]`: coefficient for prior binary confounder
#'       \eqn{X_{k-1}}, \eqn{\eta_X}.
#'   }
#' @param param_tdcov2 Numeric vector of length `length(base_cov) + 3` giving
#'   the coefficients \eqn{\alpha} for the continuous time-varying confounder
#'   mean model for `U`. If `p = length(base_cov)`, the required order is:
#'   \itemize{
#'     \item `[1]`: intercept \eqn{\alpha_0};
#'     \item `[2:(p+1)]`: baseline-covariate coefficients
#'       \eqn{\alpha_L}, in the same order as `base_cov`;
#'     \item `[p+2]`: coefficient for prior treatment \eqn{A_{k-1}},
#'       \eqn{\alpha_A};
#'     \item `[p+3]`: coefficient for prior continuous confounder
#'       \eqn{U_{k-1}}, \eqn{\alpha_U}.
#'   }
#' @param param_sw Numeric vector of length `length(base_cov) + 4` giving the
#'   coefficients \eqn{\delta} for the generic treatment-switch model. If
#'   `p = length(base_cov)`, the required order is:
#'   \itemize{
#'     \item `[1]`: intercept \eqn{\delta_0};
#'     \item `[2:(p+1)]`: baseline-covariate coefficients
#'       \eqn{\delta_L}, in the same order as `base_cov`;
#'     \item `[p+2]`: coefficient for prior treatment \eqn{A_{k-1}},
#'       \eqn{\delta_A};
#'     \item `[p+3]`: coefficient for current binary confounder \eqn{X_k},
#'       \eqn{\delta_X};
#'     \item `[p+4]`: coefficient for current continuous confounder \eqn{U_k},
#'       \eqn{\delta_U}.
#'   }
#' @param param_haz Numeric vector of length `length(base_cov) + 7` specifying
#'   the piecewise exponential survival model. If `p = length(base_cov)`, the
#'   required order is:
#'   \itemize{
#'     \item `[1]`: baseline hazard scale \eqn{\lambda_0}. This is a hazard
#'       scale, not a log-hazard coefficient, and should be positive;
#'     \item `[2]`: log hazard ratio for current treatment \eqn{A_k},
#'       \eqn{\beta_A};
#'     \item `[3:(p+2)]`: baseline-covariate log hazard ratios
#'       \eqn{\beta_L}, in the same order as `base_cov`;
#'     \item `[p+3]`: log hazard ratio for current binary confounder
#'       \eqn{X_k}, \eqn{\beta_X};
#'     \item `[p+4]`: log hazard ratio for current continuous confounder
#'       \eqn{U_k}, \eqn{\beta_U};
#'     \item `[p+5]`: log hazard ratio for `ES` status,
#'       \eqn{\beta_{ES}};
#'     \item `[p+6]`: log hazard ratio for `CE` status,
#'       \eqn{\beta_{CE}};
#'     \item `[p+7]`: log hazard ratio for `CS` status,
#'       \eqn{\beta_{CS}}.
#'   }
#' @param param_cens Optional numeric vector of length `length(base_cov) + 4`
#'   giving the coefficients \eqn{\gamma} for the discrete-time random
#'   censoring model. If `NULL`, no random censoring is generated beyond the
#'   administrative censoring mechanism. If `p = length(base_cov)`, the required
#'   order is:
#'   \itemize{
#'     \item `[1]`: intercept \eqn{\gamma_0};
#'     \item `[2:(p+1)]`: baseline-covariate coefficients
#'       \eqn{\gamma_L}, in the same order as `base_cov`;
#'     \item `[p+2]`: coefficient for current treatment \eqn{A_k},
#'       \eqn{\gamma_A};
#'     \item `[p+3]`: coefficient for current binary confounder \eqn{X_k},
#'       \eqn{\gamma_X};
#'     \item `[p+4]`: coefficient for current continuous confounder \eqn{U_k},
#'       \eqn{\gamma_U}.
#'   }
#' @param study_end Numeric scalar giving the nominal administrative study-end
#'   time. Must be between `2` and `n_visit`. The actual administrative
#'   censoring time returned as `C.a` is `study_end - T_e`, where `T_e` is a
#'   randomly sampled integer from `0:floor(n_visit / 3)`.
#' @param lagk Logical. If `TRUE`, add lagged versions of `A`, `X`, and `U` up
#'   to three prior visits and one-step future versions of `X` and `U` to the
#'   long-format output. If `FALSE`, these auxiliary lag/lead variables are not
#'   added.
#' @param param_select Optional numeric vector of length 3 controlling the
#'   `CE` versus `CS` destination model among newly switching control-arm
#'   subjects. The required order is:
#'   \itemize{
#'     \item `[1]`: intercept \eqn{\xi_0};
#'     \item `[2]`: coefficient for current \eqn{X_k}, \eqn{\xi_X};
#'     \item `[3]`: coefficient for current \eqn{U_k}, \eqn{\xi_U}.
#'   }
#'   If `NULL`, `p_pick_CE` is used instead.
#' @param p_pick_CE Numeric scalar in `[0, 1]`. Probability that a newly
#'   switching control-arm subject switches to experimental therapy (`CE`)
#'   rather than subsequent therapy (`CS`) when `param_select = NULL`.
#' @param true_hr Logical. If `TRUE`, compute an approximate Monte Carlo
#'   benchmark for the sustained always-experimental versus always-control
#'   comparison using `truth_n` simulated subjects per sustained regime. This
#'   benchmark ignores switching and uses the non-switch components of the
#'   longitudinal confounder and hazard models.
#' @param truth_n Non-negative integer Monte Carlo sample size used only when
#'   `true_hr = TRUE`. Set `true_hr = FALSE` or `truth_n = 0` to skip this
#'   potentially time-consuming benchmark.
#'
#' @return A list with components:
#' \itemize{
#'   \item `dat`: Subject-level wide-format data. It contains `id`, selected
#'     baseline covariates, visit-indexed matrices expanded into columns for
#'     `X`, `U`, `A`, `V`, `S_ES`, `S_CE`, and `S_CS`, observed event/censoring
#'     information (`T.obs`, `D.obs`), observed switch time `T.w`, administrative
#'     censoring time `C.a`, and observed switch type `type_sw` (`1 = ES`,
#'     `2 = CE`, `3 = CS`, `NA = no observed switch before exit`).
#'   \item `dat_long`: Counting-process long-format data derived from `dat`,
#'     with one row per subject-interval while the subject is under observation.
#'     It includes `t.start`, `t.stop`, `event`, `C`, selected baseline
#'     covariates, `A`, `V`, `X`, `U`, and `S_ES`/`S_CE`/`S_CS`. If `lagk = TRUE`,
#'     it also includes the lag/lead variables described above.
#'   \item `C`: An `n` by `n_visit + 1` matrix containing the cumulative
#'     censoring process after combining random and administrative censoring.
#'   \item `stats`: Named summary statistics: `init_trt`, the proportion
#'     randomized to experimental treatment; `switched`, the proportion with an
#'     observed switch before exit; `prop_ES`, `prop_CE`, and `prop_CS`, the
#'     proportions of switch types among switchers; and `prop_event`, the
#'     observed event proportion.
#'   \item `true_sdiff`: If `true_hr = TRUE` and `truth_n > 0`, a list containing
#'     an approximate sustained-regime benchmark, including a fitted
#'     `survival::survfit` object, time grid `t`, survival curves `surv0` and
#'     `surv1`, an approximate Cox hazard ratio `cHR`, survival difference
#'     `surv_diff`, and the simulated benchmark datasets `dat_A0` and `dat_A1`.
#'     Otherwise `NULL`.
#' }
#'
#' @examples
#' set.seed(123)
#'
#' sim_obj <- simswc(
#'   n = 200,
#'   n_visit = 8,
#'   base_cov = c("L1", "L3"),
#'   trt_prob = 0.5,
#'   param_tdcov  = c(-1.5, 0.3, 0.3, -0.2, 0.5),
#'   param_tdcov2 = c( 0.05, 0.2, 0.2, -0.2, 0.7),
#'   param_sw     = c(-3.0, 0.4, 0.4, -0.3, 0.3, 0.8),
#'   param_haz    = c(0.1, log(0.8), log(1.1), log(1.1),
#'                    log(1.4), log(1.4), log(1.3), log(1), log(0.9)),
#'   param_cens   = NULL,
#'   param_select = c(0, 0.5, 0.25),
#'   lagk = TRUE,
#'   true_hr = FALSE
#' )
#'
#' head(sim_obj$dat_long)
#' sim_obj$stats
#'
#' ## Convert to the triplet interface used by multimsm().
#' dat <- sim_obj$dat_long
#' dat$rand <- stats::ave(dat$A, dat$id, FUN = function(x) rep(x[1], length(x)))
#' dat$cross <- dat$S_CE
#' dat$subseq <- pmax(dat$S_ES, dat$S_CS)
#'
#' table(dat$rand, dat$cross, dat$subseq)
#'
#' @references
#' Keogh RH, Seaman SR, Gran JM, Vansteelandt S. Simulating longitudinal data
#' from marginal structural models using the additive hazard model.
#' *Biometrical Journal*. 2021;63:1526-1541.
#'
#' @export
simswc <- function(n,
                       n_visit,
                       base_cov = paste0("L", 1:6),
                       trt_prob = 0.5,
                       param_tdcov,
                       param_tdcov2,
                       param_sw,
                       param_haz,
                       param_cens = NULL,
                       study_end = n_visit,
                       lagk = FALSE,
                       param_select = NULL,
                       p_pick_CE = 0.5,
                       true_hr = TRUE,
                       truth_n = 100000) {

  expit <- function(x) stats::plogis(x)
  n_basecov <- length(base_cov)

  if (!is.numeric(n) || length(n) != 1L || is.na(n) || n <= 49 || n != as.integer(n)) {
    stop("`n` must be an integer greater than 100.")
  }
  if (!is.numeric(n_visit) || length(n_visit) != 1L || is.na(n_visit) ||
      n_visit <= 1 || n_visit != as.integer(n_visit)) {
    stop("`n_visit` must be an integer greater than 1.")
  }
  if (!is.character(base_cov) || length(base_cov) < 1L || any(!base_cov %in% paste0("L", 1:6))) {
    stop("`base_cov` must be a non-empty character vector drawn from L1 to L6.")
  }
  if (!is.numeric(trt_prob) || length(trt_prob) != 1L || is.na(trt_prob) ||
      trt_prob <= 0 || trt_prob >= 1) {
    stop("`trt_prob` must lie strictly between 0 and 1.")
  }
  if (!is.numeric(p_pick_CE) || length(p_pick_CE) != 1L || is.na(p_pick_CE) ||
      p_pick_CE < 0 || p_pick_CE > 1) {
    stop("`p_pick_CE` must lie between 0 and 1.")
  }
  if (!is.null(param_select)) {
    if (!is.numeric(param_select) || length(param_select) != 3L || anyNA(param_select)) {
      stop("`param_select` must be NULL or a numeric vector of length 3.")
    }
  }
  if (!is.numeric(study_end) || length(study_end) != 1L || is.na(study_end) ||
      study_end < 2 || study_end > n_visit) {
    stop("`study_end` must be between 2 and `n_visit`.")
  }
  if (length(param_tdcov) != n_basecov + 3L) {
    stop("`param_tdcov` must have length length(base_cov) + 3.")
  }
  if (length(param_tdcov2) != n_basecov + 3L) {
    stop("`param_tdcov2` must have length length(base_cov) + 3.")
  }
  if (length(param_sw) != n_basecov + 4L) {
    stop("`param_sw` must have length length(base_cov) + 4.")
  }
  if (length(param_haz) != n_basecov + 7L) {
    stop("`param_haz` must have length length(base_cov) + 7.")
  }
  if (!is.null(param_cens) && length(param_cens) != n_basecov + 4L) {
    stop("`param_cens` must have length length(base_cov) + 4 when supplied.")
  }
  if (!is.logical(lagk) || length(lagk) != 1L || is.na(lagk)) {
    stop("`lagk` must be TRUE or FALSE.")
  }
  if (!is.logical(true_hr) || length(true_hr) != 1L || is.na(true_hr)) {
    stop("`true_hr` must be TRUE or FALSE.")
  }
  if (!is.numeric(truth_n) || length(truth_n) != 1L || is.na(truth_n) || truth_n < 0) {
    stop("`truth_n` must be a non-negative integer.")
  }
  truth_n <- as.integer(truth_n)

  L <- data.frame(
    id = seq_len(n),
    L1 = stats::rbinom(n, 1, 0.3),
    L2 = sample(0:2, n, replace = TRUE, prob = c(0.6, 0.35, 0.05)),
    L3 = stats::rbinom(n, 1, 0.1),
    L4 = stats::rnorm(n, 60, 35),
    L5 = stats::rbinom(n, 1, 0.8),
    L6 = stats::rbinom(n, 1, 0.6)
  )

  grp <- interaction(L[, base_cov, drop = FALSE], drop = TRUE, lex.order = TRUE)
  A0 <- integer(n)
  idx_split <- split(seq_len(n), grp)
  for (ii in idx_split) {
    n_i <- length(ii)
    n_treated <- round(trt_prob * n_i)
    A0[ii] <- sample(c(rep.int(1L, n_treated), rep.int(0L, n_i - n_treated)), size = n_i, replace = FALSE)
  }
  L$A.0 <- A0

  A <- X <- V <- matrix(0L, nrow = n, ncol = n_visit)
  U <- matrix(0, nrow = n, ncol = n_visit)
  S_ES <- S_CE <- S_CS <- matrix(0L, nrow = n, ncol = n_visit)
  colnames(A) <- paste0("A.", 0:(n_visit - 1L))
  colnames(X) <- paste0("X.", 0:(n_visit - 1L))
  colnames(U) <- paste0("U.", 0:(n_visit - 1L))
  colnames(V) <- paste0("V.", 0:(n_visit - 1L))
  colnames(S_ES) <- paste0("S_ES.", 0:(n_visit - 1L))
  colnames(S_CE) <- paste0("S_CE.", 0:(n_visit - 1L))
  colnames(S_CS) <- paste0("S_CS.", 0:(n_visit - 1L))

  A[, 1L] <- A0
  X[, 1L] <- 0L
  U[, 1L] <- stats::rnorm(n, 0, 1)

  Lmat <- as.matrix(L[, base_cov, drop = FALSE])

  for (k in 2:n_visit) {
    A[, k] <- A[, k - 1L]

    p_X <- expit(
      param_tdcov[1L] +
        Lmat %*% param_tdcov[2:(1L + n_basecov)] +
        param_tdcov[2L + n_basecov] * A[, k - 1L] +
        param_tdcov[3L + n_basecov] * X[, k - 1L]
    )
    X[, k] <- pmax(X[, k - 1L], stats::rbinom(n, 1L, p_X))

    mu_U <- as.numeric(
      param_tdcov2[1L] +
        Lmat %*% param_tdcov2[2:(1L + n_basecov)] +
        param_tdcov2[2L + n_basecov] * A[, k - 1L] +
        param_tdcov2[3L + n_basecov] * U[, k - 1L]
    )
    U[, k] <- stats::rnorm(n, mean = mu_U, sd = 1)

    p_sw <- expit(
      param_sw[1L] +
        Lmat %*% param_sw[2:(1L + n_basecov)] +
        param_sw[2L + n_basecov] * A[, k - 1L] +
        param_sw[3L + n_basecov] * X[, k] +
        param_sw[4L + n_basecov] * U[, k]
    )
    risk_switch <- (V[, k - 1L] == 0L)
    new_switch <- risk_switch & (stats::rbinom(n, 1L, p_sw) == 1L)
    V[, k] <- pmax(V[, k - 1L], as.integer(new_switch))

    idx_switch <- which(new_switch)
    if (length(idx_switch) > 0L) {
      for (i in idx_switch) {
        if (A[i, 1L] == 1L) {
          S_ES[i, k:n_visit] <- 1L
          A[i, k:n_visit] <- 0L
        } else {
          p_ce_i <- if (is.null(param_select)) {
            p_pick_CE
          } else {
            expit(param_select[1L] + param_select[2L] * X[i, k] + param_select[3L] * U[i, k])
          }
          if (stats::runif(1L) < p_ce_i) {
            S_CE[i, k:n_visit] <- 1L
            A[i, k:n_visit] <- 1L
          } else {
            S_CS[i, k:n_visit] <- 1L
          }
        }
      }
    }
  }

  type_sw <- rep(NA_integer_, n)
  type_sw[rowSums(S_ES) > 0L] <- 1L
  type_sw[rowSums(S_CE) > 0L] <- 2L
  type_sw[rowSums(S_CS) > 0L] <- 3L

  C_r_mat <- matrix(0L, nrow = n, ncol = n_visit + 1L)
  colnames(C_r_mat) <- paste0("C.", 0:n_visit)

  if (!is.null(param_cens)) {
    cens_state <- rep(FALSE, n)
    for (k in 2:(n_visit + 1L)) {
      idx_uncens <- which(!cens_state)
      if (length(idx_uncens) == 0L) {
        break
      }
      kk <- k - 1L
      p_C <- expit(
        param_cens[1L] +
          Lmat[idx_uncens, , drop = FALSE] %*% param_cens[2:(1L + n_basecov)] +
          param_cens[2L + n_basecov] * A[idx_uncens, kk] +
          param_cens[3L + n_basecov] * X[idx_uncens, kk] +
          param_cens[4L + n_basecov] * U[idx_uncens, kk]
      )
      new_c <- stats::rbinom(length(idx_uncens), 1L, p_C)
      if (any(new_c == 1L)) {
        who <- idx_uncens[new_c == 1L]
        C_r_mat[who, k] <- 1L
        cens_state[who] <- TRUE
      }
    }
  }

  C_r <- apply(C_r_mat, 1L, function(x) {
    idx <- which(x == 1L)
    if (length(idx) > 0L) min(idx) - 1L else NA_integer_
  })

  C_a_true <- study_end
  T_e <- sample(0:floor(n_visit / 3), n, replace = TRUE)
  C_a <- C_a_true - T_e
  C_obs <- pmin(C_r, C_a, na.rm = TRUE)

  C <- C_r_mat
  for (i in seq_len(n)) {
    if (is.na(C_r[i]) || C_a[i] < C_r[i]) {
      C[i, ] <- 0L
      if (C_a[i] >= 1L && C_a[i] <= n_visit) {
        C[i, C_a[i] + 1L] <- 1L
      }
    }
  }
  C <- t(apply(C, 1L, cummax))

  T_obs <- rep(NA_real_, n)
  D_obs <- rep(0L, n)
  for (i in seq_len(n)) {
    for (k in seq_len(n_visit)) {
      haz <- (param_haz[1L] * exp(0.1 * k)) * exp(
        param_haz[2L] * A[i, k] +
          sum(as.numeric(L[i, base_cov, drop = TRUE]) * param_haz[3:(2L + n_basecov)]) +
          param_haz[3L + n_basecov] * X[i, k] +
          param_haz[4L + n_basecov] * U[i, k] +
          param_haz[5L + n_basecov] * S_ES[i, k] +
          param_haz[6L + n_basecov] * S_CE[i, k] +
          param_haz[7L + n_basecov] * S_CS[i, k]
      )
      tstar <- -log(stats::runif(1L)) / haz
      if (is.na(T_obs[i]) && tstar < 1) {
        T_obs[i] <- (k - 1L) + tstar
        D_obs[i] <- 1L
        break
      }
    }

    if (is.na(T_obs[i])) {
      if (C_obs[i] >= n_visit) {
        T_obs[i] <- n_visit
        D_obs[i] <- 0L
      } else {
        T_obs[i] <- C_obs[i]
        D_obs[i] <- 0L
      }
    } else if (C_obs[i] < T_obs[i]) {
      T_obs[i] <- C_obs[i]
      D_obs[i] <- 0L
    }
  }

  T_w <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    idx <- which(V[i, ] == 1L)
    if (length(idx) > 0L) {
      sw_time <- min(idx) - 1L
      if (sw_time < T_obs[i]) T_w[i] <- sw_time
    }
  }
  type_sw[is.na(T_w)] <- NA_integer_

  dat <- data.frame(
    id = seq_len(n),
    L[, base_cov, drop = FALSE],
    X,
    U,
    A,
    V,
    type_sw = type_sw,
    S_ES,
    S_CE,
    S_CS,
    T.obs = T_obs,
    D.obs = D_obs,
    T.w = T_w,
    C.a = C_a,
    check.names = FALSE
  )

  extend_long <- function(dat, n_visit, lagk) {
    varying_list <- list(
      paste0("A.", 0:(n_visit - 1L)),
      paste0("V.", 0:(n_visit - 1L)),
      paste0("X.", 0:(n_visit - 1L)),
      paste0("U.", 0:(n_visit - 1L)),
      paste0("S_ES.", 0:(n_visit - 1L)),
      paste0("S_CE.", 0:(n_visit - 1L)),
      paste0("S_CS.", 0:(n_visit - 1L))
    )

    dat_long <- stats::reshape(
      dat,
      varying = varying_list,
      v.names = c("A", "V", "X", "U", "S_ES", "S_CE", "S_CS"),
      timevar = "time",
      times = 0:(n_visit - 1L),
      direction = "long",
      idvar = "id"
    )
    rownames(dat_long) <- NULL
    dat_long <- dat_long[order(dat_long$id, dat_long$time), , drop = FALSE]
    dat_long$t.start <- dat_long$time
    dat_long$t.stop <- dat_long$time + 1
    dat_long <- dat_long[dat_long$t.start < dat_long$T.obs, , drop = FALSE]
    dat_long$event <- as.integer(dat_long$t.start <= dat_long$T.obs & dat_long$T.obs < dat_long$t.stop & dat_long$D.obs == 1L)
    idx_event <- which(dat_long$t.stop > dat_long$T.obs & dat_long$D.obs == 1L)
    if (length(idx_event) > 0L) dat_long$t.stop[idx_event] <- dat_long$T.obs[idx_event]
    dat_long$C <- as.integer(dat_long$event == 0L & dat_long$D.obs == 0L & dat_long$t.stop == dat_long$T.obs)

    if (isTRUE(lagk)) {
      n_long <- nrow(dat_long)
      dat_long$Alag1 <- integer(n_long)
      dat_long$Xlag1 <- integer(n_long)
      dat_long$Ulag1 <- numeric(n_long)
      dat_long$Alag2 <- integer(n_long)
      dat_long$Xlag2 <- integer(n_long)
      dat_long$Ulag2 <- numeric(n_long)
      dat_long$Alag3 <- integer(n_long)
      dat_long$Xlag3 <- integer(n_long)
      dat_long$Ulag3 <- numeric(n_long)
      dat_long$Xnext1 <- integer(n_long)
      dat_long$Unext1 <- numeric(n_long)

      split_idx <- split(seq_len(nrow(dat_long)), dat_long$id)
      add_lags <- function(x, k, default) {
        if (length(x) <= k) rep(default, length(x)) else c(rep(default, k), utils::head(x, -k))
      }
      add_lead <- function(x, k, default) {
        if (length(x) <= k) rep(default, length(x)) else c(utils::tail(x, -k), rep(default, k))
      }

      for (ii in split_idx) {
        dat_long$Alag1[ii] <- add_lags(dat_long$A[ii], 1L, 0L)
        dat_long$Xlag1[ii] <- add_lags(dat_long$X[ii], 1L, 0L)
        dat_long$Ulag1[ii] <- add_lags(dat_long$U[ii], 1L, 0)
        dat_long$Alag2[ii] <- add_lags(dat_long$A[ii], 2L, 0L)
        dat_long$Xlag2[ii] <- add_lags(dat_long$X[ii], 2L, 0L)
        dat_long$Ulag2[ii] <- add_lags(dat_long$U[ii], 2L, 0)
        dat_long$Alag3[ii] <- add_lags(dat_long$A[ii], 3L, 0L)
        dat_long$Xlag3[ii] <- add_lags(dat_long$X[ii], 3L, 0L)
        dat_long$Ulag3[ii] <- add_lags(dat_long$U[ii], 3L, 0)
        dat_long$Xnext1[ii] <- add_lead(dat_long$X[ii], 1L, 0L)
        dat_long$Unext1[ii] <- add_lead(dat_long$U[ii], 1L, 0)
      }
    }

    dat_long
  }

  dat_long <- extend_long(dat = dat, n_visit = n_visit, lagk = lagk)

  stats <- c(
    init_trt = mean(A0),
    switched = mean(!is.na(type_sw)),
    prop_ES = mean(type_sw == 1L, na.rm = TRUE),
    prop_CE = mean(type_sw == 2L, na.rm = TRUE),
    prop_CS = mean(type_sw == 3L, na.rm = TRUE),
    prop_event = mean(D_obs == 1L)
  )

  true_sdiff <- NULL
  if (isTRUE(true_hr) && truth_n > 0L) {
    true_sdiff_fun <- function(n, n_visit, base_cov, param_tdcov, param_tdcov2, param_haz, t_pts) {
      expit <- function(x) stats::plogis(x)
      n_basecov <- length(base_cov)
      L <- data.frame(
        id = seq_len(n),
        L1 = stats::rbinom(n, 1, 0.3),
        L2 = sample(0:2, n, replace = TRUE, prob = c(0.6, 0.35, 0.05)),
        L3 = stats::rbinom(n, 1, 0.1),
        L4 = stats::rnorm(n, 60, 35),
        L5 = stats::rbinom(n, 1, 0.8),
        L6 = stats::rbinom(n, 1, 0.6)
      )
      Lmat <- as.matrix(L[, base_cov, drop = FALSE])

      sim_regime <- function(A_regime) {
        X <- matrix(0L, nrow = n, ncol = n_visit)
        U <- matrix(0, nrow = n, ncol = n_visit)
        X[, 1L] <- 0L
        U[, 1L] <- stats::rnorm(n, 0, 1)
        for (k in 2:n_visit) {
          p_X <- expit(
            param_tdcov[1L] +
              Lmat %*% param_tdcov[2:(1L + n_basecov)] +
              param_tdcov[2L + n_basecov] * A_regime +
              param_tdcov[3L + n_basecov] * X[, k - 1L]
          )
          X[, k] <- pmax(X[, k - 1L], stats::rbinom(n, 1L, p_X))

          mu_U <- as.numeric(
            param_tdcov2[1L] +
              Lmat %*% param_tdcov2[2:(1L + n_basecov)] +
              param_tdcov2[2L + n_basecov] * A_regime +
              param_tdcov2[3L + n_basecov] * U[, k - 1L]
          )
          U[, k] <- stats::rnorm(n, mean = mu_U, sd = 1)
        }

        T_obs <- rep(NA_real_, n)
        D_obs <- rep(0L, n)
        for (k in seq_len(n_visit)) {
          haz <- (param_haz[1L] * exp(0.1 * k)) * exp(
            param_haz[2L] * A_regime +
              Lmat %*% param_haz[3:(2L + n_basecov)] +
              param_haz[3L + n_basecov] * X[, k] +
              param_haz[4L + n_basecov] * U[, k]
          )
          t_new <- -log(stats::runif(n)) / haz
          T_obs <- ifelse(is.na(T_obs) & t_new < 1, (k - 1L) + t_new, T_obs)
        }
        D_obs <- ifelse(is.na(T_obs), 0L, 1L)
        T_obs <- ifelse(is.na(T_obs), n_visit, T_obs)
        data.frame(id = seq_len(n), T.obs = T_obs, D.obs = D_obs, A = A_regime)
      }

      dat_A1 <- sim_regime(1)
      dat_A0 <- sim_regime(0)
      dat_comb <- rbind(dat_A1, dat_A0)
      fit_hr <- survival::coxph(survival::Surv(T.obs, D.obs) ~ A, data = dat_comb)
      cHR <- exp(stats::coef(fit_hr)[[1L]])
      surv_fit <- survival::survfit(survival::Surv(T.obs, D.obs) ~ A, data = dat_comb)
      surv_sum <- summary(surv_fit, times = t_pts, extend = TRUE)
      grp <- surv_sum$strata
      times_all <- surv_sum$time
      surv_all <- surv_sum$surv
      idx0 <- which(grp == "A=0")
      idx1 <- which(grp == "A=1")
      surv0 <- stats::approx(times_all[idx0], surv_all[idx0], xout = t_pts,
                             method = "constant", yleft = 1,
                             yright = utils::tail(surv_all[idx0], 1L), ties = "ordered")$y
      surv1 <- stats::approx(times_all[idx1], surv_all[idx1], xout = t_pts,
                             method = "constant", yleft = 1,
                             yright = utils::tail(surv_all[idx1], 1L), ties = "ordered")$y
      list(
        surv_fit = surv_fit,
        t = t_pts,
        surv0 = surv0,
        surv1 = surv1,
        cHR = cHR,
        surv_diff = surv0 - surv1,
        dat_A0 = dat_A0,
        dat_A1 = dat_A1
      )
    }

    true_sdiff <- true_sdiff_fun(
      n = truth_n,
      n_visit = n_visit,
      base_cov = base_cov,
      param_tdcov = param_tdcov,
      param_tdcov2 = param_tdcov2,
      param_haz = param_haz[1:(length(param_haz) - 3L)],
      t_pts = seq(0, n_visit, by = 0.05)
    )
  }

  list(
    dat = dat,
    dat_long = dat_long,
    C = C,
    stats = round(stats, 3),
    true_sdiff = true_sdiff
  )
}
