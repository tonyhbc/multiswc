# `multiswc`: Multi-regime marginal structural models for treatment switching in oncology clinical trials

<div align="center">
<img src="multiswc_crop.png" alt="multiswc logo" width="180" height="180"/>
</div>

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/multiswc)](https://cran.r-project.org/package=multiswc)
[![CRAN checks](https://badges.cranchecks.info/worst/multiswc.svg)](https://cran.r-project.org/web/checks/check_results_multiswc.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
<!-- badges: end -->

`multiswc` is an R package for analyzing randomized clinical trials with longitudinal follow-up when participants may switch treatment after randomization. The motivating setting is oncology trials with an overall survival endpoint, where control-arm patients may cross over to the experimental therapy, and patients in either arm may later initiate subsequent or next-line therapy.

The package implements a **multi-regime marginal structural Cox model** for the five treatment-regime states

\[
G(t) \in \{C, E, CE, CS, ES\},
\]

where:

- `C`: sustained control;
- `E`: sustained experimental treatment;
- `CE`: control-to-experimental crossover;
- `CS`: control-to-subsequent therapy;
- `ES`: experimental-to-subsequent therapy.

The main estimand is the sustained experimental-versus-control contrast, while the fitted model also keeps clinically distinct post-switch pathways visible rather than collapsing them into one generic switch category.

---

## Installation

Install the released version from CRAN:

```r
install.packages("multiswc")
library(multiswc)
```

You can also install the development version from GitHub:

```r
install.packages("remotes")
remotes::install_github("tonyhbc/multiswc")
library(multiswc)
```

`multiswc` imports `survival` and `nnet`, which are used internally for weighted Cox regression and multinomial regime models.

---

## What problem does `multiswc` solve?

A standard intention-to-treat survival analysis compares randomized arms regardless of what treatment patients actually receive later. That is often the right primary analysis for a treatment policy question, but it can be difficult to interpret when extensive crossover or subsequent therapy occurs after progression.

A simple per-protocol analysis that censors patients at switching is also problematic, because switching is often related to prognosis. Patients may switch because their disease has progressed, toxicity has occurred, or the treating physician believes another therapy is clinically appropriate. Censoring them without adjustment can therefore create selection bias.

`multiswc` addresses this by modeling treatment history as a time-varying categorical regime and using stabilized inverse probability weights to estimate a marginal structural Cox model. The workflow is deliberately compact:

1. Use `simswc()` to simulate trial-like longitudinal survival data with three-way switching.
2. Use `tcoarsen()` when real start-stop data need to be mapped to a regular visit grid.
3. Use `multimsm()` to estimate the weighted multi-regime Cox model.

---

## Package functions at a glance

| Function | Role in workflow | Typical use |
|---|---|---|
| `simswc()` | Simulation engine | Generate longitudinal clinical-trial-like survival data with baseline randomization, time-varying prognostic factors, and one-time switches to `CE`, `CS`, or `ES`. |
| `tcoarsen()` | Preprocessing utility | Convert irregular start-stop follow-up records onto a coarsened visit grid before fitting the MSM. |
| `multimsm()` | Main estimator | Fit a stabilized inverse-probability-weighted multi-regime marginal structural Cox model. |

---

# 1. `simswc()`: simulate trial-like treatment-switching data

`simswc()` creates a synthetic two-arm randomized trial with longitudinal follow-up, time-varying prognostic factors, treatment switching, censoring, and survival outcomes. It is useful for demonstrations, simulation studies, package testing, and teaching the data structure expected by `multimsm()`.

## Important input arguments

| Argument | Meaning | Practical note |
|---|---|---|
| `n` | Number of subjects. | Must be an integer greater than 49 in the current implementation. |
| `n_visit` | Number of planned discrete follow-up visits. | Visit indices are `0, ..., n_visit - 1`. |
| `base_cov` | Baseline covariates included in simulation models and randomization strata. | Choose from `L1`-`L6`; categorical/binary variables such as `c("L1", "L3")` are often convenient. |
| `trt_prob` | Target probability of baseline experimental treatment. | Usually `0.5` for a balanced two-arm trial. |
| `param_tdcov` | Coefficients for the binary time-varying confounder `X`. | Length must be `length(base_cov) + 3`. |
| `param_tdcov2` | Coefficients for the continuous time-varying confounder `U`. | Length must be `length(base_cov) + 3`. |
| `param_sw` | Coefficients for the generic switch model. | Length must be `length(base_cov) + 4`. |
| `param_haz` | Coefficients for the piecewise exponential survival hazard. | Length must be `length(base_cov) + 7`. |
| `param_cens` | Optional coefficients for random censoring. | Use `NULL` for no additional random censoring. |
| `param_select` | Optional model for selecting `CE` versus `CS` among control-arm switchers. | If `NULL`, `p_pick_CE` is used instead. |
| `lagk` | Whether to add lagged treatment and covariate variables. | Use `TRUE` if the downstream IPTW denominator uses lagged history such as `Alag1`. |
| `true_hr` | Whether to compute a Monte Carlo benchmark sustained-regime hazard ratio. | Set to `FALSE` for faster examples. |

## Main output components

| Component | Description |
|---|---|
| `dat` | Subject-level wide-format data with baseline covariates, visit-indexed treatment/covariate/switch histories, observed time, event indicator, censoring, and switch type. |
| `dat_long` | Counting-process long-format data with `t.start`, `t.stop`, `event`, selected covariates, treatment history, and switch indicators. |
| `C` | Matrix for the cumulative censoring process. |
| `stats` | Simulation summary, including initial treatment proportion, switching proportion, switch-type proportions, and event proportion. |
| `true_sdiff` | Optional sustained-regime benchmark if `true_hr = TRUE`; otherwise `NULL`. |

## Tutorial

The following example is the same basic example used in the package documentation. It simulates 200 patients followed over 8 planned visits, with two baseline covariates (`L1`, `L3`), two time-varying prognostic factors (`X`, `U`), and three possible switch types.

```r
set.seed(123)

sim_obj <- simswc(
  n = 200,
  n_visit = 8,
  base_cov = c("L1", "L3"),
  trt_prob = 0.5,
  param_tdcov  = c(-1.5, 0.3, 0.3, -0.2, 0.5),
  param_tdcov2 = c( 0.05, 0.2, 0.2, -0.2, 0.7),
  param_sw     = c(-3.0, 0.4, 0.4, -0.3, 0.3, 0.8),
  param_haz    = c(0.1, log(0.8), log(1.1), log(1.1),
                   log(1.4), log(1.4), log(1.3), log(1), log(0.9)),
  param_cens   = NULL,
  param_select = c(0, 0.5, 0.25),
  lagk = TRUE,
  true_hr = FALSE
)

head(sim_obj$dat_long)
sim_obj$stats
```

This first block lets you inspect the generated long-format dataset and confirm that the simulated trial has plausible treatment assignment, switching, and event frequencies. The long-format dataset is already in start-stop form and can be used directly by `multimsm()` after converting the simulator-specific switch indicators into the main triplet interface.

```r
## Convert to the triplet interface used by multimsm().
dat <- sim_obj$dat_long
dat$rand <- stats::ave(dat$A, dat$id, FUN = function(x) rep(x[1], length(x)))
dat$cross <- dat$S_CE
dat$subseq <- pmax(dat$S_ES, dat$S_CS)

table(dat$rand, dat$cross, dat$subseq)
```

Here, `rand` is the baseline randomized arm carried forward within subject, `cross` identifies control-to-experimental crossover, and `subseq` identifies initiation of any subsequent therapy. The table is a quick structural check: experimental-arm patients should not have control-to-experimental crossover, and each subject should follow the package's one-switch structure.

---

# 2. `multimsm()`: fit the multi-regime marginal structural Cox model

`multimsm()` is the main estimation function. It takes long-format survival data and fits a weighted Cox model in which the treatment exposure is the derived five-level regime process `C`, `E`, `CE`, `CS`, and `ES`.

The required treatment-switching interface is intentionally clinical:

- `rand`: baseline randomized arm (`0` = control, `1` = experimental);
- `cross`: time-varying control-to-experimental crossover indicator;
- `subseq`: time-varying subsequent-therapy indicator.

From these three variables, `multimsm()` constructs the multi-regime state internally.

## Important input arguments

| Argument | Meaning | Practical note |
|---|---|---|
| `dat_long` | Long-format start-stop survival dataset. | One row per subject-interval. |
| `id`, `tstart`, `tstop`, `event` | Names of subject ID, interval start, interval stop, and event indicator. | Defaults match the output of `simswc()`. |
| `cens` | Optional ordinary censoring indicator. | Used only if ordinary censoring weights are requested. |
| `rand` | Baseline randomized arm. | Must be coded as binary `0`/`1`. |
| `cross` | Time-varying control-to-experimental crossover indicator. | Must be absorbing under the one-switch structure. |
| `subseq` | Time-varying subsequent-therapy indicator. | Must be absorbing under the one-switch structure. |
| `base_cov` | Baseline covariates included in the final Cox model. | Optional but often used for precision and design consistency. |
| `iptw_num` | RHS-only formula for stabilized numerator regime model. | Example: `~ regime_lag + factor(visit) + L1 + L3`. |
| `iptw_den` | RHS-only formula for denominator regime model. | Required; should include confounding history. |
| `ipcw_mod` | Optional RHS-only formula for ordinary censoring model. | Used when `cens` is supplied. |
| `wt_trunc` | Optional upper quantile for weight truncation. | Example: `0.95`. |
| `prob_bounds` | Lower and upper bounds for predicted probabilities. | Helps protect against numerical instability. |
| `normalize_weights` | Whether to normalize final weights to mean 1. | Defaults to `TRUE`. |
| `robust` | Whether to use robust cluster variance by subject. | Defaults to `TRUE`. |

## Main output components

| Component | Description |
|---|---|
| `coef_table` | Full coefficient table for the regime effects. |
| `coef` | Estimated log-hazard ratios for treatment regimes against sustained control. |
| `hr` | Hazard ratios for treatment regimes against sustained control. |
| `hr_ci` | Wald 95% confidence intervals for hazard ratios. |
| `fit` | Fitted weighted Cox model from `survival::coxph()`. |
| `p.val` | P-values for the regime effect estimates. |
| `dat_long` | Augmented analysis dataset with derived regime variables and weights. |
| `regime_models` | Fitted numerator and denominator multinomial models. |
| `censor_model` | Fitted ordinary censoring model if `ipcw_mod` is supplied; otherwise `NULL`. |
| `diagnostics` | Weight quantiles, final regime counts/proportions, cohort size, and truncation information. |
| `call` | Matched function call. |

## Tutorial

The example below simulates a trial-like dataset, converts it to the triplet interface, and fits the multi-regime MSM. This mirrors the roxygen example for `multimsm()`.

```r
set.seed(123)

sim_obj <- simswc(
  n = 300,
  n_visit = 8,
  base_cov = c("L1", "L3"),
  trt_prob = 0.5,
  param_tdcov = c(-1.5, 0.3, 0.3, -0.2, 0.5),
  param_tdcov2 = c(0.05, 0.2, 0.2, -0.2, 0.7),
  param_sw = c(-3, 0.4, 0.4, -0.3, 0.3, 0.8),
  param_haz = c(0.1, log(0.8), log(1.1), log(1.1),
                log(1.4), log(1.4), log(1.3), log(1), log(0.9)),
  param_cens = NULL,
  param_select = c(0, 0.5, 0.25),
  lagk = TRUE,
  true_hr = FALSE
)

# Convert to the required treatment switch indicator triplet
dat <- sim_obj$dat_long
dat$rand <- stats::ave(dat$A, dat$id, FUN = function(x) rep(x[1], length(x)))
dat$cross <- dat$S_CE
dat$subseq <- pmax(dat$S_ES, dat$S_CS)

fit <- multimsm(
  dat_long = dat, id = "id", tstart = "t.start", tstop = "t.stop",
  event = "event", rand = "rand", cross = "cross", subseq = "subseq",
  base_cov = c("L1", "L3"),
  iptw_num = ~ regime_lag + factor(visit) + L1 + L3,
  iptw_den = ~ regime_lag + factor(visit) + L1 + L3 + X + U + Alag1,
  wt_trunc = 0.95
)

fit
```

The numerator model defines the stabilizing part of the treatment-regime weight. The denominator model includes the richer confounding history used to estimate how likely each observed regime was, given the patient's past and current information. In this example, the denominator includes time-varying prognostic factors `X` and `U` and lagged treatment `Alag1`.

Printing the fitted object gives a compact report with:

- the number and percentage of patients ending in each switch pattern;
- log-hazard ratios and hazard ratios by regime;
- 95% confidence intervals and p-values;
- selected weight quantiles.

The coefficient for `regime.E` is the primary sustained experimental-versus-control contrast. The coefficients for `regime.CE`, `regime.CS`, and `regime.ES` summarize survival contrasts for distinct switched-regime person-time relative to sustained control.

A typical printed object looks like this:

```r
==========================================
Multi-regime marginal structural Cox model
==========================================
Switch summary (N = 300):
   *  # of always C: ...
   *  # of always E: ...
   *  # of CE switch: ...
   *  # of CS switch: ...
   *  # of ES switch: ...
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Hazard ratios (HR) by regime
           logHR    HR lower upper p.val
regime.E    ...    ...  ...   ...  ...
regime.CE   ...    ...  ...   ...  ...
regime.CS   ...    ...  ...   ...  ...
regime.ES   ...    ...  ...   ...  ...
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Weight quantiles:
          0%    5%   50%   95%  100%
S.IPTW   ...   ...   ...   ...   ...
------------------------------------------
A `multimsm` object
```

For downstream checks, inspect the returned diagnostics:

```r
fit$diagnostics
head(fit$dat_long)
summary(fit$dat_long$.w_use_trunc)
```

These are useful for checking regime support, positivity, weight behavior, and whether truncation materially affects the analysis.

---

# 3. `tcoarsen()`: coarsen irregular start-stop data onto a visit grid

Real clinical datasets rarely arrive on the clean discrete visit grid used in simulations. Patients may have visits at irregular calendar times, and covariates may be updated at different times for different patients. `tcoarsen()` helps convert such data into a regularized start-stop structure suitable for longitudinal weighting workflows.

The function snaps update times to a user-specified grid, rebuilds start-stop intervals, and carries forward the last observed values of selected variables. It can also force binary treatment indicators to be absorbing, which is helpful for variables such as crossover or subsequent therapy initiation.

## Important input arguments

| Argument | Meaning | Practical note |
|---|---|---|
| `data` | Start-stop long-format dataset. | One or more rows per subject. |
| `id`, `start`, `stop`, `event` | Names of subject ID, interval start, interval stop, and event indicator. | Required. |
| `covs` | Variables to document and check as carried-forward predictors. | Include treatment and time-varying covariates. |
| `bin_width` | Width of the coarsening grid in the same time unit as `start` and `stop`. | Example: `30` for 30-day bins. |
| `dir_coarsen` | Direction for snapping update times. | `"floor"` moves updates backward to the start of the bin; `"ceiling"` moves them forward. |
| `origin` | Grid origin. | Often `0` in trial time. |
| `absorb_vars` | Binary variables to convert to absorbing status by cumulative maximum. | For switching analyses, often `c("cross", "subseq")`. |
| `keep_terminal_time` | Whether to preserve exact event/censoring exit time. | Defaults to `TRUE`, usually preferred. |
| `gap_action` | How to handle gaps before coarsening. | One of `"stop"`, `"warn"`, or `"ignore"`. |
| `add_visit`, `visit_name` | Whether to add a discrete visit index and what to call it. | Defaults to `TRUE` and `"visit"`. |
| `diagnostics` | Whether to return preprocessing diagnostics. | Defaults to `TRUE`. |

## Main output components

| Component | Description |
|---|---|
| `dat_coarsen` | Time-coarsened long-format dataset with rebuilt start-stop intervals and carried-forward covariate values. |
| `diagnostics` | Settings and basic summaries describing the coarsening operation. |

## Tutorial

The roxygen example uses the `heart` dataset from the `survival` package. This is a convenient public start-stop survival dataset with a time-varying transplant indicator.

```r
if (requireNamespace("survival", quietly = TRUE)) {
  data("heart", package = "survival")

  heart_30 <- tcoarsen(
    data = heart,
    id = "id",
    start = "start",
    stop = "stop",
    event = "event",
    covs = c("transplant", "age", "surgery"),
    bin_width = 30,
    dir_coarsen = "floor",
    origin = 0,
    absorb_vars = "transplant",
    verbose = FALSE
  )

  utils::head(heart_30$dat_coarsen, 10)
}
```

This example maps follow-up time to a 30-day grid using floor coarsening. The transplant indicator is treated as absorbing, meaning once a patient is recorded as transplanted, the status remains active in later carried-forward rows. The same logic is useful for treatment-switching variables such as `cross` and `subseq`, because crossover and subsequent therapy initiation are usually interpreted as absorbing treatment-history events.

For a treatment-switching dataset, the call would often look like:

```r
dat_grid <- tcoarsen(
  data = dat_raw,
  id = "id",
  start = "start",
  stop = "stop",
  event = "event",
  covs = c("rand", "cross", "subseq", "age", "progression"),
  bin_width = 30,
  dir_coarsen = "floor",
  origin = 0,
  absorb_vars = c("cross", "subseq")
)

analysis_dat <- dat_grid$dat_coarsen
```

After coarsening, `analysis_dat` can be passed to `multimsm()` if it contains the required survival variables, triplet treatment-switching indicators, and covariate history needed for the numerator and denominator weight models.

---

## Recommended analysis workflow

A typical workflow for a real or simulated analysis is:

```r
library(multiswc)

# 1. Start with long-format data or generate trial-like data.
sim_obj <- simswc(
  n = 300,
  n_visit = 8,
  base_cov = c("L1", "L3"),
  trt_prob = 0.5,
  param_tdcov = c(-1.5, 0.3, 0.3, -0.2, 0.5),
  param_tdcov2 = c(0.05, 0.2, 0.2, -0.2, 0.7),
  param_sw = c(-3, 0.4, 0.4, -0.3, 0.3, 0.8),
  param_haz = c(0.1, log(0.8), log(1.1), log(1.1),
                log(1.4), log(1.4), log(1.3), log(1), log(0.9)),
  param_cens = NULL,
  param_select = c(0, 0.5, 0.25),
  lagk = TRUE,
  true_hr = FALSE
)

# 2. Create the treatment-switching triplet.
dat <- sim_obj$dat_long
dat$rand <- stats::ave(dat$A, dat$id, FUN = function(x) rep(x[1], length(x)))
dat$cross <- dat$S_CE
dat$subseq <- pmax(dat$S_ES, dat$S_CS)

# 3. Fit the multi-regime MSM.
fit <- multimsm(
  dat_long = dat,
  id = "id",
  tstart = "t.start",
  tstop = "t.stop",
  event = "event",
  rand = "rand",
  cross = "cross",
  subseq = "subseq",
  base_cov = c("L1", "L3"),
  iptw_num = ~ regime_lag + factor(visit) + L1 + L3,
  iptw_den = ~ regime_lag + factor(visit) + L1 + L3 + X + U + Alag1,
  wt_trunc = 0.95
)

fit
```

For an applied analysis, spend time checking:

- whether the treatment-history variables truly obey the one-switch assumption;
- whether each regime has enough support over follow-up;
- whether predicted regime probabilities are near zero for some histories;
- whether weight truncation changes the estimate materially;
- whether alternative numerator/denominator histories lead to similar conclusions.

---

## Interpretation of the fitted regimes

The reference regime is sustained control (`C`). The estimated hazard ratios are interpreted relative to `C`:

| Coefficient | Regime contrast | Interpretation |
|---|---|---|
| `regime.E` | `E` versus `C` | Primary sustained experimental-versus-control contrast. |
| `regime.CE` | `CE` versus `C` | Survival contrast for control-to-experimental crossover person-time relative to sustained control. |
| `regime.CS` | `CS` versus `C` | Survival contrast for control-to-subsequent therapy person-time relative to sustained control. |
| `regime.ES` | `ES` versus `C` | Survival contrast for experimental-to-subsequent therapy person-time relative to sustained control. |

The `regime.E` coefficient is usually the main causal contrast of interest. The post-switch coefficients are useful because they keep clinically distinct switching pathways visible, but their interpretation depends on the usual marginal structural model assumptions, the treatment-history coding, positivity, correct model specification, and the clinical meaning of the switch states.

---

## Assumptions and current scope

`multiswc` is intentionally focused. The current implementation assumes:

- a two-arm baseline randomized trial or trial-like study;
- start-stop survival follow-up;
- at most one treatment switch per subject;
- possible pathways `C`, `E`, `CE`, `CS`, and `ES`;
- no control-to-experimental-to-subsequent sequential path in the same subject;
- correctly measured switching history and relevant prognostic history;
- sufficient support for the observed treatment regimes;
- analyst-specified numerator and denominator models for stabilized weights.

If the clinical setting allows multiple sequential switches, a richer regime state space is required. In that situation, the current five-state model should not be interpreted as capturing the full treatment process.

---

## Reporting suggestions

When reporting an analysis with `multiswc`, we recommend including:

1. the clinical definition of crossover and subsequent therapy;
2. the visit grid or coarsening rule;
3. the numerator and denominator weight-model formulas;
4. the final counts in `C`, `E`, `CE`, `CS`, and `ES`;
5. weight summaries before and after truncation;
6. the hazard ratio for `regime.E` as the primary contrast;
7. sensitivity analyses for weight truncation and model history specification.

---

## Citation

To cite `multiswc`, use:

```r
citation("multiswc")
```

---

## Bug reports and feature requests

Please report bugs or request features at:

<https://github.com/tonyhbc/multiswc/issues>

---

## License

`multiswc` is released under the MIT license.
