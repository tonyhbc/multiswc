---
title: "multiswc: An R package for multi-way treatment switch analysis in oncology
  clinical trials"
tags:
- R
- biostatistics
- clinical trials
- causal inference
- survival analysis
- treatment switching
date: "31 May 2026"
output:
  pdf_document: default
  html_document:
    df_print: paged
authors:
- name: Haobin Chen
  affiliation: 1
  corresponding: true
- name: Yuxuan Chen
  affiliation: 2
- name: Philip He
  affiliation: 2
bibliography: paper.bib
affiliations:
- name: Program in Quantitative Biomedical Sciences, The Geisel School of Medicine at Dartmouth,
    NH, USA
  index: 1
- name: Biostatistics and Data Management, Daiichi Sankyo Inc., NJ, USA
  index: 2
---

# Summary

`multiswc` is an R package for analyzing randomized clinical trials with longitudinal follow-up in which participants may change treatment after randomization. The direct motivating setting is oncology phase-3 clinical trials with an overall survival endpoint, where patients assigned to control may later receive the experimental therapy, and patients in either arm may initiate subsequent or next-in-line therapy after progression or other clinical status updates. These post-randomization treatment changes are clinically important by altering patient's survival trajectory, but they can make a standard intention-to-treat survival analysis difficult to interpret when the survival effect contrast of interest is coomparing what would have been observed under sustained initial treatment strategies without treatment switch.

The package provides three user-interface functions. `simswc()` faithfully simulate a longitudinal survival data from a oncology phase-3 trial with baseline randomization, time-varying prognostic factors, and three possible one-time switching pathways: control to experimental treatment, control to subsequent therapy, and experimental treatment to subsequent therapy. `multimsm()` fits a multi-regime marginal structural Cox model using stabilized inverse probability weights from multinomial regime models. `tcoarsen()` is a preprocessing utility for mapping irregular start-stop follow-up records onto a coarsened, prespecified longitudinal visit grid before fitting the model. Together, these functions support simulation, data preparation, model fitting, weight diagnostics, and interpretable treatment-regime effects in a single workflow.

# Statement of need

Treatment switching is common in randomized trials with time-to-event outcomes, especially in oncology, where crossover or access to subsequent therapy may be ethically or clinically appropriate after disease progression [@Fda2025; @Latimer2016]. If patients switch from their assigned treatment in ways related to prognosis, simple analyses will answer a different question from the one intended by clinicians, regulatory authority, or health technology assessment groups. Intention-to-treat analysis compares only treatment randomization and disregards the observed switching process, while per-protocol censoring at switch can introduce selection bias [@Robins2000]; @Morden2011; @Latimer2017]. Several methods exist for the analysis of longitudinal survival data subject to treatment crossover to the alternative arms, including inverse probability of censoring weights (IPCW), marginal structural models (MSM), rank-preserving structural failure time models (RPSFTM), and two-stage estimation (TSE) [@Robins2000; @Robins1991; @Robins20002; @Latimer2017].

A practical gap remains for contemporary trials where switching is not a single binary event, which has received less attention until recenly [@Gorrod2025; @Huang2026]. A control participant may cross over to the experimental therapy (*CE*), switch to another subsequent therapy (*CS*), or remain on control (*C*); an experimental participant may remain on experimental therapy (*E*) or switch to subsequent therapy (*CS*). Statistical methods, as well as software implementation, for time-to-event data that handles such three-way switching have been limited. An example is the KEYNOTE-189 trials comparing pembrolizumab (pembro) against chemotherapy, where chemotherapy patients may crossover to pembrolizumab if symptoms worsen and all patients may opt for subsequent therapy at any time [@Gandhi2018]. `multiswc` addresses this gap by implementing a multi-regime marginal structural Cox model representing the observed post-randomization treatment status as a time-varying categorical regime, $G(t) \in \{C,E,CE,CS,ES\}$, while preserving the primary sustained experimental-versus-control causal contrast of interest for pharmaceutical developers.

# State of the field

Several R packages implement methods for treatment-switching analysis. `rpsftm` implements rank-preserving structural failure time models for two-arm trials with between-arm crossovers [@Allison2017]. `ipcwswitch` carries out artificial censoring on post-switching survival times and proceeds with survival analysis on censored data under a hypothetical setting in which patients remain on their randomized treatment [@Graffeo2019]. `trtswitch` contains a broad collection of established adjustment methods, including RPSFTM, iterative parameter estimation, IPCW, TSE, and a MSM that allows for only single switching destination [@Lu2024]. These packages are valuable for standard treatment-switching analysis, but they are not specifically designed to simultaneously distinguish control-to-experimental crossover from control-to-subsequent and experimental-to-subsequent therapy switching pathways.

`multiswc` was built as a focused extension for existing tools. Its contribution is the multi-regime parameterization for three-way switching, improved upon existing MSM that only allows for a single switching pathway [@Yamaguchi2004; @Xu2022], under standard Cox analysis framework. The design resembles the broader marginal structural model literature for multiple treatment categories [@Suarez2008], but tailored to the randomized trial setting where crossover and subsequent therapy have different timing and clinical meanings. This makes the package most useful for pharmaceutical statisticians, pharmacoepidemiologists, and health technology assessment analysts who need to report both a sustained no-switch treatment effect and descriptive effect contrasts for distinct post-switch regimes.

# Software design

The central design of `multiswc` allows users supply treatment-switching history in clinically recognizable variables typically collected in a randomized trial with treatment switching. When specifying the `multimsm()` analysis, one uses `rand` for the baseline randomized arm, `cross` for the time-varying control-to-experimental crossover indicator, and `subseq` for the time-varying subsequent-therapy indicator. The function then internally constructs the five-level regime and checks that the data satisfy the package's assumed one-switch structure. Under this structure, control-arm participants may remain in `C`, transition to `CE`, or transition to `CS`; experimental-arm participants may remain in `E` or transition to `ES`. Paths with multiple sequential switches, such as control to experimental followed by subsequent therapy, require an expanded regime set and are not considered in the current implementation.

The fitted marginal structural Cox model with regime status $G(t)$ is

\begin{equation}\label{eq:msm}
\lambda_{T^{\bar g}}(t) = \lambda_0(t)\exp\{\beta_e \mathbb{I}(G(t)=E) + \beta_{ce} \mathbb{I}(G(t)=CE) + \beta_{cs} \mathbb{I}(G(t)=CS) + \beta_{es} \mathbb{I}(G(t)=ES)\}
\end{equation}

with $\mathbb{I}(\cdot)$ as the indicator function and `C` as the reference regime. The coefficient for `E`, $\beta_e$, is the primary sustained experimental-versus-control causal contrast in log-hazard ratio; the remaining coefficients summarize survival contrasts for the observed switched-regime person-time relative to sustained control.

The stabilized regime weight at visit interval $k$ is

\begin{equation}\label{eq:weights}
W_i(k)=\prod_{j=1}^{k}\frac{\Pr\{G_{ij}=g_{ij}\mid \mathcal S_{ij}\}}{\Pr\{G_{ij}=g_{ij}\mid \mathcal H_{ij}\}}
\end{equation}

where the numerator model uses a reduced stabilizing history and the denominator model includes the analyst-specified confounding history, typically including prior regime, visit time, baseline covariates, and time-varying prognostic factors. `multimsm()` estimates these probabilities with `nnet::multinom()` and fits the weighted Cox model with `survival::coxph()` [@Ripley2009; @Therneau2001]. Optional ordinary censoring weights may be multiplied into the final weight. The function also provides probability bounds, quantile truncation, weight normalization, robust cluster variance by participant, and printed diagnostics for regime support and weight distributions.

| Function | Arguments | Purpose |
|------------------------|------------------------|------------------------|
| `simswc()` | `n`, `n_visit`, `base_cov`, `param_tdcov`, `param_tdcov2`, `param_sw`, `param_haz` | Simulate trial-like longitudinal survival data with time-varying confounders, switching, and piecewise exponential hazards. |
| `multimsm()` | `id`, `tstart`, `tstop`, `event` | Identify the counting-process survival outcome structure. |
| `multimsm()` | `rand`, `cross`, `subseq` | Define the randomized arm, crossover process, and subsequent-therapy process used to derive `C`, `E`, `CE`, `CS`, and `ES`. |
| `multimsm()` | `iptw_num`, `iptw_den`, `base_cov` | Specify stabilized numerator and denominator multinomial regime models. |
| `multimsm()` | `ipcw_mod`, `cens`, `wt_trunc`, `prob_bounds` | Add ordinary censoring weights and control extreme estimated weights. |
| `tcoarsen()` | `bin_width`, `dir_coarsen`, `absorb_vars` | Harmonize irregular start-stop data to a visit grid and carry absorbing treatment indicators forward. |

: Selected user-facing arguments in the main `multiswc` workflow. \label{tab:args}

A minimal analysis begins by simulating data with `simswc`, and the true marginal causal hazard ratio can be obtained from `sim$true_sdiff$cHR` by setting `true_hr = TRUE`.

``` r
> library(multiswc)

> set.seed(2026)
> sim <- simswc(
+   n = 500, n_visit = 8, base_cov = c("L1", "L3"), trt_prob = 0.5,
+   param_tdcov = c(-1.5, 0.3, 0.3, -0.35, 0.5),
+   param_tdcov2 = c( 0.05, 0.2, 0.2, -0.35, 0.7),
+   param_sw = c(-1.5, 0.4, 0.4, -0.3, 0.3, 0.8),
+   param_haz = c(0.1, log(0.7), log(1.1), log(1.1),
+                    log(1.5), log(1.5), log(1.3), log(1), log(0.9)),
+   param_cens = NULL, param_select = c(0, 0.5, 0.25),
+   lagk = TRUE, true_hr = TRUE
+ )

> sim$true_sdiff$cHR
[1] 0.5452675
```

Then we can convert the simulator's switch-type indicators into the triplet interface (i.e., `rand`, `cross`, `subseq`) expected by `multimsm()`, and fitting the IPTW-weighted marginal Cox model:

``` r
> dat <- sim$dat_long
> dat$rand <- ave(dat$A, dat$id, FUN = function(x) rep(x[1], length(x)))
> dat$cross <- dat$S_CE
> dat$subseq <- pmax(dat$S_ES, dat$S_CS)

> mm_fit <- multimsm(
+  dat_long = dat,
+  id = "id", tstart = "t.start", tstop = "t.stop", event = "event",
+  rand = "rand", cross = "cross", subseq = "subseq",
+  base_cov = c("L1", "L3"),
+  iptw_num = ~ regime_lag + factor(visit) + L1 + L3,
+  iptw_den = ~ regime_lag + factor(visit) + L1 + L3 + X + U + Alag1,
+  wt_trunc = 0.95
+ )

> mm_fit
==========================================
Multi-regime marginal structural Cox model
==========================================
Switch summary (N = 500):
   *  # of always C: 191 (38.2%)
   *  # of always E: 201 (40.2%)
   *  # of CE switch: 37 (7.4%)
   *  # of CS switch: 22 (4.4%)
   *  # of ES switch: 49 (9.8%)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Hazard ratios (HR) by regime
           logHR    HR lower upper p.val    
regime.E  -0.638 0.528 0.414 0.674 0.000 ***
regime.CE -0.777 0.460 0.273 0.774 0.003  **
regime.CS  0.169 1.185 0.733 1.915 0.489    
regime.ES  0.437 1.549 0.995 2.410 0.052    
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Weight quantiles:
          0%    5%   50%   95%  100%
S.IPTW 0.597 0.597 1.013 1.408 1.408
------------------------------------------
A `multimsm` object
```

The printed `multimsm` object reports log-hazard ratios, hazard ratios, 95% confidence intervals, p-values, final treatment regime counts at the end of follow-up, and selected weight quantiles. The coefficients for `regime.E` is the primary causal effect of experimental treatment compared to control should there be no switching. The returned object retains the fitted Cox model, fitted numerator and denominator regime models, augmented long-format data with derived weights, and diagnostics for downstream sensitivity analyses.

# Research impact statement

`multiswc` is available on CRAN as version 0.1.2 under an MIT license and can be installed through the standard R package ecosystem [@Chen2026]. Its public repository and reference manual expose the full implementation of simulation, preprocessing, and estimation functions. The package provides easy implementation of multi-regime MSM in complex treatment-switching setting that arises in oncology trials and related comparative-effectiveness studies.

The package is also designed to support reproducible methods research. `simswc()` enables simulation studies in which the analyst can vary switching frequency, switch destination, time-varying prognostic factors, censoring, and survival effects. `tcoarsen()` addresses a common obstacle in real clinical datasets: longitudinal updates are often measured at irregular times, whereas many MSM workflows require an analysis grid. These features make `multiswc` useful not only for a single applied analysis but also for sensitivity analyses, tutorials, simulation benchmarks, and comparison with established approaches such as IPCW and two-stage methods. In practical terms, the package provides a reusable implementation for a clinically important scenario where methods exist in the literature but software support remains limited.

# AI usage disclosure

ChatGPT was used to assist with the editing this JOSS paper. The scientific framing, software design choices, examples, and final manuscript content were completed by the human authors, who remain responsible for the accuracy, originality, and compliance of all submitted materials.

# Acknowledgements

The authors thank the developers and maintainers of the R ecosystem packages on which `multiswc` depends, especially `survival` and `nnet`.

# References

```{=html}
<!--
Suggested BibTeX keys used in this draft: Allison2017, CainCole2009,
Graffeo2019, Keogh2023, Latimer2018, Lu2026trtswitch, Morden2011,
multiswcCRAN, Robins2000, RobinsFinkelstein2000, Suarez2008,
Therneau2024, VenablesRipley2002, Ying2023.
-->
```
