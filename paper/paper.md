---
title: "multiswc: An R package for analyzing survival data with multi-way treatment switching in oncology clinical trials"
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
- name: Program in Quantitative Biomedical Sciences, Geisel School of Medicine at Dartmouth College,
    NH, USA
  index: 1
- name: Biostatistics and Data Management, Daiichi Sankyo Inc., NJ, USA
  index: 2
---

# Summary

`multiswc` is an R package for analyzing randomized clinical trials with longitudinal follow-up in which participants may change treatment after randomization. The motivating setting is oncology phase-3 trials with overall survival endpoints, where control-arm patients may later receive the experimental therapy and patients in either arm may initiate subsequent therapy after progression or other clinical updates. These post-randomization changes can alter survival trajectories and make a standard intention-to-treat survival analysis hard to interpret when the target contrast concerns sustained initial treatment strategies without switching.

The package provides three user-facing functions. `simswc()` simulates longitudinal survival data from an oncology-style randomized trial with baseline randomization, time-varying prognostic factors, and three possible one-time switching pathways: control to experimental treatment, control to subsequent therapy, and experimental treatment to subsequent therapy. `multimsm()` fits a multi-regime marginal structural Cox model using stabilized inverse probability weights from multinomial regime models. `tcoarsen()` maps irregular start-stop follow-up records onto a coarsened visit grid before model fitting. Together, these functions support simulation, data preparation, model fitting, weight diagnostics, and interpretable treatment-regime effects in a single workflow.

# Statement of need

Treatment switching is common in randomized trials with time-to-event outcomes, especially in oncology, where crossover or access to subsequent therapy may be ethically or clinically appropriate after disease progression [@Fda2025; @Latimer2016]. If switching is related to prognosis, simple analyses may answer a different question from the one intended by clinicians, regulators, or health technology assessment groups. Intention-to-treat analysis compares randomized groups while disregarding the observed switching process, whereas per-protocol censoring at switch can introduce selection bias [@Robins2000; @Morden2011; @Latimer2017]. Established adjustment approaches include inverse probability of censoring weights (IPCW), marginal structural models (MSM), rank-preserving structural failure time models (RPSFTM), and two-stage estimation (TSE) [@Robins2000; @Robins1991; @Robins20002; @Latimer2017].

A practical gap remains for contemporary trials in which switching is not a single binary event, a setting that has received growing attention only recently [@Gorrod2025; @Huang2026]. A control participant may remain on control (*C*), cross over to experimental therapy (*CE*), or switch to subsequent therapy (*CS*); an experimental participant may remain on experimental therapy (*E*) or switch to subsequent therapy (*ES*). Software for time-to-event analyses that explicitly represents these three switching pathways remains limited. For example, in KEYNOTE-189, a trial comparing pembrolizumab plus chemotherapy with placebo plus chemotherapy, patients assigned to placebo could cross over to pembrolizumab after verified progression, while subsequent anticancer therapies could also occur [@Gandhi2018]. `multiswc` addresses this gap by representing post-randomization treatment status as a time-varying categorical regime, $G(t) \in \{C,E,CE,CS,ES\}$, while preserving the primary sustained experimental-versus-control contrast.

# State of the field

Several R packages implement treatment-switching methods. `rpsftm` implements rank-preserving structural failure time models for two-arm trials with between-arm crossover [@Allison2017]. `ipcwswitch` implements artificial censoring at switch followed by IPCW survival analysis under a hypothetical setting in which patients remain on randomized treatment [@Graffeo2019]. `trtswitch` includes RPSFTM, iterative parameter estimation, IPCW, TSE, and an MSM for a single switching destination [@Lu2024]. These packages are valuable for standard switching analyses, but they are not designed to distinguish control-to-experimental crossover from control-to-subsequent and experimental-to-subsequent therapy pathways in one survival model.

`multiswc` was built as a focused extension rather than a replacement for these tools. Its main contribution is a multi-regime parameterization for three-way switching, extending existing MSM formulations that use baseline randomization and a single post-randomization treatment process [@Yamaguchi2004; @Xu2022]. The design is related to marginal structural models for multiple treatment categories [@Suarez2008], but is tailored to randomized oncology trials in which crossover and subsequent therapy have different timing, clinical meanings, and interpretive roles. The target users are pharmaceutical statisticians, pharmacoepidemiologists, and health technology assessment analysts who need to report both a sustained no-switch treatment effect and contrasts for distinct post-switch regimes.

# Software design

The design of `multiswc` lets users supply switching histories through variables commonly collected in randomized trials. In `multimsm()`, `rand` identifies baseline randomized arm, `cross` identifies time-varying control-to-experimental crossover, and `subseq` identifies time-varying subsequent-therapy initiation. The function internally constructs the five-level regime and checks the assumed one-switch structure: control-arm participants may remain in `C`, transition to `CE`, or transition to `CS`, whereas experimental-arm participants may remain in `E` or transition to `ES`. Multiple sequential switches, such as control to experimental followed by subsequent therapy, require an expanded regime set and are not covered by the current implementation.

The fitted marginal structural Cox model with regime status $G(t)$ is

\begin{align*}\label{eq:msm}
\lambda_{T^{\bar g}}(t) = \lambda_0(t)\exp\{&\beta_E \mathbb{I}(G(t)=E) + \beta_{CE} \mathbb{I}(G(t)=CE) + \\ 
&\beta_{CS} \mathbb{I}(G(t)=CS) + \beta_{ES} \mathbb{I}(G(t)=ES)\}
\end{align*}

Here `C` is the reference regime. The coefficient for `E`, $\beta_E$, is the primary sustained experimental-versus-control log hazard ratio; the remaining coefficients summarize survival contrasts for switched-regime person-time relative to sustained control.

The stabilized regime weight at visit interval $k$ is

\begin{equation}\label{eq:weights}
W_i(k)=\prod_{j=1}^{k}\frac{\Pr\{G_{ij}=g_{ij}\mid \mathcal S_{ij}\}}{\Pr\{G_{ij}=g_{ij}\mid \mathcal H_{ij}\}}
\end{equation}

where the numerator model uses a reduced stabilizing history and the denominator model includes the analyst-specified confounding history, typically prior regime, visit time, baseline covariates, and time-varying prognostic factors. `multimsm()` estimates these probabilities with `nnet::multinom()` and fits the weighted Cox model with `survival::coxph()` [@Ripley2009; @Therneau2001]. Optional ordinary censoring weights may be multiplied into the final weight. The function also supports probability bounds, quantile truncation, weight normalization, robust cluster variance by participant, and printed diagnostics for regime support and weight distributions.

A minimal analysis begins by simulating data with `simswc`; the true marginal causal hazard ratio can be obtained from `sim$true_sdiff$cHR` when `true_hr = TRUE`.

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

The simulator's switch-type indicators are then converted into the triplet interface expected by `multimsm()`, followed by the IPTW-weighted marginal Cox model.

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

The printed object reports log hazard ratios, hazard ratios, 95% confidence intervals, p-values, final treatment-regime counts, and selected weight quantiles. The coefficient for `regime.E` is the primary no-switch experimental-versus-control effect. The returned object retains the fitted Cox model, fitted numerator and denominator regime models, augmented long-format data with derived weights, and diagnostics for sensitivity analyses.

# Research impact statement

`multiswc` is available on CRAN as version 0.1.2 under an MIT license and can be installed through the standard R package ecosystem [@Chen2026]. The package provides an accessible implementation of multi-regime MSMs for complex treatment-switching settings arising in oncology trials and related comparative effectiveness studies.

The package is also designed for reproducible methods research. `simswc()` enables simulation studies varying switching frequency, switch destination, time-varying prognostic factors, censoring, and survival effects. `tcoarsen()` performs follow-up time coarsening to address a common obstacle in real clinical datasets that longitudinal updates are often measured at irregular times, whereas many MSM workflows require an analysis grid. These features make `multiswc` useful for applied analyses, tutorials, simulation benchmarks, sensitivity analyses, and comparisons with established approaches such as IPCW and two-stage methods on complex treatment-switching analysis. Practically, the package provides an easy implementation for a clinically important scenario where both specific methods and software implementations remain limited.

# AI usage disclosure

ChatGPT was used to assist with editing this JOSS paper. The scientific framing, software design choices, examples, and final manuscript content were completed and verified by the human authors, who remain responsible for the accuracy, originality, and compliance of all submitted materials.

# Acknowledgements

The authors thank the developers and maintainers of the R ecosystem packages on which `multiswc` depends, especially `survival` and `nnet`.

# References
