# `multiswc`: Instrumental variable methods for survival analysis

<img src="multiswc_crop.png" alt="surviv" width="180" height="180"/>

Instrumental variable methods for **time-to-event outcomes** under **unmeasured confounding**, with a focus on IV-based Cox model estimators under various real-world data setting.

---

## Installation

Install the offical version of the package from CRAN and load it:

```r
install.packages("surviv")
library(surviv)
```

---

## Analytic goal at a glance

In many survival follow-up studies, treatment is confounded by factors that are not fully measured. **Instrumental variables (IVs)** can enable causal estimation and mitigate bias due to unmeasured confounding when a valid instrument exists (e.g., randomization, provider preference, site-level variation, policy/eligibility enactment), provided standard IV assumptions hold.

This package collects several **IV estimators for Cox-type survival models**, spanning classic baseline IV estimators and time-varying designs for flexible real-world analytics.
