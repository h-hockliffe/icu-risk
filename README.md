[![CI](https://img.shields.io/github/actions/workflow/status/h-hockliffe/icu-risk/ci.yaml?branch=main)](https://github.com/h-hockliffe/icu-risk/actions)
[![Reproducible](https://img.shields.io/badge/reproducible-renv-blue)](#reproducibility)
# ICU Mortality — Bayesian Cause-Specific Baseline Model (INLA) with Cox & Fine–Gray Comparators


A reproducible workflow for developing and validating **competing-risk ICU mortality models** using:

-   **Bayesian cause-specific hazards** with **INLA** (primary),
-   **Cox frailty** comparators,
-   **Fine–Gray** subdistribution models (used for descriptive comparisons).

Includes a `targets` pipeline, Quarto report, GitHub Actions CI, and a test suite.

------------------------------------------------------------------------

## What’s new (Aug 2025)

-   **Creatinine vs eGFR collinearity**: Primary models use **log-eGFR only**; creatinine enters only via a **residualised** sensitivity term (orthogonal to eGFR, age, sex).
-   **Shared age centering**: A shared `age_mean` per dataset ensures aligned transforms across engines.
-   **Interval selection**: early-stops when cutpoints stabilise and uses a bounded iteration cap.
-   **Index-stay dataset**: **No administrative right-censoring** (earliest ICU stay per patient).
-   **Fine–Gray usage**: FG reported **descriptively**; primary inference is **cause-specific**.
-   **Time-dependent AUC** and **Brier score** at specified horizons, **averaged across imputations** (bootstrap CIs **optional** when enabled).

## Directory layout

``` text
├── Thesis.qmd                 # Main Quarto report
├── analysis.qmd               # (older) Quarto doc – see Thesis.qmd
├── _targets.R                 # targets pipeline
├── R/
│   ├── 00_config.R
│   ├── 00_functions_collapse.R
│   ├── 00_functions_guardrail.R
│   ├── 01_clean.R
│   ├── functions_data.R
│   ├── functions_impute.R
│   ├── functions_transform.R
│   ├── functions_models.R
│   ├── functions_validate.R
│   ├── functions_plots.R
│   ├── functions_cif_helpers.R
│   ├── functions_metrics.R
│   ├── functions_tables.R
│   ├── functions_mcmc_check.R
│   ├── functions_surv_helpers.R
│   └── functions_sensitivity.R
├── tests/
│   ├── testthat.R
│   └── testthat/...
├── renv.lock
└── README.md
```

------------------------------------------------------------------------

## Data & privacy

-   The raw workbook lives on **Google Drive** and is **not** committed to this repo.
-   Point the pipeline to the Drive-mounted path via the `RAW_FILE` environment variable (see below). The file is read **in place** from Drive.
-   The pipeline reads, cleans, and writes a **pseudonymised** Parquet (`data/clean.parquet`).
-   Records are pseudonymised; this repo operates only on pseudonymous IDs and does not use direct identifiers.
-   CI (GitHub Actions) does **not** access Google Drive or patient data.


## Provenance & timelines
- Report submitted: 18 Aug 2025. Slides submitted: 25 Aug 2025. Viva: 29 Aug 2025.
- This public repository was created after report submission for reproducibility.
- Analysis code matches the Appendix version submitted on 18 Aug (non-analytic files like README/CI may differ).
- Integrity: see [`SHA256SUMS.txt`](./SHA256SUMS.txt) at the repo root for file checksums.



------------------------------------------------------------------------

## Quick start

### Local R with renv

Point `RAW_FILE` at the **Google Drive–mounted path** and run the pipeline.

``` r
install.packages("renv")
renv::restore()  # recreate the package library from renv.lock

# Windows (G: as an example)
Sys.setenv(RAW_FILE = "G:/Shared drives/<TEAM_OR_PROJECT>/<FOLDER>/FINAL BAYESIAN_07042025 (1).xlsx")

# macOS (/Volumes/GoogleDrive)
# Sys.setenv(RAW_FILE = "/Volumes/GoogleDrive/Shared drives/<TEAM_OR_PROJECT>/<FOLDER>/FINAL BAYESIAN_07042025 (1).xlsx")

# Build the full pipeline and render the thesis
targets::tar_make(thesis_html)
# …or just build model artefacts without rendering
# targets::tar_make()
```

**Windows note**: Accents in file paths can break I/O in some setups. An early `file.exists()` guard is included; if issues arise, move or symlink the raw file to an ASCII path and update `RAW_FILE`.

------------------------------------------------------------------------

## Configuration

Set the raw workbook path once per session **to the Google Drive–mounted path**:

``` r
# Windows (G: as an example)
Sys.setenv(RAW_FILE = "G:/Shared drives/<TEAM_OR_PROJECT>/<FOLDER>/FINAL BAYESIAN_07042025 (1).xlsx")

# macOS (/Volumes/GoogleDrive)
# Sys.setenv(RAW_FILE = "/Volumes/GoogleDrive/Shared drives/<TEAM_OR_PROJECT>/<FOLDER>/FINAL BAYESIAN_07042025 (1).xlsx")
```

Defaults and thresholds live in `R/00_config.R`:

-   Study collapsing thresholds: `collapse_min_alive = 30`, `collapse_min_dead = 5` (small groups pooled into “Merged-small-sites”).
-   Validation horizons in `_targets.R`:

``` r
TIME_HORIZONS <- c(0.5, 1, 3, 7, 15, 30, 60, 90)
```

------------------------------------------------------------------------

## Running the pipeline

Common targets:

``` r
targets::tar_make()               # end-to-end model artefacts (no rendering)
targets::tar_make(thesis_html)    # full run + Quarto report (HTML)
targets::tar_visnetwork()         # view DAG
```

Recover space if figures / objects are large:

``` r
# Drop superseded cache metadata (safe)
targets::tar_prune()

# Remove stored objects to reclaim disk (forces recompute on next run)
targets::tar_destroy(destroy = "objects")
```

## Modeling overview

### Primary analysis (Bayesian, INLA)

-   Cause-specific hazards for **ICU death** vs **discharge alive**.
-   Transforms: `log(eGFR)` and `log(BMI)`.
-   **Age**: natural cubic splines (**3 df**), centered at the dataset mean (`age_mean` stored in the fit).
-   **Study-group random effects** (IID intercepts) to capture residual heterogeneity across sub-studies.
-   **Piecewise-exponential baseline** on a shared time grid with up to **5 segments** and **≥10 failures per cause per band**, smoothed with **RW2**.
-   **Priors**: weakly-informative Gaussians for fixed effects (sd ≈ 5); PC/RW2 priors for baseline segments and random-effect precisions.

### Fine–Gray engine selection

-   Primary inference relies on **cause-specific** models.
-   Fine–Gray is used **for descriptive comparisons** of CIFs and ranking performance.
-   Implementation choices:
    -   `riskRegression::FGR` when flexible interfaces/robust SEs are needed.
    -   `cmprsk::crr` when a simple interface suffices (no clustering term). In the **index-stay dataset** there is effectively **no administrative right-censoring**, so FG results are ancillary and used descriptively.

### Multicollinearity policy (eGFR & creatinine)

-   **Primary models** use **log-eGFR** only.
-   **Sensitivity**: add the studentised residual of `log(creatinine)` after regressing on age, sex, and `log(eGFR)` (`add_creatinine_residual()`), ensuring orthogonality.

### Validation & performance metrics

-   **Time-dependent AUC** and **Brier score** at specified horizons, averaged across imputations (bootstrap CIs optional).
-   **Internal validation** via apparent performance with full calibration plots; no cluster bootstrap is reported.
-   **Decision-curve analysis** evaluated over the empirical risk range (thresholds confined to the range of predicted risks).

------------------------------------------------------------------------

## INLA ⇄ MCMC agreement check

-   A JAGS mirror of the INLA fixed-effect structure is provided (log transforms, shared `age_mean`, same knots/cuts).
-   When `config=TRUE`, **INLA posterior draws** (`inla.posterior.sample`) enable contrast uncertainty checks.
-   Summaries include alignment of linear predictors and convergence diagnostics (max **R-hat**, **ESS**).

Run (wired in `_targets.R`):

``` r
targets::tar_make(mcmc_check)   # classic mapping/means
targets::tar_make(mcmc_agree) # richer summary + convergence (see functions_mcmc_check.R)
```

## Tests & CI

Run locally:

``` r
library(testthat)
testthat::test_dir("tests/testthat", reporter = "summary")
```

## Reproducibility {#reproducibility}

-   The renv lock ensures deterministic R package versions; the CI workflow pins a Linux image for reproducible runs.
-   Global seeds are set where appropriate (`_targets.R` and function scope).
-   `targets` caches immutable artefacts by hash; rebuilds occur only when inputs change.

------------------------------------------------------------------------

## Troubleshooting

-   **Raw file not found**: Check `RAW_FILE`; avoid non-ASCII characters where possible.
-   **FGR does not fit**: No administrative censoring in the index-stay dataset. Use the `crr` path for descriptive FG; primary inference remains cause-specific.
-   **Large figures (\>1–2 GB)**: Export TIFF/PNG via `export_figure()` and prune cache with `targets::tar_prune()` or remove stored objects via `targets::tar_destroy(destroy = "objects")`.
-   **Prediction NA**: Guard code returns placeholders when a model is not fitted or inputs are missing; logs indicate the reason.

------------------------------------------------------------------------

## Citations

-   Rue, H., Martino, S., & Chopin, N. (2009). **Approximate Bayesian inference for latent Gaussian models by using INLA**. *JRSS-B*.

-   Fine, J. P., & Gray, R. J. (1999). **A proportional hazards model for the subdistribution of a competing risk**. *JASA*.

-   Andersen, P. K., & Keiding, N. (2002). **Multi-state models for event history analysis**. *Statistical Methods in Medical Research*.

-   van Buuren, S., & Groothuis-Oudshoorn, K. (2011). **mice**: Multivariate Imputation by Chained Equations in R.

-   Aalen, O. O., & Johansen, S. (1978). **An empirical transition matrix for non-homogeneous Markov chains based on censored observations**. *Scand J Stat*.

-   Simpson, D., Rue, H., Riebler, A., Martins, T. G., & Sørbye, S. H. (2017). **Penalising model component complexity**. *Statistical Science*.

## License

-   Code: MIT
-   No real patient data is included in this repository.

------------------------------------------------------------------------

### Maintainers

-   PI / ICU: Dr Björn Stessel
-   Clinical supervisor: Dr Ina Callebaut
-   Academic supervisor: Dr Oswaldo Gressani
-   Author: Henry Hockliffe
