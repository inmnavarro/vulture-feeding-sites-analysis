# Supplementary feeding sites and facultative scavengers — R analysis code

**Code repository for a study of carcass deposition points and conventional feeding stations in central Spain.**

This is a **code-only** repository prepared to document an observational camera-trap analysis for peer review. It contains the R analysis, posterior prediction, figure-generation, and diagnostic scripts. **It does not include the camera-trap records, hunting-estate identifiers, estate coordinates, photographs, or fitted `.rds` objects.** Fitted `brms` objects can embed underlying data and are intentionally excluded.

## Study overview

- 20 hunting estates, 64 identifiable carcass-provisioning events.
- Two facility categories: big feeding stations and carcass deposition points.
- Five Bayesian models: facultative-scavenger detection (Bernoulli), operational site safety (Bernoulli), first arrival (Bernoulli), time to first vulture arrival (Gamma), and vulture attendance per observation event (negative binomial).
- The *operational site-safety model* in the saved analysis is based on **red-fox-versus-vulture timing**; the other facultative-scavenger responses include foxes, dogs, and mustelids as coded in the source data. This distinction should not be lost when reporting the models.
- Field records are observations, not individually independent photos. Multiple monitoring periods may be associated with a provisioning.

## Files

| Script | Purpose |
| --- | --- |
| `scripts/01_fit_bayesian_models.R` | Rebuild five Bayesian fits from the private dataset |
| `scripts/02_manuscript_figures.R` | Generate main/supplementary figure suite and posterior summaries from saved final fits |
| `scripts/03_temporal_vulture_figure.R` | Produce two-panel within-season figure from saved predictions and fitted model data |
| `scripts/04_extract_hierarchical_effects.R` | Extract posterior intervals for hunting-estate random intercept and smooth-term SD |
| `scripts/05_observed_consumption.R` | Plot descriptive presence and consumption proportions |

The plot scripts include the descriptive correction that **observations explicitly noted as eating vegetation are not counted as consumption of carcass remains**, even if the raw behaviour category is `Consumo`.

## How to run

1. Install a recent version of R and the packages used in the scripts (`brms`, `readr`, `dplyr`, `lubridate`, `ggplot2`, `bayesplot`, `patchwork`, `tidyr`, `scales`). Bayesian model fitting also requires a functional Stan/C++ toolchain.
2. Place the privately supplied `BASE_DATOS_BUITRES.csv` in a local folder. The source is a semicolon-separated file with Latin-1 encoding and contains sensitive hunting-estate information; it is **not published here**.
3. Set the environment variable `SFS_DATA_DIR` to that folder. Optionally set `SFS_MODEL_DIR` to the folder containing the saved models (default: `<SFS_DATA_DIR>/models`). On Windows/RStudio, use `Sys.setenv(SFS_DATA_DIR="C:/path/to/private/data")`.
4. Run `scripts/01_fit_bayesian_models.R` if re-fitting is intended; outputs are saved locally as `.rds` objects. The three previously divergent models use `adapt_delta = 0.99` and `max_treedepth = 15`, as used in the revised runs.
5. Run `scripts/02_manuscript_figures.R` first. It generates the posterior predictions consumed by `scripts/03_temporal_vulture_figure.R`. The two supplementary scripts may be run independently after fitting.

**Reproducibility limitation:** Without access to the curated source dataset and exact R/Stan package environment, this public code can be inspected but cannot independently reproduce the numerical analyses in the manuscript. Code sharing must not be mistaken for deposition of original research data.

## Ethical and privacy precautions

- Never upload raw data, images, location coordinates, field notes, hunting-estate names, owner names, or `.rds` fits to this repository.
- Generated tables may contain row-level hunting-estate or provisioning identifiers. The `.gitignore` intentionally excludes all CSVs and model objects.
- Confirm coauthor agreement and the journal's anonymous-review requirements before making the GitHub repository or its link public.
- To share software during double-anonymized peer review, consider an anonymous repository mirror such as [Anonymous GitHub](https://anonymous.4open.science/) and inspect the mirror for identifying material before providing the link.
- The publisher's **original-data repository** field should not be completed with this software repository URL unless the submission system explicitly accepts *code* as the deposited resource. Data availability and code availability are separate statements.

## Licensing and attribution

No open-source license is granted in this initial review package. The authors should select a license and attribution terms together before public release. Authorship details are deliberately not included in this version to facilitate anonymous peer review.
