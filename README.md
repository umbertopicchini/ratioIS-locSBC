# Ratio-corrected importance sampling and localized SBC

This repository contains the R code accompanying

> **U. Picchini** (2026). [*Amortized ratio-estimation importance sampling and localized simulation-based calibration for intractable likelihoods*](https://arxiv.org/abs/2609.39712), arXiv:2609.39712.

The code implements the methodology studied in the paper: surrogate-based importance sampling corrected by amortized likelihood-ratio estimation (ratioIS), together with localized simulation-based calibration (SBC). It also contains the simulation models and scripts used for the numerical examples and applications.

## Repository contents

The three main run scripts correspond to:

- a stochastic SIR example;
- a stochastic Lotka--Volterra example;
- the ST131 epidemiological application.

The `results/` directory is created by the run scripts and stores fitted surrogate objects, ratioIS output, posterior draws, diagnostics, localized-SBC output, and some figures. All figures can be obtained by running parts of the scripts.

## Main workflow

Each example follows the same broad workflow:

1. simulate a prior-predictive training sample and fit the initial surrogate;
2. construct focused simulation rounds and fit the final proposal surrogate;
3. draw from the final surrogate proposal;
4. estimate the required ratio correction using the classifier-based ratio estimator;
5. form the corrected importance-sampling posterior, including the clipped-weight version used for the main reported results;
6. optionally run localized SBC around the fitted posterior;
7. produce the example-specific posterior and posterior-predictive plots.

The run scripts keep the settings that are useful to vary between experiments explicit, while stable implementation choices are handled inside the inference routines.

## Software requirements

The code is written in R.

The user should install the `xLLiM` package via devtools::install\_github("epertham/xLLiM", ref = "master")
The paper used version 2.3.3 of `xLLiM`.

The examples use the following packages:

```r
install.packages(c(
  "deSolve",
  "glmnet",
  "loo",
  "MASS",
  "mixtools",
  "mvtnorm",
  "readxl",
  "rjson"
))
```

The Lotka--Volterra example additionally requires `smfsb`, and the surrogate fitting code requires `xLLiM`.

Please install these packages separately if they are not already available in your R installation.

## Running the examples

Run scripts from the repository root so that their relative `src/`, `models/`, and `results/` paths resolve correctly.

### SIR example

```r
source("sir_run.R")
```

### Stochastic Lotka--Volterra example

```r
source("lotka_volterra_run.R")
```

### ST131 epidemiological application

```r
source("ecoli_st131_scaledDt.R")
```

For the ST131 application, set the clade near the top of the script:

```r
clade_to_fit <- "A"
```

with `"C1"` and `"C2"` available for the corresponding analyses.

The ST131 script performs surrogate fitting in the numerically scaled parameterization

```text
(tau, R0, Dt_scaled),   Dt_scaled = 0.001 * Dt,
```

while the epidemiological simulator receives `Dt` on its original scale in weeks.

## Important run settings

The main scripts expose the experiment-level choices directly, including:

```r
K_start
cov_structure
eps_defensive_stage1
eps_defensive_stage2
local_sbc_N

feature_args <- list(
  include_linear = TRUE,
  include_quadratic = FALSE,
  include_within_block_interactions = FALSE,
  include_cross_interactions = TRUE
)
```

The reported classifier fits use 10-fold cross-validation and the `lambda.1se` rule. Main RatioIS inference does not perform an additional observation-specific component-pruning step. Localized SBC separately uses an active-component threshold when constructing the local surrogate.

The primary corrected posterior reported in the paper uses deterministic clipped importance weights with clipping level

```text
M_c = floor(sqrt(M)),
```

where `M` is the number of final importance samples.

## Simulation budgets used in the paper

The cleaned scripts retain the simulation budgets used for the reported experiments:

| Example | Initial prior-predictive | Focused round | Refreshed round | Final IS draws |
|---|---:|---:|---:|---:|
| SIR | 20,000 | 20,000 | 10,000 | 10,000 |
| Lotka--Volterra | 20,000 | 20,000 | 10,000 | 10,000 |
| ST131 | 50,000 | 30,000 | 20,000 | 10,000 |

Localized SBC uses the inflation factors

```r
c(1, 1.05, 1.1, 1.25, 1.5, 2)
```

with the example-specific settings given in the corresponding run scripts.

## Output

Depending on the example, the scripts save objects such as

```text
surrogate_out.rds
ratio_out.rds
proposal_draws_q2.csv
posterior_draws_corrected.csv
posterior_draws_corrected_clipped.csv
diagnostics.csv
settings.json
```

as well as localized-SBC results and posterior/posterior-predictive figures.

The clipped posterior draws are the primary ratioIS posterior sample used in the reported comparisons and posterior predictive checks.

## Reproducibility

Random-number seeds are set explicitly in the run scripts. The localized-SBC code also uses deterministic seed conventions for the separate inflation factors (and, where applicable, observations or ST131 clades).

Because the examples contain stochastic simulation, exact reproduction requires using the same code version, data files, package versions, seeds, and simulation settings.

## Data for the ST131 application

The ST131 script uses data from [Ojala et al 2025](https://doi.org/10.1038/s41467-025-65301-1) available at https://github.com/fanoja/Ecoli-R0

```text
models/ecoli_st131/data/NORM_incidence_updated.csv
models/ecoli_st131/data/ST131_clades_OR_E_coli_carriage_disease_collapsed_wBSAC2.xlsx
```

The run script reads these files and constructs the clade-specific model inputs before calling the simulator.

## Citation

If you use this code, please cite:

```bibtex
@article{picchini2026ratiois,
  author  = {Picchini, Umberto},
  title   = {Amortized ratio-estimation importance sampling and localized simulation-based calibration for intractable likelihoods},
  year    = {2026},
  journal = {arXiv preprint arXiv:2609.39712},
  url     = {https://arxiv.org/abs/2609.39712}
}
```

## Paper

The preprint is available at:

https://arxiv.org/abs/2609.39712

