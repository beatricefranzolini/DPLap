# DPLap

R code implementing Laplace and skew-Laplace approximations for posterior density estimation in Dirichlet process mixtures of univariate Gaussian kernels. 
The repository also includes the variational Bayes and slice-sampling methods used for comparison, together with the scripts that reproduce all numerical results in the paper.

**Authors**: [Beatrice Franzolini](https://beatricefranzolini.github.io/) and [Francesco Pozza](https://francesco16p.github.io/)

## Citation

If you use this repository in your research, please cite:

> Franzolini, B. and Pozza, F. (2026). Laplace and skew-Laplace approximations for Dirichlet process mixture posterior density. *Measurement: Interdisciplinary Research and Perspectives*, in press.

A preprint version is available at [arXiv:2604.25410](https://beatricefranzolini.github.io/FranzoliniPozza2026.pdf).

```bibtex
@article{franzolini2026laplace,
  author  = {Franzolini, Beatrice and Pozza, Francesco},
  title   = {Laplace and skew-{L}aplace approximations for {D}irichlet process mixture posterior density},
  journal = {Measurement: Interdisciplinary Research and Perspectives},
  year    = {2026},
  note    = {In press}
}
```

## Repository structure

```
DPLap/
├── dplap.R                       # Laplace and skew-Laplace approximations
├── utils_dplap.R                 # utilities for the Laplace-based methods
├── dpvb.R                        # variational Bayes (Blei and Jordan, 2006)
├── dpslice.R                     # slice sampler (Franzolini and Gaffi, 2026)
├── utils_slice.R                 # utilities for the slice sampler (Franzolini and Gaffi, 2026)
├── plot_top2_atoms.R             # plots of the atoms with the two largest weights
├── sim1/, ..., sim4/             # simulation scenarios 1–4, α ~ Gamma(3, 6)
├── sim1_logn/, ..., sim4_logn/   # simulation scenarios 1–4, α ~ Gamma(3, 3 log n)
└── real_data/                    # real-data examples
```

Each experiment folder contains its main script together with the outputs reported in the paper.

### Main scripts

The following scripts reproduce the numerical studies in the paper.

| Script | Setting |
| --- | --- |
| `sim1/_main_sim1.R` | Scenario 1: four-component Gaussian mixture |
| `sim2/_main_sim2.R` | Scenario 2: 100-component Gaussian mixture with Zipf-type weights |
| `sim3/_main_sim3.R` | Scenario 3: as Scenario 1, with Student-t kernels (5 degrees of freedom) |
| `sim4/_main_sim4.R` | Scenario 4: as Scenario 2, with Student-t kernels (5 degrees of freedom) |
| `sim1_logn/_main_sim1_logn.R`, ..., `sim4_logn/_main_sim4_logn.R` | Scenarios 1–4 under the alternative prior α ~ Gamma(3, 3 log n) |
| `real_data/_main_real_data.R` | `faithful`, `galaxies`, `iris`, `rock` and `quakes` datasets |

The `simX/` scripts use α ~ Gamma(3, 6), the setting reported in the main text, while the `simX_logn/` scripts repeat the same scenarios with α ~ Gamma(3, 3 log n), reported in the Supplementary Material.

Each simulation script generates data with n ∈ {20, 50, 100, 200, 500, 1000, 1500, 2000}, fits the variational Bayes, Laplace and skew-Laplace approximations, runs the slice sampler, and writes the following files to its folder.

| File | Content | In the paper (`simX/`) | In the paper (`simX_logn/`) |
| --- | --- | --- | --- |
| `time.csv` | Run times, in seconds, of the four methods | Table 1 | Supplement S.2 |
| `tv_point_estimate.csv` | Total variation distances of the posterior mean densities from the true density and from the slice-sampler estimate | Table 2 | — |
| `Boxplots_TVerror_densityordinatesposteriors.pdf` | Pointwise total variation distances between the posterior distributions of the density ordinates under each approximation and under the slice sampler | Figure 2 | Supplement S.2 |
| `bivariate_atom_n100.pdf`, `bivariate_atom_n1500.pdf` | Joint posterior of the atoms with the two largest weights, for n = 100 and n = 1500 | Figures 3–6 | Supplement S.2 |
| `ESS.csv` | Effective sample sizes of the slice sampler, based on the log-likelihood | Supplement S.4 | Supplement S.4 |
| `densities_point_estimate.pdf` | Posterior mean densities over the data histograms | Figure 1, Supplement S.1 | Supplement S.3 |

The density plots are drawn in the active graphics device at the end of each simulation script, and the full workspace is saved as an `.RData` file in the same folder.

The real-data script analyses the standardised `faithful`, `galaxies`, `iris`, `rock` and `quakes` data, prints the run times and total variation distances of Table 3 to the console, and draws Figure 7 (`real_data/density_real.pdf`) and Figure 8 (`real_data/Boxplot_real.pdf`) in the active graphics device.

### Core functions

- `dplap.R` — `dp_lap_mixmodel_normal_normal()` locates the posterior mode of the truncated DPM by multistart BFGS on the unconstrained scale (logit stick-breaking variables, atom locations, log α), merges atoms closer than `tol_merge` and refits, then samples from the Laplace (`skew = FALSE`) or skew-Laplace (`skew = TRUE`) approximation.
- `utils_dplap.R` — log-posterior and gradient of the truncated DPM with random α, the skew-symmetric sampler (`sim_eqc_skew_sym()`, Algorithm 1 in the paper), posterior density summaries on a grid (`posterior_density_dpm()`) and total variation utilities.
- `dpvb.R` — `dp_vb_mixmodel_normal_normal()` implements truncated mean-field variational Bayes (Blei and Jordan, 2006) with a Gamma variational factor for α; `posterior_density_dpm_vb()` returns the corresponding density summaries.
- `dpslice.R` — `dp_slice_mixmodel_normal_normal()` implements the slice sampler used as benchmark, which targets the untruncated DPM posterior with random α.
- `utils_slice.R` — utilities for the slice sampler, including posterior density estimation with completion of the unoccupied stick-breaking tail (`estimate_dpmm_density_ss()`) and the log-likelihood trace used for effective sample sizes.
- `plot_top2_atoms.R` — `plot_top2_atoms()` plots the joint posterior of the atoms with the two largest weights (ranked within each draw) for any combination of variational Bayes, Laplace, skew-Laplace and slice-sampler fits.

### External code

- `dpslice.R` and `utils_slice.R` — adapted from the code of Franzolini, B. and Gaffi, F. (2026), *Complexity bounds for Dirichlet process slice samplers*. 43rd International Conference on Machine Learning (ICML 2026).

## Required R packages

The experiments in the paper were run with R 4.6.1 and the package versions below.

### Packages used by the code

- `optimx` (version 2025-4.9) — multistart BFGS search for the posterior mode
- `mvtnorm` — sampling from the Gaussian Laplace approximation
- `matrixStats` (version 1.5.0) — loaded by `utils_dplap.R`
- `LaplacesDemon` (version 16.1.6) — Dirichlet draws in the slice sampler and effective sample sizes
- `progress` (version 1.2.3) — progress bar of the slice sampler
- `ggplot2` (version 4.0.1) — boxplots and bivariate atom plots
- `MASS` — `galaxies` dataset and two-dimensional density contours in `plot_top2_atoms()` (included in standard R installations)

### Additional packages loaded by the main scripts

- `rstudioapi` (version 0.17.1) — sets the working directory to the folder of the script
- `patchwork` (version 1.3.2), `fossil` (version 0.4.0) and `salso` (version 0.3.53)

## Installation

Clone the repository

```bash
git clone https://github.com/beatricefranzolini/DPLap.git
```

and install the required packages, for example with

```r
install.packages(c("optimx", "mvtnorm", "matrixStats", "LaplacesDemon", "progress",
                   "ggplot2", "rstudioapi", "patchwork", "fossil", "salso"))
```

## Getting started

To reproduce an experiment, open its main script (for example `sim1/_main_sim1.R`) in RStudio and run it directly or source it. 
To run a script outside RStudio, replace the line `setwd(dirname(rstudioapi::getActiveDocumentContext()$path))` with a call to `setwd()` pointing to the folder of the script.

Simulated data are generated within each script with fixed seeds, and the real datasets are loaded from the `datasets` and `MASS` packages, so no data download is needed.

Most of the run time is spent in the slice sampler (100,000 iterations for each dataset). 
On the laptop used for the paper (Intel Core Ultra 7 268V, 32 GB of RAM), the slice sampler alone took between about 35 and 60 minutes per simulation script and about 20 minutes for the real-data script. 
For a quick check, reduce `MCMC_draws` and `Burn_in` at the top of a script, keeping `MCMC_draws - Burn_in` at least equal to `MC_samples`.

## Using the approximations on your own data

The Laplace and skew-Laplace approximations only require `utils_dplap.R` and `dplap.R`. From the root of the repository:

```r
library(optimx)  # provides multistart()
source("utils_dplap.R")
source("dplap.R")

y <- scale(faithful$eruptions)[, 1]  # standardised data

fit <- dp_lap_mixmodel_normal_normal(
  y,
  K      = floor(5 * log(length(y))),  # initial truncation level
  sigma  = 0.5,                        # kernel standard deviation
  m0     = 0, s0 = 0.25,               # base measure N(m0, s0^2)
  a      = 3, b = 3,                   # alpha ~ Gamma(a, b), with b a rate
  skew   = TRUE,                       # FALSE gives the Gaussian Laplace approximation
  n_samp = 2000,
  seed   = 0
)

xgrid <- seq(min(y) - 0.5, max(y) + 0.5, length.out = 400)
dens  <- posterior_density_dpm(xgrid, fit$samples, sigma = 0.5)

plot(xgrid, dens$mean_density, type = "l", ylim = c(0, max(dens$quantiles)),
     xlab = "y", ylab = "density")
lines(xgrid, dens$quantiles[1, ], lty = 2)  # pointwise 2.5% quantile
lines(xgrid, dens$quantiles[3, ], lty = 2)  # pointwise 97.5% quantile
```

The draws in `fit$samples` are on the unconstrained scale (logit stick-breaking variables, atom locations and log α), and `fit$K` is the truncation level after merging.

## Questions or bug reports

For questions or bug reports, please contact:

- [Beatrice Franzolini](https://beatricefranzolini.github.io/) — franzolini@pm.me
- [Francesco Pozza](https://francesco16p.github.io/) — francesco.pozza2@unibocconi.it
