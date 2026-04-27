# DPLap

R code implementing Laplace and Skew-Laplace approximations for posterior density estimation in Dirichlet process mixtures of univariate Gaussian kernels.

**Authors**: [Beatrice Franzolini](https://beatricefranzolini.github.io/) and [Francesco Pozza](https://francesco16p.github.io/)

## Citation

If you use this repository in your research, please cite:

> Franzolini, B. and Pozza, F. (2026). *Laplace and skew-Laplace sampling for posterior approximations in Dirichlet process mixtures*.

## Repository structure

### Main scripts

The following scripts reproduce the numerical studies in the paper.

> **Important**  
> These scripts should be run directly and **not sourced**.

- `_main_real_data.R` — runs the numerical experiments on the real datasets.
- `_main_sim1.R` — runs the numerical experiment for Simulation Scenario 1.
- `_main_sim2.R` — runs the numerical experiment for Simulation Scenario 2.
- `_main_sim3.R` — runs the numerical experiment for Simulation Scenario 3.
- `_main_sim4.R` — runs the numerical experiment for Simulation Scenario 4.

### Core functions

- `dplap.R` — main functions implementing the Laplace and Skew-Laplace approximations.
- `utils_dplap.R` — utility functions used throughout the repository.

### External code

- `DPalg/` — code from Franzolini, B. and Gaffi, F. (2026), *Complexity bounds for Dirichlet process slice samplers*.  
  [arXiv:2602.00878](https://arxiv.org/abs/2602.00878)

## Required R packages

### Packages used by this repository

- `rstudioapi` (version 0.17.1)
- `ggplot2` (version 4.0.1)
- `matrixStats` (version 1.5.0)

### Additional packages required by `DPalg/`

- `LaplacesDemon` (version 16.1.6) — sampling from finite Dirichlet distributions
- `progress` (version 1.2.3) — progress bars
- `rstudioapi` (version 0.17.1) — working directory handling
- `ggplot2` (version 4.0.1) — plotting
- `patchwork` (version 1.3.2) — combining plots
- `fossil` (version 0.4.0) — Rand index computation
- `salso` (version 0.3.53) — posterior similarity matrices and clustering point estimates

## Installation

Clone the repository and make sure all required packages are installed before running the scripts.

## Getting started

To reproduce the experiments, open one of the main scripts (for example `_main_sim1.R`) in RStudio and run it directly.

## Questions or bug reports

For questions or bug reports, please contact:

- [Beatrice Franzolini](https://beatricefranzolini.github.io/) — franzolini@pm.me
- [Francesco Pozza](https://francesco16p.github.io/) — francesco.pozza2@unibocconi.it
