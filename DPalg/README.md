## Full repository available at [https://github.com/beatricefranzolini/DPalg](https://github.com/beatricefranzolini/DPalg)

# DPalg
R codes implementing MCMC algorithms for DP mixture of univariate Normals

| Level | Distribution |
|---|---|
| Kernel | `y_i \| μ_i, σ² ~ Normal(μ_i, σ²)` |
| Mixing distribution | `μ_i \| G ~ G` |
| Dirichlet process prior | `G \| α ~ DP(α, G₀)` |
| Base measure | `G₀ = Normal(μ₀, τ₀²)` |
| Concentration parameter | `α ~ p(α)` |

## Citation 
Please cite the following publication if you use this repository in your research:

> [Franzolini, B.](https://beatricefranzolini.github.io/) & [Gaffi, F.](https://francescogaffi.github.io/), (2026). *Complexity bounds for Dirichlet process slice samplers*. 43rd International Conference on Machine Learning (ICML 2026)

## Contents 

`_main_to_run1.R` runs the numerical experiment with 3 equally-sized clusters

`_main_to_run2.R` runs the numerical experiment with perturbed zipf

`BGS.R`           MCMC function to run the block Gibbs sampler

`CRPnoAtoms.R`    MCMC function to run the CRP with atoms

`CRPwithAtoms.R`  MCMC function to run the CRP with atoms

`DPalg.R`         this code sources the required functions and packages it is called by _main_to_run1.R and _main_to_run2.R

`Slice.R`         MCMC function to run the slice sampler

`SlicenoAtoms.R`  MCMC function to run the slice sampler with no atoms

`utils`           contains helper functions for the MCMCs

## R libraries needed to run the scripts

- `LaplacesDemon` (version 16.1.6) — sampling from finite Dirichlet distributions
- `progress` (version 1.2.3) — progress bars
- `rstudioapi` (version 0.17.1) — working directory handling
- `ggplot2` (version 4.0.1) — plotting
- `patchwork` (version 1.3.2) — combining plots
- `fossil` (version 0.4.0) — Rand index computation
- `salso` (version 0.3.53) — posterior similarity matrices and clustering point estimates

## Questions or bugs

For bug reporting purposes, e-mail [Beatrice Franzolini](https://beatricefranzolini.github.io/) (franzolini@pm.me)

