# RBORM: Robust Bayesian Ordinal Regression Model via Divergence Approach

This repository provides the R code implementing the robust bayesian inference for the ordinal response model via the divergence approach, as proposed by the following paper.

Momozaki, T. and Nakagawa, T. (2023). Robustness of Bayesian ordinal response model against outliers via divergence approach. https://arxiv.org/abs/2305.07553

The repository includes the following files.
- `robustORM_WLB.R`: The script implementing the proposed bayesian inference methods for the ordinal response model with the density-power and $\gamma$-divergences using the weighted likelihood bootstrap (Newton and Raftery, 1994, JRSSB)
- `robustORM_stan.R`: The script implementing the proposed bayesian inference methods for the ordinal response model with the density-power and $\gamma$-divergences using rstan
- `links.R`: The script containing functions that compose several density and distribution functions
- robustORM_stan folder: The folder containing the Stan files needed to run the robustORM_stan.R file
