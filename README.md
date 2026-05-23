# Reproducible simulation code

This repository contains reproducible R code for a simulation study involving synthetic family-structured survival data.

The objective of the simulation is to illustrate the performance of the proposed statistical model and to show how the estimation procedure can be implemented in practice. The proposed model accounts for intra-family dependence using copulas, incorporates time-dependent covariate effects through B-splines, and uses conditional likelihood contributions to correct for selection bias induced by family ascertainment.

The simulated data are synthetic and do not correspond to real clinical or genetic data. They are designed to reproduce a family-structured survival setting in which individuals within the same family may exhibit dependent time-to-event outcomes.

## Objectives of the simulation

The simulation study aims to demonstrate the following components of the proposed methodology:

1. The ability of the proposed model to estimate survival-related parameters from family-structured data.

2. The use of a copula-based dependence structure to model intra-family association between event times.

3. The integration of B-spline functions to model a time-dependent covariate.

4. The implementation of an iterative optimization procedure for estimating marginal and dependence parameters.

5. The use of conditional likelihood contributions to account for selection bias due to ascertainment through affected probands.

6. The calculation of asymptotic variance estimates for statistical inference.

## Repository structure

- `R/`: R functions used for data generation, likelihood evaluation, copula calculations, optimization, and variance estimation.
- `reproduce_one_simulation.R`: main script used to generate one reproducible simulated data set and estimate the model.
- `output/`: folder where generated results are saved locally.
- `sessionInfo.txt`: information about the R session used to run the code.

## Computational time

For each value of the penalty parameter gamma, parameter estimation and robust asymptotic variance calculation take approximately 20 minutes on average for one simulated data set. The total running time depends on the computer and R environment.

## How to run the code

To reproduce the simulation example, run:

```r
source("reproduce_one_simulation.R")


