
############################################################
## Reproducible simulation and estimation script
## 
## This script generates one simulated data set and estimates
## the proposed model for several penalty parameters:
## gamma = 0, 0.05, 0.2, 0.8.
##
## The objective is to provide reviewers with a complete
## reproducible example of the simulation, estimation and
## variance calculation procedures.
############################################################

rm(list = ls())

############################################################
## 0. Required packages
############################################################

library(survival)
library(numDeriv)
library(Matrix)
library(fda)

############################################################
## 1. Source functions
############################################################

source("R/Distribution_Marginales_Simulation_Data.R")
source("R/Distribution_marginales_10.R")
source("R/Generate.Survival_g.R")
source("R/kinship.R")
source("R/indices.R")
source("R/DerivCop.R")
source("R/LogLik.Uni_2026.R")
source("R/LogLik.Bivariete.R")
source("R/Optimisation_itt_2026.R")
source("R/LogLik.Uni_vec.R")
source("R/LogLik.Bivariete_vec.R")
source("R/D2.R")

############################################################
## 2. Create output folders
############################################################

if (!dir.exists("output")) {
  dir.create("output")
}

if (!dir.exists("output/data")) {
  dir.create("output/data")
}

if (!dir.exists("output/estimates")) {
  dir.create("output/estimates")
}

if (!dir.exists("output/variances")) {
  dir.create("output/variances")
}

############################################################
## 3. Reproducibility seed
############################################################

set.seed(12345)

############################################################
## 4. Model parameters and simulation settings
############################################################

lambda <- 1 / 1000
l.lambda <- log(lambda)

alpha <- 4
l.alpha <- log(alpha)

beta1 <- 1
beta2 <- -1

h <- 0.5

m.a <- 4
v.a <- 2

p <- 0.5

I <- 200

min.c <- 0
max.c <- 10

beta_true <- c(
  0.5443045, 0.8185666, 1.0198131,
  0.2207210, 0.5481451, 0.9364897,
  0.4976894, 0.5799883, 0.7871138,
  0.7813054
)

params <- c(
  l.alpha,
  l.lambda,
  beta1
)

params0 <- c(
  l.alpha,
  l.lambda,
  beta1,
  beta_true
)

############################################################
## 5. Create indices
############################################################

indices.proband <- creer.indices.proband(I)

indices.non.proband <- creer.indices.non.proband(I)

indices.bivariee <- creer.indices.bivariee(I)

############################################################
## 6. Penalty matrix for P-splines
############################################################

pspline_penalty_core <- function(K, d = 1) {
  
  if (K <= d) {
    stop("K must be greater than d.")
  }
  
  D <- Diagonal(K)
  
  for (i in seq_len(d)) {
    D <- diff(D)
  }
  
  P <- crossprod(D)
  
  return(P)
}

P1 <- pspline_penalty_core(K = 10, d = 2)

############################################################
## 7. Utility functions
############################################################

drop_col_by_name <- function(df, colname) {
  
  if (!is.null(df) && colname %in% names(df)) {
    df[, setdiff(names(df), colname), drop = FALSE]
  } else {
    df
  }
}

extract_par <- function(obj) {
  
  if (is.null(obj)) {
    stop("compute.mle returned NULL.")
  }
  
  if (is.numeric(obj)) {
    return(obj)
  }
  
  if (is.list(obj)) {
    
    if (!is.null(obj$par) && is.numeric(obj$par)) {
      return(obj$par)
    }
    
    if (!is.null(obj$theta) && is.numeric(obj$theta)) {
      return(obj$theta)
    }
    
    if (!is.null(obj$params) && is.numeric(obj$params)) {
      return(obj$params)
    }
  }
  
  stop("The output of compute.mle is not recognized.")
}

safe_gamma_label <- function(gamma_value) {
  
  label <- switch(
    as.character(gamma_value),
    "0" = "000",
    "0.05" = "005",
    "0.2" = "020",
    "0.8" = "080",
    gsub("\\.", "", as.character(gamma_value))
  )
  
  return(label)
}

############################################################
## 8. Generate one simulated data set
############################################################

cat("Generating one simulated data set...\n")

data <- Generate.data(
  I = I,
  params = params,
  m.a = m.a,
  v.a = v.a,
  min.c = min.c,
  max.c = max.c,
  p = p,
  h = h,
  kinship = kinship
)

saveRDS(
  data,
  file = "output/data/simulated_data.rds"
)

############################################################
## 9. Construct analysis data sets
############################################################

cat("Constructing analysis data sets...\n")

data.proband <- create.data.proband(
  data = data,
  indices.proband = indices.proband
)

data.non.proband <- create.data.non.proband(
  data = data,
  indices.non.proband = indices.non.proband,
  kin = kin,
  I = I
)

data.non.proband.delta.0 <- subset(
  data.non.proband,
  delta.non.proband == 0
)

data.non.proband.delta.1 <- subset(
  data.non.proband,
  delta.non.proband == 1
)

data.non.proband.delta.0 <- drop_col_by_name(
  data.non.proband.delta.0,
  "delta.non.proband"
)

data.non.proband.delta.1 <- drop_col_by_name(
  data.non.proband.delta.1,
  "delta.non.proband"
)

Data_Biv <- create.data.bivarie(
  data = data,
  indices.bivariee = indices.bivariee,
  kinJK = kinJK,
  I = I
)

Data_Biv_00 <- subset(
  Data_Biv,
  delta.non.probandIJ == 0 & delta.non.probandIK == 0
)

Data_Biv_01 <- subset(
  Data_Biv,
  delta.non.probandIJ == 0 & delta.non.probandIK == 1
)

Data_Biv_10 <- subset(
  Data_Biv,
  delta.non.probandIJ == 1 & delta.non.probandIK == 0
)

Data_Biv_11 <- subset(
  Data_Biv,
  delta.non.probandIJ == 1 & delta.non.probandIK == 1
)

saveRDS(
  data.proband,
  file = "output/data/data_proband.rds"
)

saveRDS(
  data.non.proband,
  file = "output/data/data_non_proband.rds"
)

saveRDS(
  Data_Biv,
  file = "output/data/data_bivariate.rds"
)

############################################################
## 10. Function to compute survival quantities for variance
############################################################

compute_survival_blocks <- function(paramsb) {
  
  S1_00 <- survie(
    t = Data_Biv_00$Y.proband,
    paramsb,
    X1 = Data_Biv_00$X1.proband,
    X2 = Data_Biv_00$X2.proband
  )
  
  SJ_00 <- survie(
    t = Data_Biv_00$Y.non.probandIJ,
    paramsb,
    X1 = Data_Biv_00$X1.non.probandIJ,
    X2 = Data_Biv_00$X2.non.probandIJ
  )
  
  SK_00 <- survie(
    t = Data_Biv_00$Y.non.probandIK,
    paramsb,
    X1 = Data_Biv_00$X1.non.probandIK,
    X2 = Data_Biv_00$X2.non.probandIK
  )
  
  S1_01 <- survie(
    t = Data_Biv_01$Y.proband,
    paramsb,
    X1 = Data_Biv_01$X1.proband,
    X2 = Data_Biv_01$X2.proband
  )
  
  SJ_01 <- survie(
    t = Data_Biv_01$Y.non.probandIJ,
    paramsb,
    X1 = Data_Biv_01$X1.non.probandIJ,
    X2 = Data_Biv_01$X2.non.probandIJ
  )
  
  SK_01 <- survie(
    t = Data_Biv_01$Y.non.probandIK,
    paramsb,
    X1 = Data_Biv_01$X1.non.probandIK,
    X2 = Data_Biv_01$X2.non.probandIK
  )
  
  S1_10 <- survie(
    t = Data_Biv_10$Y.proband,
    paramsb,
    X1 = Data_Biv_10$X1.proband,
    X2 = Data_Biv_10$X2.proband
  )
  
  SJ_10 <- survie(
    t = Data_Biv_10$Y.non.probandIJ,
    paramsb,
    X1 = Data_Biv_10$X1.non.probandIJ,
    X2 = Data_Biv_10$X2.non.probandIJ
  )
  
  SK_10 <- survie(
    t = Data_Biv_10$Y.non.probandIK,
    paramsb,
    X1 = Data_Biv_10$X1.non.probandIK,
    X2 = Data_Biv_10$X2.non.probandIK
  )
  
  S1_11 <- survie(
    t = Data_Biv_11$Y.proband,
    paramsb,
    X1 = Data_Biv_11$X1.proband,
    X2 = Data_Biv_11$X2.proband
  )
  
  SJ_11 <- survie(
    t = Data_Biv_11$Y.non.probandIJ,
    paramsb,
    X1 = Data_Biv_11$X1.non.probandIJ,
    X2 = Data_Biv_11$X2.non.probandIJ
  )
  
  SK_11 <- survie(
    t = Data_Biv_11$Y.non.probandIK,
    paramsb,
    X1 = Data_Biv_11$X1.non.probandIK,
    X2 = Data_Biv_11$X2.non.probandIK
  )
  
  return(
    list(
      S1_00 = S1_00,
      SJ_00 = SJ_00,
      SK_00 = SK_00,
      S1_01 = S1_01,
      SJ_01 = SJ_01,
      SK_01 = SK_01,
      S1_10 = S1_10,
      SJ_10 = SJ_10,
      SK_10 = SK_10,
      S1_11 = S1_11,
      SJ_11 = SJ_11,
      SK_11 = SK_11
    )
  )
}

############################################################
## 11. Function to estimate model and variance for one gamma
############################################################

estimate_for_gamma <- function(gamma_value) {
  
  cat("Estimating model for gamma =", gamma_value, "\n")
  
  obj <- compute.mle(
    params = params0,
    h = h,
    data.proband = data.proband,
    data.non.proband.delta.0 = data.non.proband.delta.0,
    data.non.proband.delta.1 = data.non.proband.delta.1,
    Data_Biv_00 = Data_Biv_00,
    Data_Biv_01 = Data_Biv_01,
    Data_Biv_10 = Data_Biv_10,
    Data_Biv_11 = Data_Biv_11,
    gamma = gamma_value,
    P = P1
  )
  
  par_hat <- extract_par(obj)
  
  paramsb <- par_hat[1:13]
  
  survival_blocks <- compute_survival_blocks(paramsb)
  
  var_hat <- Var_asymp(
    Theta = par_hat,
    data.proband = data.proband,
    data.non.proband.delta.0 = data.non.proband.delta.0,
    data.non.proband.delta.1 = data.non.proband.delta.1,
    Data_Biv_00 = Data_Biv_00,
    Data_Biv_01 = Data_Biv_01,
    Data_Biv_10 = Data_Biv_10,
    Data_Biv_11 = Data_Biv_11,
    S1_00 = survival_blocks$S1_00,
    SJ_00 = survival_blocks$SJ_00,
    SK_00 = survival_blocks$SK_00,
    S1_01 = survival_blocks$S1_01,
    SJ_01 = survival_blocks$SJ_01,
    SK_01 = survival_blocks$SK_01,
    S1_10 = survival_blocks$S1_10,
    SJ_10 = survival_blocks$SJ_10,
    SK_10 = survival_blocks$SK_10,
    S1_11 = survival_blocks$S1_11,
    SJ_11 = survival_blocks$SJ_11,
    SK_11 = survival_blocks$SK_11,
    gamma = gamma_value,
    P = P1
  )
  
  return(
    list(
      gamma = gamma_value,
      par_hat = par_hat,
      var_hat = var_hat
    )
  )
}

############################################################
## 12. Estimate model for gamma = 0, 0.05, 0.2, 0.8
############################################################

gamma_values <- c(0, 0.05, 0.2, 0.8)

results_list <- vector("list", length(gamma_values))

names(results_list) <- paste0(
  "gamma_",
  vapply(gamma_values, safe_gamma_label, character(1))
)

for (i in seq_along(gamma_values)) {
  
  g <- gamma_values[i]
  g_label <- safe_gamma_label(g)
  
  fit_g <- estimate_for_gamma(gamma_value = g)
  
  results_list[[i]] <- fit_g
  
  saveRDS(
    fit_g$par_hat,
    file = paste0(
      "output/estimates/parameter_estimates_gamma_",
      g_label,
      ".rds"
    )
  )
  
  saveRDS(
    fit_g$var_hat,
    file = paste0(
      "output/variances/asymptotic_variance_gamma_",
      g_label,
      ".rds"
    )
  )
}

############################################################
## 13. Save combined results
############################################################

saveRDS(
  results_list,
  file = "output/all_results_by_gamma.rds"
)

parameter_estimates <- do.call(
  rbind,
  lapply(results_list, function(x) x$par_hat)
)

rownames(parameter_estimates) <- names(results_list)

saveRDS(
  parameter_estimates,
  file = "output/estimates/all_parameter_estimates.rds"
)

############################################################
## 14. Save session information
############################################################

writeLines(
  capture.output(sessionInfo()),
  con = "sessionInfo.txt"
)

############################################################
## 15. End
############################################################

cat("\nReproducible simulation completed successfully.\n")
cat("Results are saved in the output/ folder.\n")
cat("Session information saved in sessionInfo.txt.\n")
