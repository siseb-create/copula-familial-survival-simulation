############################################################
## B simulations script
##
## This script runs B reproducible simulations.
## For each simulated data set, the proposed model is estimated
## for gamma = 0, 0.05, 0.2, and 0.8.
##
## All outputs are saved in the output_B/ folder.
############################################################

rm(list = ls())

############################################################
## 0. Required packages
############################################################

library(survival)
library(numDeriv)
library(Matrix)
library(fda)
library(VineCopula)
library(condMVNorm)
library(statmod)

############################################################
## 1. Source functions
############################################################

source("R/Distribution_Marginales_Simulation_Data.R")
source("R/Distribution_marginales_10.R")
source("R/Generate.Covariates.R")
source("R/Proba.R")
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
## 2. Create output_B folders
############################################################

dir.create("output_B", recursive = TRUE, showWarnings = FALSE)
dir.create("output_B/data", recursive = TRUE, showWarnings = FALSE)
dir.create("output_B/estimates", recursive = TRUE, showWarnings = FALSE)
dir.create("output_B/variances", recursive = TRUE, showWarnings = FALSE)
dir.create("output_B/summary", recursive = TRUE, showWarnings = FALSE)
dir.create("output_B/errors", recursive = TRUE, showWarnings = FALSE)

############################################################
## 3. Reproducibility seed
############################################################

set.seed(678910)

############################################################
## 4. Monte Carlo settings
############################################################

B <- 20

gamma_values <- c(0, 0.05, 0.2, 0.8)

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

gamma_labels <- vapply(
  gamma_values,
  safe_gamma_label,
  character(1)
)

names(gamma_values) <- paste0("gamma_", gamma_labels)

############################################################
## 5. Model parameters and simulation settings
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
## 6. Create indices
############################################################

indices.proband <- creer.indices.proband(I)

indices.non.proband <- creer.indices.non.proband(I)

indices.bivariee <- creer.indices.bivariee(I)

############################################################
## 7. Penalty matrix for P-splines
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
## 8. Utility functions
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

############################################################
## 9. Function to construct analysis data sets
############################################################

build_analysis_data <- function(data) {
  
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
  
  return(
    list(
      data.proband = data.proband,
      data.non.proband = data.non.proband,
      data.non.proband.delta.0 = data.non.proband.delta.0,
      data.non.proband.delta.1 = data.non.proband.delta.1,
      Data_Biv = Data_Biv,
      Data_Biv_00 = Data_Biv_00,
      Data_Biv_01 = Data_Biv_01,
      Data_Biv_10 = Data_Biv_10,
      Data_Biv_11 = Data_Biv_11
    )
  )
}

############################################################
## 10. Function to compute survival quantities for variance
############################################################

compute_survival_blocks <- function(paramsb, analysis_data) {
  
  Data_Biv_00 <- analysis_data$Data_Biv_00
  Data_Biv_01 <- analysis_data$Data_Biv_01
  Data_Biv_10 <- analysis_data$Data_Biv_10
  Data_Biv_11 <- analysis_data$Data_Biv_11
  
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

estimate_for_gamma <- function(gamma_value, analysis_data) {
  
  obj <- compute.mle(
    params = params0,
    h = h,
    data.proband = analysis_data$data.proband,
    data.non.proband.delta.0 = analysis_data$data.non.proband.delta.0,
    data.non.proband.delta.1 = analysis_data$data.non.proband.delta.1,
    Data_Biv_00 = analysis_data$Data_Biv_00,
    Data_Biv_01 = analysis_data$Data_Biv_01,
    Data_Biv_10 = analysis_data$Data_Biv_10,
    Data_Biv_11 = analysis_data$Data_Biv_11,
    gamma = gamma_value,
    P = P1
  )
  
  par_hat <- extract_par(obj)
  
  paramsb <- par_hat[1:13]
  
  survival_blocks <- compute_survival_blocks(
    paramsb = paramsb,
    analysis_data = analysis_data
  )
  
  vcov_hat <- Var_asymp(
    Theta = par_hat,
    data.proband = analysis_data$data.proband,
    data.non.proband.delta.0 = analysis_data$data.non.proband.delta.0,
    data.non.proband.delta.1 = analysis_data$data.non.proband.delta.1,
    Data_Biv_00 = analysis_data$Data_Biv_00,
    Data_Biv_01 = analysis_data$Data_Biv_01,
    Data_Biv_10 = analysis_data$Data_Biv_10,
    Data_Biv_11 = analysis_data$Data_Biv_11,
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
  
  variance_diag <- diag(vcov_hat)
  se_hat <- sqrt(variance_diag)
  
  return(
    list(
      gamma = gamma_value,
      par_hat = par_hat,
      vcov_hat = vcov_hat,
      variance_diag = variance_diag,
      se_hat = se_hat
    )
  )
}

############################################################
## 12. Monte Carlo simulation loop
############################################################

all_results <- vector("list", B)
error_log <- vector("list", B)

for (b in seq_len(B)) {
  
  cat("\n====================================\n")
  cat("Simulation b =", b, "of", B, "\n")
  cat("====================================\n")
  
  error_log[[b]] <- character(0)
  
  simulated_data <- tryCatch(
    Generate.data(
      I = I,
      theta = params,
      m.a = m.a,
      v.a = v.a,
      min.c = min.c,
      max.c = max.c,
      p = p,
      h = h,
      kinship = kinship
    ),
    error = function(e) {
      error_log[[b]] <<- c(
        error_log[[b]],
        paste0("Generate.data: ", e$message)
      )
      return(NULL)
    }
  )
  
  if (is.null(simulated_data)) {
    next
  }
  
  saveRDS(
    simulated_data,
    file = paste0("output_B/data/simulated_data_b", b, ".rds")
  )
  
  analysis_data <- tryCatch(
    build_analysis_data(simulated_data),
    error = function(e) {
      error_log[[b]] <<- c(
        error_log[[b]],
        paste0("build_analysis_data: ", e$message)
      )
      return(NULL)
    }
  )
  
  if (is.null(analysis_data)) {
    next
  }
  
  all_results[[b]] <- vector("list", length(gamma_values))
  names(all_results[[b]]) <- names(gamma_values)
  
  for (g_name in names(gamma_values)) {
    
    g <- gamma_values[[g_name]]
    g_label <- safe_gamma_label(g)
    
    cat("Estimating model for", g_name, "\n")
    
    fit_g <- tryCatch(
      estimate_for_gamma(
        gamma_value = g,
        analysis_data = analysis_data
      ),
      error = function(e) {
        error_log[[b]] <<- c(
          error_log[[b]],
          paste0("gamma = ", g, ": ", e$message)
        )
        return(NULL)
      }
    )
    
    if (is.null(fit_g)) {
      next
    }
    
    all_results[[b]][[g_name]] <- fit_g
    
    saveRDS(
      fit_g$par_hat,
      file = paste0(
        "output_B/estimates/parameter_estimates_b",
        b,
        "_gamma_",
        g_label,
        ".rds"
      )
    )
    
    saveRDS(
      fit_g$vcov_hat,
      file = paste0(
        "output_B/variances/asymptotic_variance_covariance_b",
        b,
        "_gamma_",
        g_label,
        ".rds"
      )
    )
    
    saveRDS(
      fit_g$variance_diag,
      file = paste0(
        "output_B/variances/asymptotic_variance_b",
        b,
        "_gamma_",
        g_label,
        ".rds"
      )
    )
    
    saveRDS(
      fit_g$se_hat,
      file = paste0(
        "output_B/variances/standard_errors_b",
        b,
        "_gamma_",
        g_label,
        ".rds"
      )
    )
  }
  
  saveRDS(
    all_results,
    file = "output_B/all_results_by_replication.rds"
  )
  
  saveRDS(
    error_log,
    file = "output_B/errors/error_log.rds"
  )
}

############################################################
## 13. Build parameter matrices by gamma
############################################################

parameter_estimates_by_gamma <- vector(
  "list",
  length(gamma_values)
)

names(parameter_estimates_by_gamma) <- names(gamma_values)

for (g_name in names(gamma_values)) {
  
  par_list_g <- lapply(
    all_results,
    function(x) {
      if (is.null(x)) return(NULL)
      if (is.null(x[[g_name]])) return(NULL)
      x[[g_name]]$par_hat
    }
  )
  
  par_list_g <- Filter(Negate(is.null), par_list_g)
  
  if (length(par_list_g) > 0) {
    parameter_estimates_by_gamma[[g_name]] <- do.call(
      rbind,
      par_list_g
    )
  } else {
    parameter_estimates_by_gamma[[g_name]] <- NULL
  }
}

saveRDS(
  parameter_estimates_by_gamma,
  file = "output_B/summary/parameter_estimates_by_gamma.rds"
)

############################################################
## 14. Mean parameter estimates by gamma
############################################################

mean_parameter_estimates_by_gamma <- lapply(
  parameter_estimates_by_gamma,
  function(mat) {
    if (is.null(mat)) return(NULL)
    colMeans(mat, na.rm = TRUE)
  }
)

saveRDS(
  mean_parameter_estimates_by_gamma,
  file = "output_B/summary/mean_parameter_estimates_by_gamma.rds"
)

############################################################
## 15. Build variance diagonal matrices by gamma
############################################################

variance_diag_by_gamma <- vector(
  "list",
  length(gamma_values)
)

names(variance_diag_by_gamma) <- names(gamma_values)

for (g_name in names(gamma_values)) {
  
  var_list_g <- lapply(
    all_results,
    function(x) {
      if (is.null(x)) return(NULL)
      if (is.null(x[[g_name]])) return(NULL)
      x[[g_name]]$variance_diag
    }
  )
  
  var_list_g <- Filter(Negate(is.null), var_list_g)
  
  if (length(var_list_g) > 0) {
    variance_diag_by_gamma[[g_name]] <- do.call(
      rbind,
      var_list_g
    )
  } else {
    variance_diag_by_gamma[[g_name]] <- NULL
  }
}

saveRDS(
  variance_diag_by_gamma,
  file = "output_B/summary/variance_diag_by_gamma.rds"
)

############################################################
## 16. Mean variance diagonals by gamma
############################################################

mean_variance_diag_by_gamma <- lapply(
  variance_diag_by_gamma,
  function(mat) {
    if (is.null(mat)) return(NULL)
    colMeans(mat, na.rm = TRUE)
  }
)

saveRDS(
  mean_variance_diag_by_gamma,
  file = "output_B/summary/mean_variance_diag_by_gamma.rds"
)

############################################################
## 17. Standard errors by gamma
############################################################

standard_errors_by_gamma <- lapply(
  variance_diag_by_gamma,
  function(mat) {
    if (is.null(mat)) return(NULL)
    sqrt(mat)
  }
)

saveRDS(
  standard_errors_by_gamma,
  file = "output_B/summary/standard_errors_by_gamma.rds"
)

############################################################
## 18. Mean standard errors by gamma
############################################################

mean_standard_errors_by_gamma <- lapply(
  standard_errors_by_gamma,
  function(mat) {
    if (is.null(mat)) return(NULL)
    colMeans(mat, na.rm = TRUE)
  }
)

saveRDS(
  mean_standard_errors_by_gamma,
  file = "output_B/summary/mean_standard_errors_by_gamma.rds"
)

############################################################
## 19. Build summary tables by gamma
############################################################

summary_tables_by_gamma <- vector(
  "list",
  length(gamma_values)
)

names(summary_tables_by_gamma) <- names(gamma_values)

for (g_name in names(gamma_values)) {
  
  mean_par <- mean_parameter_estimates_by_gamma[[g_name]]
  mean_var <- mean_variance_diag_by_gamma[[g_name]]
  mean_se <- mean_standard_errors_by_gamma[[g_name]]
  
  if (is.null(mean_par)) {
    summary_tables_by_gamma[[g_name]] <- NULL
    next
  }
  
  summary_tables_by_gamma[[g_name]] <- data.frame(
    parameter = seq_along(mean_par),
    mean_estimate = as.numeric(mean_par),
    mean_variance_diag = as.numeric(mean_var),
    mean_standard_error = as.numeric(mean_se)
  )
}

saveRDS(
  summary_tables_by_gamma,
  file = "output_B/summary/summary_tables_by_gamma.rds"
)

for (g_name in names(summary_tables_by_gamma)) {
  
  tab_g <- summary_tables_by_gamma[[g_name]]
  
  if (is.null(tab_g)) {
    next
  }
  
  write.csv(
    tab_g,
    file = paste0("output_B/summary/summary_table_", g_name, ".csv"),
    row.names = FALSE
  )
}

############################################################
## 20. Save session information
############################################################

writeLines(
  capture.output(sessionInfo()),
  con = "output_B/sessionInfo_B.txt"
)

############################################################
## 21. End
############################################################

cat("\nB-simulation script completed successfully.\n")
cat("Number of requested simulations:", B, "\n")
cat("Results are saved in the output_B/ folder.\n")
cat("Session information saved in output_B/sessionInfo_B.txt.\n")
