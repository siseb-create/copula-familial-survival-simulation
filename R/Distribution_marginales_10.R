
#library(statmod)
GL50 <- gauss.quad(50, kind = "legendre")
GL50_nodes   <- 0.5 * GL50$nodes + 0.5
GL50_weights <- 0.5 * GL50$weights   # longueur 50

# Base spline + pénalité (une fois)

bspline_basis <- create.bspline.basis(c(0, 4), nbasis = 10, norder = 4)


Taux.h <- function(t, theta, X1, X2) {
  
  alpha  <- exp(theta[1])
  lambda <- exp(theta[2])
  beta_g <- theta[3]
  
  base <- alpha * lambda * t^(alpha - 1) * exp(beta_g * X1)
  
  eta <- numeric(length(t))
  idx <- t > X2
  if (any(idx)) {
    u <- pmax(0, pmin(t[idx] - X2[idx], 4))
    B <- eval.basis(u, bspline_basis)#[, 4:7, drop = FALSE]   # <-- essentiel
    eta[idx] <- drop(B %*% theta[4:13])
  }
  
  base * exp(eta)
}


Taux.H <- function(t, theta, X1, X2) {
  alpha  <- exp(theta[1])
  lambda <- exp(theta[2])
  beta_g <- theta[3]
  beta_true <- theta[4:13]  # <-- 4 coeffs cohérents avec 4:7
  
  n <- length(t)
  H <- numeric(n)
  
  pos <- t > 0
  if (!any(pos)) return(H)
  
  # Partie analytique: ∫0^{min(t,X2)} h0(u)du
  H[pos] <- lambda * pmin(t[pos], X2[pos])^alpha * exp(beta_g * X1[pos])
  
  masque <- (t > X2) & (t > 0)
  if (!any(masque)) return(H)
  
  a <- pmax(X2[masque], 0)
  b <- t[masque]
  width <- b - a
  if (any(width <= 0)) return(H)
  
  n_mask <- length(a)
  
  # u = a + width * node
  u_mat <- a + tcrossprod(width, GL50_nodes)   # (n_mask x 50)
  
  # âge post-événement
  age_post <- u_mat - X2[masque]
  age_post[age_post < 0] <- 0
  
  # clip au support de la base (ici supposé [0,65])
  age_vec <- as.numeric(age_post)
  age_vec <- pmax(0, pmin(age_vec, 4))
  
  # exp(eta) aux nœuds
  B_all <- eval.basis(age_vec, bspline_basis)#[, 4:7, drop = FALSE]  # (n_mask*50 x 4)
  exp_eta <- exp(drop(B_all %*% beta_true))                         # (n_mask*50)
  
  exp_mat <- matrix(exp_eta, nrow = n_mask, ncol = 50, byrow = FALSE)
  
  # facteur Weibull dans l'intégrande
  base_mat <- alpha * lambda * (u_mat^(alpha - 1)) * exp(beta_g * X1[masque])
  
  # intégrale GL
  integr <- width * as.vector((base_mat * exp_mat) %*% GL50_weights)
  
  H[masque] <- H[masque] + integr
  H
}

# Clamp dans (0,1) pour éviter 0 et 1 exacts
clamp01 <- function(u, eps = 1e-12) {
  pmin(pmax(u, eps), 1 - eps)
}

survie <- function(t, theta, X1, X2, eps = 1e-12,res=0) {
  # cumulative hazard
  H <- Taux.H(t = t, theta = theta, X1 = X1, X2 = X2)
  
  # survie brute (protection contre underflow)
  S_raw <- exp(-H)
  
  # clamp final
  S <- clamp01(S_raw, eps)
  
  S-res
}


densite <- function(t, theta, X1, X2) {
  Taux.h(t = t, theta = theta, X1 = X1, X2 = X2) *
    survie(t = t, theta = theta, X1 = X1, X2 = X2)
}

l.densite<-function(t, theta, X1, X2){
 ll<- Taux.h(t = t, theta = theta, X1 = X1, X2 = X2) *
    survie(t = t, theta = theta, X1 = X1, X2 = X2)
 return(log(ll)) 
}


