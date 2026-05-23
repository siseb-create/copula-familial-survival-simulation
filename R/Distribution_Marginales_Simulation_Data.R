############################################################
## Marginal functions for data simulation
##
## This file contains the marginal distribution functions required for the synthetic data
## generation procedure used in the simulation study.
############################################################

# Pré-calcul une seule fois au début du script


# Pré-calcul GL(50) sur [0,1]

GL50 <- gauss.quad(50, kind = "legendre")
GL50_nodes   <- 0.5 * GL50$nodes + 0.5      # transforme [-1,1] → [0,1]
GL50_weights <- 0.5 * GL50$weights          # ajuste les poids

# forme paramétrique de la fonction B-spline utilisée

g <- function(t) {  0.55 +  0.35 * sin(2*pi*t/1.8) * exp(-0.35*t) + 0.15 * (t/4)}

# ──────────────────────────────
# Cumul.H() finale — intégrale de 0 à t, avec spline active seulement après st1
# ──────────────────────────────


Taux.h_param <- function(t, theta, X1, X2) {
  alpha  <- exp(theta[1])
  lambda <- exp(theta[2])
  beta_g <- theta[3]
  
  base <- alpha * lambda * t^(alpha - 1) * exp(beta_g * X1)
  
  eta <- numeric(length(t))
  idx <- t > X2
  if (any(idx)) {
    eta[idx] <- g(t[idx] - X2[idx])
  }
  
  base * exp(eta)
}


# =========================
# 1) Fabrique d'intégrateur GL composite (pré-calcul des nœuds/poids)
# =========================

Taux.H_param <- function(t, theta, X1, X2) {
  alpha  <- exp(theta[1])
  lambda <- exp(theta[2])
  beta_g <- theta[3]
  
  n <- length(t)
  H <- numeric(n)
  
  pos <- t > 0
  if (!any(pos)) return(H)
  
  # ---- Partie analytique : [0, min(t, X2)]
  H[pos] <- lambda * pmin(t[pos], X2[pos])^alpha * exp(beta_g * X1[pos])
  
  # ---- Partie spline : [X2, t]
  masque <- (t > X2) & (t > 0)
  if (!any(masque)) return(H)
  
  a <- pmax(X2[masque], 0)
  b <- t[masque]
  width <- b - a
  if (any(width <= 0)) return(H)
  
  # Nœuds GL
  u_mat <- a + tcrossprod(width, GL50_nodes)  # (n_mask x 50)
  
  # âge post-événement
  age_post <- u_mat - X2[masque]
  age_post[age_post < 0] <- 0
  
  # exp(g(u - X2))
  exp_eta_mat <- exp(g(age_post))
  
  # facteur Weibull
  base_mat <- alpha * lambda * (u_mat^(alpha - 1)) *
    exp(beta_g * X1[masque])
  
  # intégrale GL
  integr <- width * as.vector((base_mat * exp_eta_mat) %*% GL50_weights)
  
  H[masque] <- H[masque] + integr
  H
}




survie_param <- function(t, theta, X1, X2,res=0) {exp(-Taux.H_param(t = t, theta = theta, X1 = X1, X2 = X2))-res}

Densite_param <- function(t, theta, X1, X2) {Taux.h_param(t = t, theta = theta,  X1 = X1, X2 = X2) * survie_param(t = t, theta = theta, X1 = X1, X2 = X2)}

 

