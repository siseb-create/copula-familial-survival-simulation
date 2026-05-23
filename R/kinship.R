kinship <- matrix(
  c(1,1/2,1/8,1/8,1/8,1/8,
    1/2,1,1/8,1/8,1/8,1/8,
    1/8,1/8,1,1/2,1/8,1/8,
    1/8,1/8,1/2,1,1/8,1/8,
    1/8,1/8,1/8,1/8,1,1/2,
    1/8,1/8,1/8,1/8,1/2,1),
  ncol=6)
kin <- kinship[1,2:6]

kinJK <- cbind(kinship[1, combn(2:6, 2)[1, ]],
               kinship[1, combn(2:6, 2)[2, ]],
               mapply(function(j, k) kinship[j, k], combn(2:6, 2)[1, ], combn(2:6, 2)[2, ])
)


#data.simul2<-creer.data.simul_famille_sujet(data,familles_resample) 

# Extraire les coefficients de parenté (kinship)
extract_kinship <- function(kinship, data.liste) {
  # Créer une liste d'éléments pour chaque famille
  elements <- lapply(data.liste, function(famille) {
    data.frame(
      proband = 1,
      non.proband = famille$sujet[-1]
    )
  })
  
  # Appliquer pour chaque famille et chaque non-proband l'extraction du coefficient de parenté
  kinship_values <- unlist(lapply(elements, function(element) {
    sapply(element$non.proband, function(s) kinship[1, s])
  }))
  
  # Retourner le vecteur de kinship concaténé
  #names(kinship_values)<-NULL
  return(kinship_values)
}



extract_kinship_biv <- function(kinship, data.liste) {
  
  # Cette fonction interne extrait la matrice 3 colonnes pour une famille donnée
  extraire_pour_une_famille <- function(famille) {
    sujets <- famille$sujet
    non_proband <- sujets[-1]  # On retire le proband (premier sujet)
    
    # Toutes les combinaisons (j, k) entre non-probands
    paires <- t(combn(non_proband, 2))
    
    # Extraire les valeurs de parenté
    K_1j <- kinship[1, paires[, 1]]  # Parenté entre le proband (1) et j
    K_1k <- kinship[1, paires[, 2]]  # Parenté entre le proband (1) et k
    K_jk <- mapply(function(j, k) kinship[j, k], paires[, 1], paires[, 2])  # Parenté entre j et k
    
    # Retourne une matrice 3 colonnes
    cbind(K_1j, K_1k, K_jk)
  }
  
  # Appliquer à chaque famille de data.liste
  kinship_matrices <- lapply(data.liste, extraire_pour_une_famille)
  
  # Combiner toutes les matrices des familles en une seule
  kinship_matrix_finale <- do.call(rbind, kinship_matrices)
  
  return(kinship_matrix_finale)
}
 

