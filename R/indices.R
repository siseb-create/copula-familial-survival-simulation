creer.indices.proband <- function(I)
{
seq(from=1,to=6*I-5,by=6)
}

creer.indices.non.proband <- function(I)
{
indices.non.proband <- NULL
for(i in 1:I)
{
for(j in 2:6)
{
indices.non.proband <- rbind(indices.non.proband,c(6*i-5,6*i+j-6))
}
}
indices.non.proband
}

creer.indices.bivariee<-function(I)
{
indices.bivariee <- NULL
for (i in 1:I) 
{
elements <- (6*i - 4):(6*i)
combs <- t(combn(elements, 2))
indices.bivariee <- rbind(indices.bivariee, cbind(6*i - 5, combs))
}
indices.bivariee
}


create.data.proband <- function(data,indices.proband) {data[indices.proband,-c(3,6)]}

create.data.non.proband <- function(data,indices.non.proband,kin,I)
{
data.non.proband <- cbind(data[indices.non.proband[,1],-c(1,3,6,7)],data[indices.non.proband[,2],-1])
data.non.proband <- cbind(data.non.proband,rep(kin,I))
colnames(data.non.proband) <- c("Y.proband","X1.proband","X2.proband","Y.non.proband","delta.non.proband","X1.non.proband","X2.non.proband","sujet","ID","kin")
data.non.proband
}

create.data.non.proband.s <- function(data,indices.non.proband,kin,I)
{
  data.non.proband <- cbind(data[indices.non.proband[,1],-c(1,3,6,7)],data[indices.non.proband[,2],-1])
  data.non.proband <- cbind(data.non.proband,kin)
  colnames(data.non.proband) <- c("Y.proband","X1.proband","X2.proband","Y.non.proband","delta.non.proband","X1.non.proband","X2.non.proband","sujet","ID","kin")
  data.non.proband
}


create.data.bivarie<-function(data, indices.bivariee,kinJK,I)
{
data.bivarie<- cbind(data[indices.bivariee[,1],c(-1,-3,-6,-7)],
                     data[indices.bivariee[,2], c(-1,-6,-7)],
                     data[indices.bivariee[,3], -c(1,6)],
                     do.call(rbind, replicate(I, kinJK, simplify = FALSE))
                     )
colnames(data.bivarie)<-c("Y.proband","X1.proband","X2.proband",
                          "Y.non.probandIJ","delta.non.probandIJ","X1.non.probandIJ","X2.non.probandIJ",
                          "Y.non.probandIK","delta.non.probandIK","X1.non.probandIK","X2.non.probandIK","ID",
                          "kin1J","kin1K","kinJK")
data.bivarie
}


creer.data.simul<-function(data, familles_resample)
{
data <- do.call(rbind, lapply(familles_resample, function(fam) {
data[data$ID == fam, ]
  }))
  return(data)
}


creer.data.simul_famille_sujet <- function(data, familles_resample) {
  
  data_resample<-NULL
  
  # Étape 2 - Pour chaque famille tirée, tirage avec remise des membres non-proband
  for (famille in familles_resample) {
    
    membres <- data[data$ID == famille, ]
    
    # On suppose que le proband est sujet 1, donc on garde seulement sujets 2 à 6
    membres_non_proband <- membres[2:6, ]
    
    # Tirage avec remise des membres non-proband (autant que le nombre original)
    sujets_resample <- sample(2:6, size = 5, replace = TRUE)
    
    # Créer la nouvelle famille rééchantillonnée
    
    data_famille_resample = data[data$ID==famille,][c(1, sujets_resample),]  # on garde le proband (Sujet 1) et les non-probands tirés
    
    # Ajouter à la base finale
    data_resample <- rbind(data_resample, data_famille_resample)
  }
  
  return(data_resample)
}


create.data.bivarie_simul<-function(kin_biv, data.liste) {
  
  # Cette fonction interne extrait la matrice 3 colonnes pour une famille donnée
  extraire_pour_une_famille <- function(famille) {
    sujets <- famille$sujet
    non_proband <- sujets[-1]  # On retire le proband (premier sujet)
    
    # Toutes les combinaisons (j, k) entre non-probands
    paires <- t(combn(non_proband, 2))
    
    # Extraire les valeurs de parenté
    d.proband<- famille[rep(1,10) ,c(-1,-3,-6,-7)]
    d.non.proband.j <- famille[paires[, 1],c(-1,-6,-7)] 
    d.non.proband.k <- famille[paires[, 2],c(-1,-6)]
    
    # Retourne une matrice 3 colonnes
    cbind(d.proband, d.non.proband.j, d.non.proband.k)
  }
  
  # Appliquer à chaque famille de data.liste
  data.bivarie_simul.l <- lapply(data.liste, extraire_pour_une_famille)
  
  # Combiner toutes les matrices des familles en une seule
  data.bivarie_simul<-cbind(do.call(rbind, data.bivarie_simul.l),kin_biv)
 
  colnames(data.bivarie_simul)<-c("Y.proband","X1.proband","X2.proband",
                            "Y.non.probandIJ","delta.non.probandIJ","X1.non.probandIJ","X2.non.probandIJ",
                            "Y.non.probandIK","delta.non.probandIK","X1.non.probandIK","X2.non.probandIK","ID",
                            "kin1J","kin1K","kinJK")
  
  return(data.bivarie_simul)
}



