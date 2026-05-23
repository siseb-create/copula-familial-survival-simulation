############################################################
## Three-dimensional copula functions and derivatives
##
## This file contains the functions used to evaluate
## three-dimensional copulas and their partial derivatives.
## These quantities are used to model intra-family dependence
## in the proposed survival model.
############################################################

C011.i <- function(i,u1,u2,u3,r12,r13,r23)
{
sigma = matrix(c(1,r12[i],r13[i],r12[i],1,r23[i],r13[i],r23[i],1),ncol=3)
mean = rep(0,3)
pcmvnorm(upper=qnorm(u1[i]),mean=mean,sigma=sigma,dependent.ind=1,given.ind=c(2,3),X.given=c(qnorm(u2[i]),qnorm(u3[i])))
}

C011 <- function(u1,u2,u3,r12,r13,r23)
{
len <- length(u1)
sapply(1:len,C011.i,u1=u1,u2=u2,u3=u3,r12=r12,r13=r13,r23=r23)
}

log.C011.bis <- function(u1,u2,u3,r12,r13,r23)
{
Z1 <- qnorm(u1)
Z2 <- qnorm(u2)
Z3 <- qnorm(u3)
D <- 1-((r23)^2)
Moyenne <- (Z2*(r12-r13*r23)+Z3*(r13-r12*r23))/D
Variance <- 1-(r12^2+r13^2-2*r12*r13*r23)/D
Ecart.Type <- sqrt(Variance)
pnorm(Z1,mean=Moyenne,sd=Ecart.Type,log.p=TRUE)
}

#########################################################

C001.i <- function(i,u1,u2,u3,r12,r13,r23)
{
mean <- rep(0,3)
sigma <- matrix(c(1,r12[i],r13[i],r12[i],1,r23[i],r13[i],r23[i],1),ncol=3)
pcmvnorm(upper=c(qnorm(u1[i]),qnorm(u2[i])),mean=mean,sigma=sigma,dependent.ind=c(1,2),given.ind=3,X.given=qnorm(u3[i]))
}

C001 <- function(u1,u2,u3,r12,r13,r23)
{
len <- length(u1)
sapply(1:len,C001.i,u1=u1,u2=u2,u3=u3,r12=r12,r13=r13,r23=r23)
}

C001.i.bis <- function(i,Z1,Z2,M1,M2,V11,V12,V22)
{
mean <- c(M1[i],M2[i])
sigma <- matrix(c(V11[i],V12[i],V12[i],V22[i]),ncol=2)
pmvnorm(upper=c(Z1[i],Z2[i]),mean=mean,sigma=sigma)
}

C001.bis <- function(u1,u2,u3,r12,r13,r23)
{
len <- length(u1)
Z1 <- qnorm(u1)
Z2 <- qnorm(u2)
Z3 <- qnorm(u3)
M1 <- Z3*r13
M2 <- Z3*r23
V11 <- 1-r13^2
V12 <- r12-r13*r23
V22 <- 1-r23^2
sapply(1:len,C001.i.bis,Z1=Z1,Z2=Z2,M1=M1,M2=M2,V11=V11,V12=V12,V22=V22)
}

###########################################

log.C111.i <- function(i,r12, r13, r23, u1, u2, u3) {
  
  # Matrice de corrélation 3x3
  Sigma <- matrix(c(1,r12[i],r13[i],
                    r12[i],1,r23[i],
                    r13[i],r23[i], 1), 
                  nrow = 3, byrow = TRUE)
  
  # Calcul de l'inverse et du déterminant de Sigma
  Sigma_inv <- solve(Sigma)
  det_Sigma <- det(Sigma)
  
  # Transformation des u en quantiles gaussiens
  z <- qnorm(c(u1[i], u2[i], u3[i]))
  
  # Calcul de la densité de la copule gaussienne
##  densite <- (1/sqrt(det_Sigma))*exp(-0.5*(t(z)%*%(Sigma_inv-diag(3))%*%z))
  
log.densite <- -0.5*(t(z)%*%(Sigma_inv-diag(3))%*%z) - 0.5*log(det_Sigma)

  return(as.numeric(log.densite))
}

log.C111<- function(r12, r13, r23, u1, u2, u3) 
{
Densite3D<-sapply(1:length(u1),log.C111.i,r12, r13, r23, u1, u2, u3)
return(Densite3D)
}
#densite3 <- C111(r12=h*r12, r13=h*r13, r23=h*r23, u1=u1, u2=u2, u3=u3)


