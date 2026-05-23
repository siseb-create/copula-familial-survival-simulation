############################################################
## Conditional and data generation functions
##
## This file contains the conditional functions used in the
## data-generating mechanism, together with the procedures
## required to simulate synthetic family-structured survival
## data. These functions are used to generate event times
## conditionally on the family structure and model parameters.
############################################################

survie.cond <- function(t,a,theta,X1,X2,res=0) 
{
S.t <- survie_param(t,theta,X1,X2)
S.a <- survie_param(a,theta,X1,X2)
(S.t-S.a)/(1-S.a)-res
}

Generate.survival.proband <- function(a,theta,X1.proband,X2.proband)
{
T.proband <- uniroot(survie.cond,lower=0,upper=a,a=a,
                     theta=theta,
                    X1=X1.proband,X2=X2.proband,res=runif(1))$root
U.proband <- survie_param(T.proband,theta,X1.proband,X2.proband)
Z.proband <- qnorm(U.proband)
return(data.frame(T.proband=T.proband,Z.proband=Z.proband))
}

Generate.survival.non.proband <- function(Z.proband,h,kinship,theta,X1.non.proband,X2.non.proband)
{
Sigma <- h*kinship + (1-h)*diag(rep(1,6))
condDist <- condMVN(mean=rep(0,6),sigma= Sigma,dependent.ind=(2:6),given.ind=1,X.given=Z.proband)
Z.non.proband <- rmvnorm(1,mean=condDist$condMean,sigma=condDist$condVar)
U.non.proband <- pnorm(Z.non.proband)
T.non.proband <- rep(-999,5)
for(i in 1:5)
{
T.non.proband[i] <- uniroot(survie_param,lower=0,upper=20,theta=theta,
                            X1=X1.non.proband[i],X2=X2.non.proband[i],res=U.non.proband[i])$root
}
T.non.proband
}

Generate.data.family <- function(theta,m.a,v.a,min.c,max.c,p,h,kinship)
{
X1 <- Generate.X1(p)
X1.proband <- X1[1]
X1.non.proband <- X1[-1]
X2 <- Generate.X2()
X2.proband <- X2[1]
X2.non.proband <- X2[-1]
a <- Generate.a(m.a,v.a)
C <- c(a,Generate.C(min.c,max.c))

survival.proband <- Generate.survival.proband(a,theta,X1.proband,X2.proband)
T.proband <- survival.proband$T.proband
Z.proband <- survival.proband$Z.proband
T.non.proband <- Generate.survival.non.proband(Z.proband,h,kinship,
                                               theta,X1.non.proband,X2.non.proband)
T <- c(T.proband,T.non.proband)
Y <- pmin(T,C)
delta <- as.numeric(T<C)
res <- data.frame(a=a,Y=Y,delta=delta,X1=X1,X2=X2)
res
}

Generate.data <- function(I,theta,m.a,v.a,min.c,max.c,p,h,kinship)
{
data <- NULL
for(i in 1:I)
{
data.family <- Generate.data.family(theta,m.a,v.a,min.c,max.c,p,h,kinship)
data.family$sujet<-1:6
data <- rbind(data,data.family)
}
#data$sujet<-1:(6*I)
data$ID<-rep(1:I, each = 6)
data
}

