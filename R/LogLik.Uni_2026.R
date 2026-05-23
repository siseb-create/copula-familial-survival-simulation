NegLogLik.proband <- function(params,data)
{

terme1 <- l.densite(data$Y,params,data$X1,data$X2)
terme2 <- log(1-survie(data$a,params,data$X1,data$X2))
-sum(terme1-terme2)
}

NegLogLik.non.proband.delta.0 <- function(params,h,data)
{

S.proband <- survie(data$Y.proband,params,data$X1.proband,data$X2.proband)
S.non.proband <- survie(data$Y.non.proband,params,data$X1.non.proband,data$X2.non.proband)
-sum(log(BiCopHfunc2(u1=S.non.proband,u2=S.proband,family=1,par=h*data$kin)))
}

NegLogLik.non.proband.delta.1 <- function(params,h,data)
{

S.proband <- survie(data$Y.proband,params,data$X1.proband,data$X2.proband)
S.non.proband <- survie(data$Y.non.proband,params,data$X1.non.proband,data$X2.non.proband)
terme1 <- -sum(log(BiCopPDF(u1=S.non.proband,u2=S.proband,family=1,par=h*data$kin)))
terme2 <- -sum(l.densite(data$Y.non.proband,params,data$X1.non.proband,data$X2.non.proband))
terme1+terme2
}

NegLogLik.Univariee <- function(params,h,data.proband,data.non.proband.delta.0,data.non.proband.delta.1,gamma,P)
{
penalite <- 0.5*gamma*as.numeric(t(params[4:13]) %*% P %*% params[4:13])
NegLogLik.proband (params,data.proband)+
NegLogLik.non.proband.delta.0(params,h,data.non.proband.delta.0)+
NegLogLik.non.proband.delta.1(params,h,data.non.proband.delta.1)+penalite
}
