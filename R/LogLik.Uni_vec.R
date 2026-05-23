
NegLogLik.proband_vec <- function(params,data)
{

terme1 <- l.densite(data$Y,params,data$X1,data$X2)
terme2 <- log(1-survie(data$a,params,data$X1,data$X2))
cbind(-(terme1-terme2),data$ID)
}

NegLogLik.non.proband.delta.0_vec <- function(params,h,data)
{

S.proband <- survie(data$Y.proband,params,data$X1.proband,data$X2.proband)
S.non.proband <- survie(data$Y.non.proband,params,data$X1.non.proband,data$X2.non.proband)
cbind(-(log(BiCopHfunc2(u1=S.non.proband,u2=S.proband,family=1,par=h*data$kin))),data$ID)
}

NegLogLik.non.proband.delta.1_vec <- function(params,h,data)
{

S.proband <- survie(data$Y.proband,params,data$X1.proband,data$X2.proband)
S.non.proband <- survie(data$Y.non.proband,params,data$X1.non.proband,data$X2.non.proband)
terme1 <- -(log(BiCopPDF(u1=S.non.proband,u2=S.proband,family=1,par=h*data$kin)))
terme2 <- -(l.densite(data$Y.non.proband,params,data$X1.non.proband,data$X2.non.proband))
cbind(terme1+terme2,data$ID)
}

NegLogLik.Univariee_vec <- function(params,h,data.proband,data.non.proband.delta.0,data.non.proband.delta.1,gamma,P)
{
penalite <-   0.5*gamma*as.numeric(t(params[4:13]) %*% P %*% params[4:13])
vect<-rbind(NegLogLik.proband_vec(params,data.proband),
NegLogLik.non.proband.delta.0_vec(params,h,data.non.proband.delta.0),
NegLogLik.non.proband.delta.1_vec(params,h,data.non.proband.delta.1))
vec_ag<-aggregate(vect[,1]~vect[,2],FUN=sum)+penalite/200
return(vec_ag[,2])
}


