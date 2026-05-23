
compute.mle<-function(params, h,data.proband,data.non.proband.delta.0,data.non.proband.delta.1,Data_Biv_00,Data_Biv_01,Data_Biv_10,Data_Biv_11,gamma,P){
  

 repeat {
    
 Theta0 <- c(params,h)
 
 R <- optim(params,NegLogLik.Univariee,h=h,
                            data.proband=data.proband,data.non.proband.delta.0=data.non.proband.delta.0,
                            data.non.proband.delta.1=data.non.proband.delta.1,
            gamma=gamma,P=P)$par

 # Function de survive
 
 S1_00=survie(t=Data_Biv_00$Y.proband,params,X1=Data_Biv_00$X1.proband,X2=Data_Biv_00$X2.proband)
 SJ_00=survie(t=Data_Biv_00$Y.non.probandIJ,params,X1=Data_Biv_00$X1.non.probandIJ,X2=Data_Biv_00$X2.non.probandIJ)
 SK_00=survie(t=Data_Biv_00$Y.non.probandIK,params,X1=Data_Biv_00$X1.non.probandIK,X2=Data_Biv_00$X2.non.probandIK)
 
 
 S1_01=survie(t=Data_Biv_01$Y.proband,params,X1=Data_Biv_01$X1.proband,X2=Data_Biv_01$X2.proband)
 SJ_01=survie(t=Data_Biv_01$Y.non.probandIJ,params,X1=Data_Biv_01$X1.non.probandIJ,X2=Data_Biv_01$X2.non.probandIJ)
 SK_01=survie(t=Data_Biv_01$Y.non.probandIK,params,X1=Data_Biv_01$X1.non.probandIK,X2=Data_Biv_01$X2.non.probandIK)
 
 
 S1_10=survie(t=Data_Biv_10$Y.proband,params,X1=Data_Biv_10$X1.proband,X2=Data_Biv_10$X2.proband)
 SJ_10=survie(t=Data_Biv_10$Y.non.probandIJ,params,X1=Data_Biv_10$X1.non.probandIJ,X2=Data_Biv_10$X2.non.probandIJ)
 SK_10=survie(t=Data_Biv_10$Y.non.probandIK,params,X1=Data_Biv_10$X1.non.probandIK,X2=Data_Biv_10$X2.non.probandIK)
 
 
 S1_11=survie(t=Data_Biv_11$Y.proband,params,X1=Data_Biv_11$X1.proband,X2=Data_Biv_11$X2.proband)
 SJ_11=survie(t=Data_Biv_11$Y.non.probandIJ,params,X1=Data_Biv_11$X1.non.probandIJ,X2=Data_Biv_11$X2.non.probandIJ)
 SK_11=survie(t=Data_Biv_11$Y.non.probandIK,params,X1=Data_Biv_11$X1.non.probandIK,X2=Data_Biv_11$X2.non.probandIK)
  
h.op<-optimise(f=NegLogLik.Biv,lower = 0,upper = 1,params=params,
                      Data_Biv_00=Data_Biv_00,Data_Biv_01=Data_Biv_01,Data_Biv_10=Data_Biv_10,Data_Biv_11=Data_Biv_11,
                      S1_00=S1_00,SJ_00=SJ_00,SK_00=SK_00,
                      S1_01=S1_01,SJ_01=SJ_01,SK_01=SK_01,
                      S1_10=S1_10,SJ_10=SJ_10,SK_10=SK_10,
                      S1_11=S1_11,SJ_11=SJ_11,SK_11=SK_11)$minimum
 Theta <- c(R,h.op)
    if (max(abs(Theta-Theta0))<1e-3) {
      break
    }
 params <-R
 h<-h.op
 }
 return(Theta)
}
