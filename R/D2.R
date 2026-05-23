Creer.m<-function(Theta,data.proband,data.non.proband.delta.0,data.non.proband.delta.1,
                  Data_Biv_00,Data_Biv_01,Data_Biv_10,Data_Biv_11,S1_00,SJ_00,SK_00,S1_01,SJ_01,SK_01,S1_10,SJ_10,SK_10,S1_11,SJ_11,SK_11,gamma,P){
  
  m1<-jacobian(NegLogLik.Univariee,
               x=Theta[1:13],h=Theta[14],
               data.proband=data.proband,data.non.proband.delta.0=data.non.proband.delta.0,
               data.non.proband.delta.1=data.non.proband.delta.1,gamma=gamma,P=P)
  m2<-jacobian(NegLogLik.Biv,
               x=Theta[14],params = Theta[1:13],
               Data_Biv_00=Data_Biv_00,Data_Biv_01=Data_Biv_01,Data_Biv_10=Data_Biv_10,Data_Biv_11=Data_Biv_11,
               S1_00=S1_00,SJ_00=SJ_00,SK_00=SK_00,
               S1_01=S1_01,SJ_01=SJ_01,SK_01=SK_01,
               S1_10=S1_10,SJ_10=SJ_10,SK_10=SK_10,
               S1_11=S1_11,SJ_11=SJ_11,SK_11=SK_11)
  m=cbind(m1,m2)
  m=t(m)
  return(as.vector(m))
}

Var_asymp<-function(Theta,data.proband,data.non.proband.delta.0,data.non.proband.delta.1,
                    Data_Biv_00,Data_Biv_01,Data_Biv_10,Data_Biv_11,S1_00,SJ_00,SK_00,S1_01,SJ_01,SK_01,S1_10,SJ_10,SK_10,S1_11,SJ_11,SK_11,
                    gamma,P
)
{
  m1<-jacobian(NegLogLik.Univariee_vec,
               x=Theta[1:13],
               h=Theta[14],data.proband=data.proband,
               data.non.proband.delta.0=data.non.proband.delta.0,
               data.non.proband.delta.1=data.non.proband.delta.1,
               gamma=gamma,P=P)
  
  ## 2.2  Dérivée première de la vraisemblance bivariée (m_h)
  
  m2<-jacobian(NegLogLik.Biv_vec,
               x=Theta[14],
               params=Theta[1:13],
               Data_Biv_00=Data_Biv_00,Data_Biv_01=Data_Biv_01,Data_Biv_10=Data_Biv_10,Data_Biv_11=Data_Biv_11,
               S1_00=S1_00,SJ_00=SJ_00,SK_00=SK_00,
               S1_01=S1_01,SJ_01=SJ_01,SK_01=SK_01,
               S1_10=S1_10,SJ_10=SJ_10,SK_10=SK_10,
               S1_11=S1_11,SJ_11=SJ_11,SK_11=SK_11
               )
  
  ## 2.3 Expression de B chapeau
  
  m.Deriv1<-cbind(m1,m2)
  B.hat.list <-list()
  for (i in 1:I) {
    B.hat.list[[i]]<-m.Deriv1[i,]%*%t(m.Deriv1[i,])
  }
  
  B.hat <-matrix(0, nrow = 14,ncol = 14)
  for (j in 1:I) {
    B.hat<-B.hat+B.hat.list[[j]]
  }
  B.hat<-B.hat/I
  
  # 3 Calcul de A chapeau
  
  ## 3.1 Fonction  de la dérivée première (m_thet, m_h) afin d'appliquer à nouveau la fonction "jacobian"
  
  ## 3.2 fonction de la dérivée seconde (d.m_thet, d.m_h)
  
  m.Deriv<-jacobian(Creer.m,
                    x=Theta,
                    data.proband=data.proband,
                    data.non.proband.delta.0=data.non.proband.delta.0,
                    data.non.proband.delta.1=data.non.proband.delta.1,
                    Data_Biv_00=Data_Biv_00,Data_Biv_01=Data_Biv_01,Data_Biv_10=Data_Biv_10,Data_Biv_11=Data_Biv_11,
                    S1_00=S1_00,SJ_00=SJ_00,SK_00=SK_00,
                    S1_01=S1_01,SJ_01=SJ_01,SK_01=SK_01,
                    S1_10=S1_10,SJ_10=SJ_10,SK_10=SK_10,
                    S1_11=S1_11,SJ_11=SJ_11,SK_11=SK_11,
                    gamma=gamma,P=P)
  
  ## 3.3 Expression de A chapeau 
  
  A.chapeau<-m.Deriv/I
  
  # 4 expression de la variance
  
  var.hat<-(solve(A.chapeau)%*%B.hat%*%solve(A.chapeau))/I
  
  return(var.hat)
}


