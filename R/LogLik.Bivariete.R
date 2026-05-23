############################################################
## Bivariate likelihood for intra-family dependence
##
## This file contains the bivariate likelihood contributions
## used to estimate the dependence parameter h. These likelihood
## components are based on pairs of family members and are used
## to capture intra-family dependence through the copula-based
## structure of the proposed survival model.
############################################################

NegLogLik.Biv <- function(params,h,Data_Biv_00,Data_Biv_01,Data_Biv_10,Data_Biv_11,S1_00,SJ_00,SK_00,S1_01,SJ_01,SK_01,S1_10,SJ_10,SK_10,S1_11,SJ_11,SK_11)
{
Biva00<-sum(C001.bis(u1=SJ_00,u2=SK_00,u3=S1_00,r12=h*Data_Biv_00$kinJK,r13=h*Data_Biv_00$kin1J,r23=h*Data_Biv_00$kin1K))
Biva01<-sum(log.C011.bis(u1=SJ_01,u2=SK_01,u3=S1_01,r12=h*Data_Biv_01$kinJK,r13=h*Data_Biv_01$kin1J,r23=h*Data_Biv_01$kin1K))
Biva10<-sum(log.C011.bis(u1=SK_10,u2=SJ_10,u3=S1_10,r12=h*Data_Biv_10$kinJK,r13=h*Data_Biv_10$kin1K,r23=h*Data_Biv_10$kin1J))
Biva11<-sum(log.C111(u1=SJ_11,u2=SK_11,u3=S1_11,r12=h*Data_Biv_11$kinJK,r13=h*Data_Biv_11$kin1J,r23=h*Data_Biv_11$kin1K))
return(-sum(log(Biva00),Biva01,Biva10,Biva11))
}

