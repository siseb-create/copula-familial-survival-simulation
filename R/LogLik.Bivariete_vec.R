NegLogLik.Biv_vec <- function(params,h,Data_Biv_00,Data_Biv_01,Data_Biv_10,Data_Biv_11,S1_00,SJ_00,SK_00,S1_01,SJ_01,SK_01,S1_10,SJ_10,SK_10,S1_11,SJ_11,SK_11)
{
Biva00<-cbind(-log(C001.bis(u1=SJ_00,u2=SK_00,u3=S1_00,r12=h*Data_Biv_00$kinJK,r13=h*Data_Biv_00$kin1J,r23=h*Data_Biv_00$kin1K)),Data_Biv_00$ID)
Biva01<-cbind(-log.C011.bis(u1=SJ_01,u2=SK_01,u3=S1_01,r12=h*Data_Biv_01$kinJK,r13=h*Data_Biv_01$kin1J,r23=h*Data_Biv_01$kin1K),Data_Biv_01$ID)
Biva10<-cbind(-log.C011.bis(u1=SK_10,u2=SJ_10,u3=S1_10,r12=h*Data_Biv_10$kinJK,r13=h*Data_Biv_10$kin1K,r23=h*Data_Biv_10$kin1J),Data_Biv_10$ID)
Biva11<-cbind(-log.C111(u1=SJ_11,u2=SK_11,u3=S1_11,r12=h*Data_Biv_11$kinJK,r13=h*Data_Biv_11$kin1J,r23=h*Data_Biv_11$kin1K),Data_Biv_11$ID)
Vect.mean.h <-aggregate(rbind(Biva00,Biva01,Biva10,Biva11)[,1]~rbind(Biva00,Biva01,Biva10,Biva11)[,2],FUN=sum)
return(Vect.mean.h[,2])
}

