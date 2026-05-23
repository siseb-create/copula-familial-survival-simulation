Generate.X1 <- function(p) {sample(x=c(-1,1),size=6,replace=TRUE,prob=c(1-p,p))}

Generate.X2 <- function() {runif(6,min=0,max=4)}

Generate.a <- function(m,v)
{
ab <- calculer.a.b.from.m.v(m,v)
shape <- ab$a
scale <- ab$b
rweibull(1,shape=shape,scale=scale)
} 

Generate.C <- function(min,max) {runif(5,min=min,max=max)} 
