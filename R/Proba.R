
ff <- function(a,valeur=0) {(gamma(1+2/a)/((gamma(1+1/a))**2))-valeur}

calculer.a.b.from.m.v <- function(m,v)
{
a <- uniroot(ff,lower=0.1,upper=1000,valeur=1+(v/(m^2)))$root
b <- m/gamma(1+1/a)
return(data.frame(a=a,b=b))
}


calculer.alpha.lambda.from.a.b <- function(a,b)
{
alpha <- a
lambda <- 1/(b^a)
return(data.frame(alpha=alpha,lambda=lambda))
}

calculer.m.v.from.a.b <-function(a,b)
{
m <- b*gamma(1+1/a)
v <- b^2*(gamma(1+2/a)-gamma(1+1/a)**2)
return(data.frame(m=m,v=v))
}

calculer.a.b.from.alpha.lambda <- function(alpha,lambda)
{
a <- alpha
b <- (1/lambda)^(1/a)
return(data.frame(a=a,b=b))
}


calculer.alpha.lambda.from.m.v <- function(m,v)
{
ab  <- calculer.a.b.from.m.v(m,v)
a <- ab$a
b <- ab$b
calculer.alpha.lambda.from.a.b(a,b)
}


calculer.m.v.from.alpha.lambda <- function(alpha,lambda)
{
ab <- calculer.a.b.from.alpha.lambda(alpha,lambda)
a <- ab$a
b <- ab$b
calculer.m.v.from.a.b(a,b)
}


