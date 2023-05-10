require(evd)

# log-Weibull distribution ------------------------------------------------
## probability density function ------------------------------------------------
dlogWeibull = function(x, location=0, scale=1, log=FALSE) {
  if (any(scale < 0)) 
    stop(paste("scale parameter must be positive", "\n", ""))
  log.lik <- -log(scale) + ((x - location)/scale) - exp((x - location)/scale)
  if (log == FALSE) 
    fy <- exp(log.lik)
  else fy <- log.lik
  fy
}
## cumulative distribution function ------------------------------------------------
plogWeibull = function(q, location=0, scale=1, lower.tail=TRUE, log.p=FALSE) {
  if (any(scale < 0)) 
    stop(paste("scale parameter must be positive", "\n", ""))
  cdf <- 1 - exp(-exp((q - location)/scale))
  if (lower.tail == TRUE) 
    cdf <- cdf
  else cdf <- 1 - cdf
  if (log.p == FALSE) 
    cdf <- cdf
  else cdf <- log(cdf)
  cdf
}
## quantile function ------------------------------------------------
qlogWeibull = function(p, location=0, scale=1, lower.tail = TRUE, log.p = FALSE) {
  if (any(scale < 0)) 
    stop(paste("scale parameter must be positive", "\n", ""))
  if (log.p == TRUE) 
    p <- exp(p)
  else p <- p
  if (lower.tail == TRUE) 
    p <- p
  else p <- 1 - p
  if (any(p < 0) | any(p > 1)) 
    stop(paste("p must be between 0 and 1", "\n", ""))
  q <- location + scale * log(-log(1 - p))
  q
}
## random sampler function ------------------------------------------------
rlogWeibull = function(n, location=0, scale=1) {
  if (any(scale < 0)) 
    stop(paste("scale parameter must be positive", "\n", ""))
  if (any(n <= 0)) 
    stop(paste("n must be a positive integer", "\n", ""))
  n <- ceiling(n)
  p <- runif(n)
  r <- qlogWeibull(p, location=location, scale=scale)
  r
}

# link functions ------------------------------------------------
pfunc = function(x, method, df=7){
  if(method == "logistic"){
    return( plogis(x) )
  }else if(method == "probit"){
    return( pnorm(x,mean=0,sd=1) )
  }else if(method == "cauchit"){
    return( pcauchy(x) )
  }else if(method == "tlink"){
    return( pt(x,df) )
  }else if(method == "cloglog"){
    return( plogWeibull(x) )
  }else if(method == "loglog"){
    return( pgumbel(x) )
  }else if(method == "lptn"){
    return( plptn(x) )
  }
}

# density functions ------------------------------------------------
dfunc = function(x, method, df=7){
  if(method == "logistic"){
    return( dlogis(x) )
  }else if(method == "probit"){
    return( dnorm(x,mean=0,sd=1) )
  }else if(method == "cauchit"){
    return( dcauchy(x) )
  }else if(method == "tlink"){
    return( dt(x,df) )
  }else if(method == "cloglog"){
    y = dlogWeibull(x)
    y[which(x==Inf)] = 0
    return( y )
  }else if(method == "loglog"){
    y = dgumbel(x)
    y[which(x==-Inf)] = 0
    return( y )
  }else if(method == "lptn"){
    return( dlptn(x) )
  }
}

# derivative of density functions ------------------------------------------------
pd_dfunc = function(x, method, df=7){
  if(method == "logistic"){
    y = - dlogis(x) * plogis(x)
    return( y )
  }else if(method == "probit"){
    y = -x * dnorm(x,mean=0,sd=1)
    y[which(x==Inf | x==-Inf)] = 0
    return( y )
  }else if(method == "tlink"){
    y = -dt(x,df) * (1+x^2/df)^{-1} * (1+1/df)
    return( y )
  }else if(method == "cauchit"){
    y = -dt(x,df=1) * (1+x^2)^{-1} * 2
    return( y )
  }else if(method == "cloglog"){
    y = exp(x-exp(x)) - exp(2*x-exp(x))
    y[which(x==Inf)] = 0
    return( y )
  }else if(method == "loglog"){
    y = exp(-x-exp(-x)) * {-1 + exp(-x)}
    y[which(x==-Inf)] = 0
  }
}

# random sampler functions ------------------------------------------------
rfunc = function(n, method, df=7){
  if(method == "logistic"){
    return( rlogis(n) )
  }else if(method == "probit"){
    return( rnorm(n,mean=0,sd=1) )
  }else if(method == "cauchit"){
    return( rcauchy(n) )
  }else if(method == "tlink"){
    return( rt(n,df) )
  }else if(method == "cloglog"){
    return( rlogWeibull(n) )
  }else if(method == "loglog"){
    return( rgumbel(n) )
  }
}


