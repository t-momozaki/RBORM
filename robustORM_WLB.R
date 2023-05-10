## Robust Bayesian Ordinal Regression Model via Divergence Approach with weighted likelihood bootstrap

## INPUT
# y (factor): an ordered categorical response vector
# X (matrix): (n,p)-matrix of covariates without intercept term (n: total sample size; p: the number of covariates)
# method (character): name of link function to be used
# divergence (character): name of posterior to be used; "KL": standard posterior with Kullback-Liebler, "DP": density-power (DP) posterior with DP, "g_syn": $\gamma$-synthetic posterior with $\gamma$-divergence, "g_gen": $\gamma$-general posterior with $\gamma$-divergence
# tnp (numeric): the value of tuning parameter for DP and $\gamma$-divergences; defalt is 0.5
# init_parameters (numeric): (p+M-1) vector of initial value of parameters; Input initial values for cutpoints for the first M-1 elements, and initial values for coefficients for the rest.
# Priors (list): list including the prior name, "b_prior", for coefficient parameters and its hyper-parameter, "A" and "B"; defalt is $\beta \sim N(A,B)$ where "A" is rep(0, dim(X)[2] and "B" is diag(rep(100^2, dim(X)[2]))
# MCs (list): list for MCMC iteration including total MCMC length "mc", burn-in "bn", and the number of thinnings "thin"; defalt is "mc" is 2000, "bn" is 0, and "thin" is 1
# parallel (logical): whether to generate posterior samples in parallel; defalt is TRUE

## OUTPUT
# (data.frame): posterior samples of coefficient and cutoff parameters


# source
source("links.R")

# library
list.of.packages <- c("ordinal","mvtnorm","MASS","MCMCpack", "doSNOW", "doRNG", "evd", "parallel", "foreach")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(doSNOW)


robustORM_WLB = function(y, X, method=c("probit", "logistic", "cloglog", "loglog", "cauchit"), 
                         divergence=c("KL", "DP", "g_syn", "g_gen"), tnp, init_parameters, 
                         Priors, MCs, parallel) {
  # Preparation -------------------------------------------------------
  method <- match.arg(method)
  divergence <- match.arg(divergence)
  if (missing(init_parameters)) init_parameters <- NULL
  if(missing(tnp)) tnp = 0.5
  if(missing(Priors)) {
    Priors = list(
      b_prior = "normal",
      A=rep(0, dim(X)[2]), 
      B=diag(rep(100^2, dim(X)[2]))
    )
  }
  if(missing(MCs)) {
    MCs = list(
      mc = 2000,
      bn = 0,
      thin = 1
    )
  }
  if(missing(parallel)) parallel = ifelse(Priors$b_prior=="normal", TRUE, FALSE)
  if(parallel==TRUE & Priors$b_prior!="normal" & Priors$b_prior!="uniform") stop("parallel is only available using the normal and uniform of the coefficients prior!")
  if(tnp==0.0) divergence <- "KL"
  ## Parameters setting -------------------------------------------------------
  n = dim(X)[1] # sample size
  p = dim(X)[2] # the number of coefficient parameters
  M = max( as.numeric(y) ) # the number of categories
  q = p + M-1 # the number of parameters
  ## Setting the initial values -------------------------------------------------------
  if (is.null(init_parameters)) {
    ### Setting the start value for parameters "start <- clm.start(y.levels = rorm.struct$y.levels, threshold = threshold, X = rorm.struct$X, NOM = rorm.struct$NOM, has.intercept = TRUE)" ----------
    #### Setting the start value for cutpoints "st <- start.threshold(y.levels, threshold)" ----------
    init_delta <- as.vector( qlogis((1:(M-1))/((M-1) + 1)) )
    init_delta_tilde <- c(init_delta[1])
    for (m in 2:(M-1)) {
      init_delta_tilde[m] <- sqrt(init_delta[m] - init_delta[m-1])
    }
    #### Setting the start value for coefficients "start.beta(X, has.intercept = TRUE)" ----------
    init_beta <- rep(0, p)
    # init_parameters = c(init_beta, init_delta_tilde)
  }else if(all(init_parameters=="MLE")) {
    d = data.frame(y=as.ordered(y),X)
    fit = ordinal::clm(y~., data=d, link=method)
    init_beta = as.vector(fit$beta)
    init_delta = as.vector(fit$Theta)
    init_delta_tilde <- c(init_delta[1])
    for (m in 2:(M-1)) {
      init_delta_tilde[m] <- sqrt(init_delta[m] - init_delta[m-1])
    }
  }else{
    init_delta = init_parameters[1:(M-1)]
    init_delta_tilde <- c(init_delta[1])
    for (m in 2:(M-1)) {
      init_delta_tilde[m] <- sqrt(init_delta[m] - init_delta[m-1])
    }
    init_beta = init_parameters[-(1:(M-1))]
  }
  ## Setting initial values of optimaization "control <- do.call(clm.control, c(control, list(...)))" ----------------------------------------------
  maxIter = 100L
  gradTol = 1e-06
  relTol = 1e-06
  tol = sqrt(.Machine$double.eps)
  maxLineIter = 15L
  maxModIter = 5L
  ## Setting objective function -------------------------------------------------------
  obj_rborm = switch(divergence,
                     "KL" = function(beta,delta_tilde,y,X,link,Priors,w) {
                       n = dim(X)[1]
                       p = dim(X)[2]
                       M = max( as.numeric(y) )
                       delta = c(delta_tilde[1], delta_tilde[1]+cumsum(delta_tilde[-1]^2), Inf)
                       delta_0_M = c(-Inf, delta[-M])
                       Xb = X%*%beta
                       ll = drop(pfunc(delta[y]-Xb,method=link)-pfunc(delta_0_M[y]-Xb,method=link))
                       -sum(w*log(ll)) +
                         switch (Priors$b_prior,
                           "normal" = - mvtnorm::dmvnorm(beta,mean=Priors$A,sigma=Priors$B,log=T),
                           "uniform" = 0
                         )
                     },
                     "DP" = function(beta,delta_tilde,y,X,link,tnp,Priors,w) {
                       n = dim(X)[1]
                       p = dim(X)[2]
                       M = max( as.numeric(y) )
                       delta = c(delta_tilde[1], delta_tilde[1]+cumsum(delta_tilde[-1]^2), Inf)
                       delta_0_M = c(-Inf, delta[-M])
                       Delta = matrix( rep(delta,n),n,M,T )
                       Delta_0_M = matrix( rep(delta_0_M,n),n,M,T )
                       Xb = X%*%beta
                       XB = matrix( rep(X%*%beta,M),n,M )
                       ll = drop(pfunc(delta[y]-Xb,method=link)-pfunc(delta_0_M[y]-Xb,method=link))
                       LL = pfunc(Delta-XB,method=link)-pfunc(Delta_0_M-XB,method=link)
                       -sum( w * { (1/tnp)*ll^tnp - (1/(1+tnp))*rowSums(LL^(1+tnp)) } ) + 
                         switch (Priors$b_prior,
                                 "normal" = - mvtnorm::dmvnorm(beta,mean=Priors$A,sigma=Priors$B,log=T),
                                 "uniform" = 0
                         )
                     },
                     "g_syn" = function(beta,delta_tilde,y,X,link,tnp,Priors,w) {
                       n = dim(X)[1]
                       p = dim(X)[2]
                       M = max( as.numeric(y) )
                       delta = c(delta_tilde[1], delta_tilde[1]+cumsum(delta_tilde[-1]^2), Inf)
                       delta_0_M = c(-Inf, delta[-M])
                       Delta = matrix( rep(delta,n),n,M,T )
                       Delta_0_M = matrix( rep(delta_0_M,n),n,M,T )
                       Xb = X%*%beta
                       XB = matrix( rep(X%*%beta,M),n,M )
                       ll = drop(pfunc(delta[y]-Xb,method=link)-pfunc(delta_0_M[y]-Xb,method=link))
                       LL = pfunc(Delta-XB,method=link)-pfunc(Delta_0_M-XB,method=link)
                       -(n/tnp)*log( sum( w*{ ll/rowSums(LL^(1+tnp))^(1/(1+tnp)) }^tnp ) ) +
                         switch (Priors$b_prior,
                                 "normal" = - mvtnorm::dmvnorm(beta,mean=Priors$A,sigma=Priors$B,log=T),
                                 "uniform" = 0
                         )
                     },
                     "g_gen" = function(beta,delta_tilde,y,X,link,tnp,Priors,w) {
                       n = dim(X)[1]
                       p = dim(X)[2]
                       M = max( as.numeric(y) )
                       delta = c(delta_tilde[1], delta_tilde[1]+cumsum(delta_tilde[-1]^2), Inf)
                       delta_0_M = c(-Inf, delta[-M])
                       Delta = matrix( rep(delta,n),n,M,T )
                       Delta_0_M = matrix( rep(delta_0_M,n),n,M,T )
                       Xb = X%*%beta
                       XB = matrix( rep(X%*%beta,M),n,M )
                       ll = drop(pfunc(delta[y]-Xb,method=link)-pfunc(delta_0_M[y]-Xb,method=link))
                       LL = pfunc(Delta-XB,method=link)-pfunc(Delta_0_M-XB,method=link)
                       -sum( w*(1/tnp)*{ll/rowSums(LL^(1+tnp))^(1/(1+tnp))}^tnp ) +
                         switch (Priors$b_prior,
                                 "normal" = - mvtnorm::dmvnorm(beta,mean=Priors$A,sigma=Priors$B,log=T),
                                 "uniform" = 0
                         )
                     } )
  
  # Optimization -------------------------------------------------------
  time = proc.time()[3]
  mc = MCs$mc
  bn = MCs$bn
  thin = MCs$thin
  if(parallel) {
    f <- function(iterator){
      pb <- txtProgressBar(min = 1, max = iterator - 1, style = 3)
      count <- 0
      function(...) {
        count <<- count + length(list(...)) - 1
        setTxtProgressBar(pb, count)
        flush.console()
        rbind(...) # this can feed into .combine option of foreach
      }
    }
    cores <- getOption("mc.cores",parallel::detectCores())
    cl <- parallel::makeCluster(cores)
    doSNOW::registerDoSNOW(cl)
    parallel::clusterSetRNGStream(cl, 1)
    post = foreach(r=1:mc,.combine=f(mc),.export=c("pfunc","dfunc"),.packages=c("evd")) %dopar% {
      beta = init_beta
      delta_tilde = init_delta_tilde
      b_diff = 1
      d_diff = 1
      nIter = 0
      w = as.vector( MCMCpack::rdirichlet(n=1, alpha=rep(1,n)) )
      while( (b_diff>tol&d_diff>tol) & nIter < maxIter) {
        old_beta = beta
        old_delta_tilde = delta_tilde
        
        opt_beta = switch(divergence,
                          "KL" = optim(par=beta, fn=obj_rborm, 
                                       delta_tilde=delta_tilde,
                                       y=y, X=X, link=method,
                                       Priors=Priors, w=w),
                          "DP" = optim(par=beta, fn=obj_rborm, 
                                       delta_tilde=delta_tilde,
                                       y=y, X=X, link=method, tnp=tnp,
                                       Priors=Priors, w=w),
                          "g_syn" = optim(par=beta, fn=obj_rborm, 
                                          delta_tilde=delta_tilde,
                                          y=y, X=X, link=method, tnp=tnp,
                                          Priors=Priors, w=w),
                          "g_gen" = optim(par=beta, fn=obj_rborm, 
                                          delta_tilde=delta_tilde,
                                          y=y, X=X, link=method, tnp=tnp,
                                          Priors=Priors, w=w) )
        beta = opt_beta$par
        
        opt_delta_tilde = switch(divergence,
                                 "KL" = optim(par=delta_tilde, fn=obj_rborm,
                                              beta=beta,
                                              y=y, X=X, link=method,
                                              Priors=Priors, w=w),
                                 "DP" = optim(par=delta_tilde, fn=obj_rborm,
                                              beta=beta,
                                              y=y, X=X, link=method, tnp=tnp,
                                              Priors=Priors, w=w),
                                 "g_syn" = optim(par=delta_tilde, fn=obj_rborm,
                                                 beta=beta,
                                                 y=y, X=X, link=method, tnp=tnp,
                                                 Priors=Priors, w=w),
                                 "g_gen" = optim(par=delta_tilde, fn=obj_rborm,
                                                 beta=beta,
                                                 y=y, X=X, link=method, tnp=tnp,
                                                 Priors=Priors, w=w) )
        delta_tilde = opt_delta_tilde$par
        
        b_diff = sqrt(crossprod(beta-old_beta)/p)
        d_diff = sqrt(crossprod(delta_tilde-old_delta_tilde)/{M-1})
        nIter = nIter+1
        
        if(nIter%%20==0){
          cat(nIter,"th iteration", "\n",
              proc.time()[3]-time, "seconds gone.", "\n",
              "convergence criterion values: \n", 
              " for beta: ", b_diff, "\n", 
              " for delta: ", d_diff, "\n", "\n")
        } 
      }
      delta = c(delta_tilde[1], delta_tilde[1]+cumsum(delta_tilde[-1]^2))
      c(delta, beta)
    }
    stopCluster(cl)
    cat("The total time required is", proc.time()[3]-time, ".\n")
  }else {
    post = matrix(NA, mc, q)
    for (r in 1:mc) {
      beta = init_beta
      delta_tilde = init_delta_tilde
      b_diff = 1
      d_diff = 1
      nIter = 0
      w = as.vector( MCMCpack::rdirichlet(n=1, alpha=rep(1,n)) )
      while( (b_diff>tol&d_diff>tol) & nIter < maxIter) {
        old_beta = beta
        old_delta_tilde = delta_tilde
        
        opt_beta = switch(divergence,
                          "KL" = optim(par=beta, fn=obj_rborm, 
                                       delta_tilde=delta_tilde,
                                       y=y, X=X, link=method,
                                       Priors=Priors, w=w),
                          "DP" = optim(par=beta, fn=obj_rborm, 
                                       delta_tilde=delta_tilde,
                                       y=y, X=X, link=method, tnp=tnp,
                                       Priors=Priors, w=w),
                          "g_syn" = optim(par=beta, fn=obj_rborm, 
                                          delta_tilde=delta_tilde,
                                          y=y, X=X, link=method, tnp=tnp,
                                          Priors=Priors, w=w),
                          "g_gen" = optim(par=beta, fn=obj_rborm, 
                                          delta_tilde=delta_tilde,
                                          y=y, X=X, link=method, tnp=tnp,
                                          Priors=Priors, w=w) )
        beta = opt_beta$par
        
        opt_delta_tilde = switch(divergence,
                                 "KL" = optim(par=delta_tilde, fn=obj_rborm,
                                              beta=beta,
                                              y=y, X=X, link=method,
                                              Priors=Priors, w=w),
                                 "DP" = optim(par=delta_tilde, fn=obj_rborm,
                                              beta=beta,
                                              y=y, X=X, link=method, tnp=tnp,
                                              Priors=Priors, w=w),
                                 "g_syn" = optim(par=delta_tilde, fn=obj_rborm,
                                                 beta=beta,
                                                 y=y, X=X, link=method, tnp=tnp,
                                                 Priors=Priors, w=w),
                                 "g_gen" = optim(par=delta_tilde, fn=obj_rborm,
                                                 beta=beta,
                                                 y=y, X=X, link=method, tnp=tnp,
                                                 Priors=Priors, w=w) )
        delta_tilde = opt_delta_tilde$par
        
        b_diff = sqrt(crossprod(beta-old_beta)/p)
        d_diff = sqrt(crossprod(delta_tilde-old_delta_tilde)/{M-1})
        nIter = nIter+1
        
        if(nIter%%20==0 & r==1){
          cat(nIter,"th iteration", "\n",
              proc.time()[3]-time, "seconds gone.", "\n",
              "convergence criterion values: \n", 
              " for beta: ", b_diff, "\n", 
              " for delta: ", d_diff, "\n", "\n")
        } 
      }
      delta = c(delta_tilde[1], delta_tilde[1]+cumsum(delta_tilde[-1]^2))
      post[r,] = c(delta, beta)
      ## progress bar -----
      if(any(r*100/mc==seq(10,100,by=10))){
        cat(r*100/mc, "% of MCMC done", "\n")
        if(r*100/mc==10) cat("The time taken for the MCMC for 10% of is", proc.time()[3]-time, "\n")
        cat("The total time required up to this point is", proc.time()[3]-time, "\n")
      }
    }
  }
  post = post[seq(bn+thin,mc,thin),]
  colnames(post) = c(paste0("delta",1:(M-1)), paste0("beta",1:p))
  
  # results -------------------------------------------------------
  data.frame(post)
}

