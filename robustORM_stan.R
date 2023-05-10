## Robust Bayesian Ordinal Regression Model via Divergence Approach with Stan

## INPUT
# y (factor): an ordered categorical response vector
# X (matrix): (n,p)-matrix of covariates without intercept term (n: total sample size; p: the number of covariates)
# method (character): name of link function to be used
# divergence (character): name of posterior to be used; "KL": standard posterior with Kullback-Liebler, "DP": density-power (DP) posterior with DP, "g_syn": $\gamma$-synthetic posterior with $\gamma$-divergence, "g_gen": $\gamma$-general posterior with $\gamma$-divergence
# tnp (numeric): the value of tuning parameter for DP and $\gamma$-divergences; defalt is 0.5
# init_parameters (numeric): (p+M-1) vector of initial value of parameters; Input initial values for cutpoints for the first M-1 elements, and initial values for coefficients for the rest.
# Priors (list): list including the prior name, "b_prior", for coefficient parameters and its hyper-parameter, "A" and "B"; defalt is $\beta \sim N(A,B)$ where "A" is rep(0, dim(X)[2] and "B" is diag(rep(100^2, dim(X)[2]))
# MCs (list): list for MCMC iteration including total MCMC length "mc", burn-in "bn", and the number of thinnings "thin"; defalt is "mc" is 2500, "bn" is 500, and "thin" is 1

## OUTPUT
# object of stan 


# library
list.of.packages <- c("rstan")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(rstan)


robustORM_stan = function(y, X, method=c("probit","logistic","cloglog","loglog"), 
                          divergence=c("KL","DP","g_syn","g_gen"), tnp, init_parameters, Priors, MCs) {
  # Setting parameters -------------------------------------------------------
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
      mc = 2500,
      bn = 500,
      thin = 1,
      chains = 1,
      seed = 1234
    )
  }
  if(tnp==0.0) divergence <- "KL"
  
  # Preparation for data of stan ---------------------------------------------------
  data = switch (divergence,
    "KL" = list(M = max(as.numeric(y)),
                n = dim(X)[1],
                p = dim(X)[2],
                y = as.numeric(y),
                x = as.matrix(X) ),
    "DP" = list(M = max(as.numeric(y)),
                n = dim(X)[1],
                p = dim(X)[2],
                y = as.numeric(y),
                x = as.matrix(X),
                alpha = tnp ),
    "g_syn" = list(M = max(as.numeric(y)),
                   n = dim(X)[1],
                   p = dim(X)[2],
                   y = as.numeric(y),
                   x = as.matrix(X),
                   gamma = tnp ),
    "g_gen" = list(M = max(as.numeric(y)),
                   n = dim(X)[1],
                   p = dim(X)[2],
                   y = as.numeric(y),
                   x = as.matrix(X),
                   gamma = tnp ),
    stop("Please select the appropriate 'divergence'!", "\n",
         "if you'd like to the 'standard posterior', then you should take div='KL',", "\n",
         "if you'd like to the 'DP-posterior', then you should take div='DP',", "\n", 
         "if you'd like to the 'gamma-synthetic posterior', then you should take div='g_syn',", "\n",
         "and, if you'd like to the 'gamma-general posterior', then you should take div='g_gen'.", "\n")
  )
  
  # Setting the initial values -------------------------------------------------------
  if (is.null(init_parameters)) {
    ## Setting the start value for parameters "start <- clm.start(y.levels = rorm.struct$y.levels, threshold = threshold, X = rorm.struct$X, NOM = rorm.struct$NOM, has.intercept = TRUE)" ----------
    ### Setting the start value for cutpoints "st <- start.threshold(y.levels, threshold)" ----------
    init_delta <- as.vector( qlogis((1:(data$M-1))/((data$M-1) + 1)) )
    #### Setting the start value for coefficients "start.beta(X, has.intercept = TRUE)" ----------
    init_beta <- rep(0, data$p)
    init_parameters = c(init_beta, init_delta)
  }else if(all(init_parameters=="MLE")) {
    d = data.frame(y=as.ordered(data$y),data$x)
    fit = ordinal::clm(y~., data=d, link=method)
    init_beta = as.vector(fit$beta)
    init_delta = as.vector(fit$Theta)
    init_parameters = c(init_beta, init_delta)
  }
  
  # Preparation for stan file ---------------------------------------------------
  file = switch (divergence,
    "KL" = paste0("orm_",method),
    "DP" = paste0("rorm_dp_",method),
    "g_syn" = paste0("rorm_g_syn_",method),
    "g_gen" = paste0("rorm_g_gen_",method)
  )
  stan_file = paste0("robustORM_stan/", file, ".stan")
  
  # Fitting model ---------------------------------------------------
  model = rstan::stan_model(file=stan_file)
  fit = rstan::sampling(
    model,
    data = data,
    seed = MCs$seed,
    chains = MCs$chains,
    iter = MCs$mc,
    warmup = MCs$bn,
    thin = MCs$thin,
    init = init_parameters,
    cores = getOption("mc.cores", parallel::detectCores())
  )
  return(fit)
}
