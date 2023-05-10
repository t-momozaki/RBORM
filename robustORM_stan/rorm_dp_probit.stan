//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
  int<lower=2> M;
  int<lower=0> n;
  int<lower=1> p;
  int<lower=1,upper=M> y[n];
  row_vector[p] x[n];
  real<lower=0> alpha;
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[p] beta;
  ordered[M-1] delta;
}

// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  vector[M] pi;
  vector[M] pi_alpha;
  for(i in 1:n) {
    real eta;
    eta = x[i] * beta;
    pi[1] = Phi(delta[1] - eta);
    pi_alpha[1] = pow(pi[1], 1+alpha);
    for (m in 2:(M-1)) {
      pi[m] = Phi(delta[m] - eta) - Phi(delta[m-1] - eta);
      pi_alpha[m] = pow(pi[m], 1+alpha);
    }
    pi[M] = 1 - Phi(delta[M-1] - eta);
    pi_alpha[M] = pow(pi[M], 1+alpha);
    
    if(y[i]==1){
      target += (1/alpha)*pi[1]^alpha - (1/(1+alpha))*sum(pi_alpha);
      // increment_log_prob(
      //   (1/alpha)*pi[1]^alpha - (1/(1+alpha))*sum(pi_alpha)
      // );
    }else if(y[i]==M){
      target += (1/alpha)*pi[M]^alpha - (1/(1+alpha))*sum(pi_alpha);
      // increment_log_prob(
      //   (1/alpha)*pi[M]^alpha - (1/(1+alpha))*sum(pi_alpha)
      // );
    }else{
      target += (1/alpha)*pi[y[i]]^alpha - (1/(1+alpha))*sum(pi_alpha);
      // increment_log_prob(
      //   (1/alpha)*pi[y[i]]^alpha - (1/(1+alpha))*sum(pi_alpha)
      // );
    }
  }
}
