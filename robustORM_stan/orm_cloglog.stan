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
  for(i in 1:n) {
    real eta;
    eta = x[i] * beta;
    pi[1] = 1 - exp(-exp(delta[1] - eta));
    for (m in 2:(M-1)) {
      pi[m] = - exp(-exp(delta[m] - eta)) + exp(-exp(delta[m-1] - eta));
    }
    pi[M] = exp(-exp(delta[M-1] - eta));
    y[i] ~ categorical(pi);
  }
}

