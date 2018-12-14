//Dynamic Linearized Ricker Model

//Notes:
//  matrix[3, 3] m[6, 7] - m to be a two-dimensional array of size 6 x 7,
//    containing values that are 3 x 3 matrices. 

data {
  int<lower=0> N;
  vector[N] ln_rps;
  vector[N] spawn;
  int<lower=0> K; // number of covariates
  matrix[N, K] covars; // matrix of covariates
  //int log[N];
}
parameters {
  //Coefficients
  vector[K] coef0;  //Initial Coefficient Values
  vector[K] pro_dev[N-1]; // elements accessed [N-1,K] - A vector of length N-1, each element of which is a vector of lenth K
  
  //Ricker Params
  real<lower=0> alpha;
  real<lower=0> beta;
  //Variances
  real<lower=0> sigma_pe[K];
  real<lower=0> sigma_oe;
}
transformed parameters {
  vector[N] pred;
  vector[K] coef[N]; // elements accessed [N,K] - vector of length N, each element is a vector of length K
  for(k in 1:K) {  //Loop through coefficients
    coef[1,k] = coef0[k];
    for(n in 2:N) {
      coef[n,k] = coef[n-1,k] + pro_dev[n-1,k];
    }
  }
  for(n in 1:N) {
    pred[n] = alpha - beta*spawn[n] + covars[n] * coef[n];
  }
}
model {
  //Priors
  alpha ~ normal(0,10);
  beta ~ normal(0,0.001);
  sigma_oe ~ normal(0,5);//cauchy(0,5);
  for(k in 1:K) {
    coef0[k] ~ normal(0,1);
    sigma_pe[k] ~ normal(0,5);//cauchy(0,5);
    for(n in 1:(N-1)) {
      pro_dev[n,k] ~ normal(0, sigma_pe[k]);
    }//next n
  }
  //Likelihood
  ln_rps ~ normal(pred, sigma_oe);
}
generated quantities {
  vector[N] log_lik;
  for(n in 1:N){ 
    log_lik[n] = normal_lpdf(ln_rps[n] | pred[n], sigma_oe);
  }
}

