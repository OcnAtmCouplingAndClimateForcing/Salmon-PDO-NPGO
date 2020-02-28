//Dynamic Linearized Ricker Model - Autoregressive
//  This model is fit to multiple populations and estimates separate time-varying coefficients each region. 

//Notes:
//  matrix[3, 3] m[6, 7] - m to be a two-dimensional array of size 6 x 7,
//    containing values that are 3 x 3 matrices. 

data {
  
  int<lower=0> S; //number of stocks
  int N[S];  //number of SR observations for each stock
  int<lower=0> maxN;  //maximum number of observations across stocks
  int years[S, maxN]; //Pointer vector for brood year references.
  matrix[S,maxN] ln_rps;
  matrix[S,maxN] spawn;
  int<lower=0> K; // number of covariates
  // vector[K] covars[S,maxN]; //[S,N,K]
  // matrix[S,maxN] covars[K]; //[K,S,N]
  //real[S,maxN,K] covars;
  //int log[N];
  
  matrix[S,maxN] PDO;
  matrix[S,maxN] NPGO;
  
  int<lower=0> R; //number of regions
  int region[S]; //region pointer
  
}
parameters {
  //Coefficients
  //Initial Coefficient Values
  vector[R] coef0_PDO;
  vector[R] coef0_NPGO;
  
  vector[maxN-1] pro_dev_PDO[R]; // elements accessed [R,N-1] - A vector of length R, each element of which is a vector of lenth maxN-1
  vector[maxN-1] pro_dev_NPGO[R];
  
  //Ricker Params
  real<lower=0> alpha[S];
  real<lower=0> beta[S];
  //Variances
  real<lower=0> sigma_pe_PDO[R];
  real<lower=0> sigma_pe_NPGO[R];
  real<lower=0> sigma_oe[S];
  //ARMA(1) Components
  vector[R] phi_PDO;  //AR(1) Component
  vector[R] phi_NPGO;
  vector[R] theta_PDO;  //MA(1) Component
  vector[R] theta_NPGO;
}
transformed parameters {
  vector[maxN] pred[S]; // [S,maxN] - vector of length S, each element of which is a vector of length maxN
  vector[maxN] coef_PDO[R]; // elements accessed [R,maxN] - vector of length R, each element is a vector of length maxN
  vector[maxN] coef_NPGO[R];
  
  //Coefficient Time Series
  for(r in 1:R) {  //Loop through regions
    coef_PDO[r,1] = coef0_PDO[r];  //Initial time step
    coef_NPGO[r,1] = coef0_NPGO[r];
    for(n in 2:maxN) {
      //phi is for AR(1) process - Default to 1 if not used
      coef_PDO[r,n] = phi_PDO[r]*coef_PDO[r,n-1] + pro_dev_PDO[r,n-1];
      coef_NPGO[r,n] = phi_NPGO[r]*coef_NPGO[r,n-1] + pro_dev_NPGO[r,n-1];
    }//next n
  }//next r
  
  //Generate Predicted values
  for(s in 1:S) {
    // for(n in 1:N[s]) {  //Only loop through observed values
    // for(n in 1:N[s]) {  //Only loop through observed values
    //   //NOTE: We use the "region" pointer to identify which coefficient time-series to use.
    //   pred[s,n] = alpha[s] - beta[s]*spawn[s,n] + PDO[s,n] * coef_PDO[region[s],n] + NPGO[s,n] * coef_NPGO[region[s],n];
    // }//next n
    // //Fill in Missing Values
    // // pred[s,(N[s]+1):maxN] = 0;
    // for(n in (N[s]+1):maxN) {
    //   pred[s,n] = 0;
    // }//next n
    
    //Alternatively, we may be able to do this with an if statement.... But it makes the JAGSer inside of me cry :(
    for(n in 1:maxN) {
      if(n<=N[s]) {
        // pred[s,n] = alpha[s] - beta[s]*spawn[s,n] + PDO[s,n] * coef_PDO[region[s],n] + NPGO[s,n] * coef_NPGO[region[s],n];
        pred[s,n] = alpha[s] - beta[s]*spawn[s,n] + PDO[s,n] * coef_PDO[region[s],years[s,n]] + NPGO[s,n] * coef_NPGO[region[s],years[s,n]]; //Updated with years pointer vector
      }else {
        pred[s,n] = 0;
      }
    }// next n
  }//next s
}
model {
  //Priors
  for(s in 1:S) {
    alpha[s] ~ normal(0,10);
    beta[s] ~ normal(0,0.001);
    sigma_oe[s] ~ normal(0,5);//cauchy(0,5);
  }//next s
  
  for(r in 1:R) { //Regions
    coef0_PDO[r] ~ normal(0,1);
    coef0_NPGO[r] ~ normal(0,1);
    sigma_pe_PDO[r] ~ normal(0,5);//cauchy(0,5);
    sigma_pe_NPGO[r] ~ normal(0,5);//cauchy(0,5);
    
    //Priors for ARMA(1) Components
    phi_PDO[r] ~ uniform(-0.99,0.99);#normal(0,2);
    phi_NPGO[r] ~ uniform(-0.99,0.99);#normal(0,2);
    theta_PDO[r] ~ uniform(-0.99,0.99);#normal(0,2);
    theta_NPGO[r] ~ uniform(-0.99,0.99);#normal(0,2);
    
    //Deviations - theta is for MA(1) process
    //Initial
    pro_dev_PDO[r,1] ~ normal(0,sigma_pe_PDO[r]);
    pro_dev_NPGO[r,1] ~ normal(0,sigma_pe_NPGO[r]);
    for(n in 2:(maxN-1)) {
      pro_dev_PDO[r,n] ~ normal(theta_PDO[r]*pro_dev_PDO[r,n-1], sigma_pe_PDO[r]);
      pro_dev_NPGO[r,n] ~ normal(theta_NPGO[r]*pro_dev_NPGO[r,n-1], sigma_pe_NPGO[r]);
    }//next n
  }
  //Likelihood
  for(s in 1:S) {
    // ln_rps[s,1:N[s]] ~ normal(pred[1:N[s]], sigma_oe[s]); //DOESN'T WORK!!!!
    for(n in 1:N[s]) {
      ln_rps[s,n] ~ normal(pred[s,n], sigma_oe[s]);
    }//next n
  }//next s
}
generated quantities {
  vector[maxN] log_lik[S];
  for(s in 1:S) {
    // for(n in 1:N[s]){
    for(n in 1:maxN) {
      if(n<=N[s]) {
        log_lik[s,n] = normal_lpdf(ln_rps[s,n] | pred[s,n], sigma_oe[s]);
      }else {
        log_lik[s,n] = 0;
      }
    }//next n
  }//next s
}

