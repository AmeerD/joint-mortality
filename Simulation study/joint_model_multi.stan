data {
  int n;       //number of observations
  int n_state; //number of states
  int n_code;  //number of county codes
  int n_race;   //number of races to model
  int n_year;  //number of years to be modelled
  int n_ages;  //number of age groups in the mortality curve 
  
  int state[n_code]; //state corresponding to counties
  int code[n];       //observed county code
  int race[n];        //observed race
  int year[n];       //observed year 
  int ages[n];       //observed age category
  
  vector[n] deaths;    //observed deaths
  vector[n] pop;       //observed population
  
  int n_pcs; //number of principal components
  matrix[n_ages, n_pcs] pcs; //principal components
}

parameters {
  matrix[n_pcs, n_race] eps[n_code, n_year];      // one pc x sex matrix per county-year
  matrix[n_ages, n_race] u[n_code, n_year];       // one age x sex error matrix per county-year
  
  matrix[n_pcs, n_race] mu_beta[n_state, n_year]; //matrix[n_pcs, n_race] mu_beta[n_state, n_year]; // one pc vector of mu_beta intercepts for year 1
  matrix<lower=0>[n_pcs, n_year] sigma_beta[n_state]; // one pc x year matrix per state (sigmas shared by sexes)
  vector<lower=0>[n_ages] sigma_x;
  cholesky_factor_corr[n_race] beta_chol[n_pcs, n_year];
  cholesky_factor_corr[n_race] chol[n_ages, n_year];
  vector<lower=0>[n_pcs] sigma_mu_beta;

//Comments
// - Added state grouping for future use
// - Added a sex index into the matrices, eventually we will want that structure in place to change the generating
//   procedure from univariate to bivariate (maybe the beta generation should be along all pcs and have different matrices for each sex?)
}

transformed parameters {
  matrix[n_pcs, n_race] beta[n_code, n_year];    // one pc x sex coefficient matrix per county-year
  matrix[n_ages, n_race] logmx[n_code, n_year];  // one age x sex log-mort matrix per county-year  
  
  for (co in 1:n_code) {
    for (yr in 1:n_year) {
      for (p in 1:n_pcs) {
        beta[co, yr, p] = mu_beta[state[co], yr, p] + sigma_beta[state[co], p, yr] * eps[co, yr, p] * beta_chol[p, yr]';
      }
      // beta[co, yr] = mu_beta[state[co], yr] + diag_pre_multiply(sigma_beta[state[co]][:, yr],eps[co, yr]);
      for (a in 1:n_ages) {
        logmx[co, yr, a, :] = pcs[a, :] * beta[co, yr] + sigma_x[a] * u[co, yr, a, :]  * chol[a, yr]';
      }
    }
  }
  
//Comments
// - This is going to seem slower than the model we were given but structurally it makes sense to have the inner matrices be
//   sex x pc instead of pc x year so that joint data generation is easier. Also, if we do decide to introduce correlated
//   coefficients (which I kind of like the idea of), I don't think we can sum over the entire matrix in one go regardless.
//   (slower because of the nested year for loop)
}

model {
// SAMPLE LIKELIHOOD (POISSON: sum over i of -Ni*mi + yi*log mi )  
// this is equiv. to yi ~ poisson(Ni*mi) i=1...N
// but more efficient

  {
    vector[n] logmx_ord;
    for (i in 1:n){
      logmx_ord[i] = logmx[code[i], year[i]][ages[i], race[i]];
    }
    target += -dot_product(pop, exp(logmx_ord)) + dot_product(deaths, logmx_ord);
  }

// PRIORS
  
  for (co in 1:n_code) {
    for (yr in 1:n_year) {
      to_vector(eps[co, yr]) ~ std_normal();
      to_vector(u[co, yr]) ~ std_normal(); 
    }
  }
  
  sigma_x ~ normal(0, 0.25);
  sigma_mu_beta ~ lognormal(-1.5, 0.5);
  
  for (a in 1:n_ages){
    for (yr in 1:n_year) {
      chol[a, yr] ~ lkj_corr_cholesky(1);
    }
  }
  
  
  for (s in 1:n_state) {
    for (p in 1:n_pcs) {
      // for (yr in 2:n_year) {
      //   mu_beta[s, yr, p] ~ multi_normal_cholesky(mu_beta[s, yr-1, p], sigma_mu_beta[p] * mu_chol[s,p]);
      // }
      for (x in 1:n_race) {
        to_vector(mu_beta[s,3:n_year,p,x]) ~ normal(2 * to_vector(mu_beta[s,2:(n_year-1),p,x]) - to_vector(mu_beta[s,1:(n_year-2),p,x]), sigma_mu_beta[p]);
      }
      for (yr in 1:n_year) {
        beta_chol[p, yr] ~ lkj_corr_cholesky(1);
      }
      
    }
    
    to_vector(sigma_beta[s]) ~ std_normal();
  }
}

generated quantities {
  corr_matrix[n_race] cor[n_ages, n_year];
  corr_matrix[n_race] beta_cor[n_pcs, n_year];
  
  for (a in 1:n_ages) {
    for (yr in 1:n_year) {
      cor[a, yr] = chol[a, yr] * chol[a, yr]';
    }
  }
  
  for (p in 1:n_pcs) {
    for (yr in 1:n_year) {
      beta_cor[p, yr] = beta_chol[p, yr] * beta_chol[p, yr]';
    }
  }
  
  // vector[n] yrepl;
  // 
  // for (i in 1:n) {
  //   yrepl[i] = poisson_rng(exp(logmx[code[i], year[i]][ages[i], race[i]]) * pop[i]);
  // }
}
