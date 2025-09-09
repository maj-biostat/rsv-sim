data{
  int N;
  array[N] int y;
  array[N] int n;
  
  // intervention
  array[N] int trt;
  
  // region (Alice, Darwin)
  array[N] int reg;
  // locality (urban, remote)
  array[N] int loc;
  
  // severity?
  
  // priors
  vector[2] pri_b_0;
  vector[2] pri_b_trt;
  vector[2] pri_b_reg;
  vector[2] pri_b_loc;
  
  int prior_only;
}
parameters{
  real b_0;
  real b_trt_raw;
  real b_reg_raw;
  real b_loc_raw;
  
}
transformed parameters{
  
  vector[2] b_trt;
  vector[2] b_reg;
  vector[2] b_loc;
  vector[N] eta;
  
  b_trt[1] = 0.0;
  b_reg[1] = 0.0;
  b_loc[1] = 0.0;
  
  b_trt[2] = b_trt_raw;
  b_reg[2] = b_reg_raw;
  b_loc[2] = b_loc_raw;
  
  eta = b_0 + b_trt[trt] + b_reg[reg] + b_loc[loc];
  
}
model{
  
  target += logistic_lpdf(b_0 | pri_b_0[1], pri_b_0[2]);
  target += normal_lpdf(b_trt_raw | pri_b_trt[1], pri_b_trt[2]);
  target += normal_lpdf(b_reg_raw | pri_b_reg[1], pri_b_reg[2]);
  target += normal_lpdf(b_loc_raw | pri_b_loc[1], pri_b_loc[2]);
  
  if(!prior_only){
    target += binomial_logit_lpmf(y | n, eta);  
  }
  
}
generated quantities{
  
  vector[N] w = dirichlet_rng(to_vector(n));
  
  vector[N] mu_1 = inv_logit(b_0 + b_trt[1] + b_reg[reg] + b_loc[loc]);
  vector[N] mu_2 = inv_logit(b_0 + b_trt[2] + b_reg[reg] + b_loc[loc]);
  
  real p_1 = w' * mu_1   ;
  real p_2 = w' * mu_2   ;
  
  // start with the marginal effect of treatment
  real rd_2_1 = p_2 - p_1;
  
}
