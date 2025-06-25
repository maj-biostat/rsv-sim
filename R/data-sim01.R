library(data.table)




get_sim01_enrol_time <- function(N = 2500, lambda = 1.52,
                               rho = function(t) pmin(t/360, 1)){
  
  c(0, poisson::nhpp.event.times(lambda, N - 1, rho))
}



get_sim01_trial_data <- function(
    l_spec
){
  
  if(is.null(l_spec$ia)){
    ia <- 1
  } else {
    ia <- l_spec$ia
  }
  if(is.null(l_spec$is)){
    is <- 1
  } else {
    is <- l_spec$is
  }
  if(is.null(l_spec$ie)){
    ie <- 1
  } else {
    ie <- l_spec$ie
  }
  if(is.null(l_spec$t0)){
    t0 <- 1
  } else {
    t0 <- l_spec$t0
  }
  
  d <- data.table(
    ia = ia,
    id = is:ie,
    t0 = t0,
    reg = rbinom(l_spec$N[ia], 1, prob = l_spec$p_reg_alloc) + 1
  )
  d[reg == 1, loc := rbinom(.N, 1, prob = l_spec$p_loc_alloc_a) + 1]
  d[reg == 2, loc := rbinom(.N, 1, prob = l_spec$p_loc_alloc_d) + 1]
  
  # alice, urban
  d[reg == 1 & loc == 1, trt := sample(rep(1:2, length = .N), size = .N, replace = F)]
  # alice, remote
  d[reg == 1 & loc == 2, trt := sample(rep(1:2, length = .N), size = .N, replace = F)]
  # darwin, urban
  d[reg == 2 & loc == 1, trt := sample(rep(1:2, length = .N), size = .N, replace = F)]
  # darwin, remote
  d[reg == 2 & loc == 2, trt := sample(rep(1:2, length = .N), size = .N, replace = F)]
  
  
  # Given the simplicity of the model we can specify the linear predictor 
  # directly in terms of risk increments.
  d[, p := l_spec$bmu + l_spec$breg[reg] + l_spec$bloc[loc] + l_spec$btrt[trt]]
  d[, y := rbinom(.N, 1, p)]
  
  d
}



get_sim01_stan_data <- function(d_all){
  
  # convert from binary representation to binomial (successes/trials)
  d_mod <- d_all[, .(y = sum(y), n = .N, eta = unique(p)), 
                 keyby = .(reg, loc, trt)]
  
  d_mod[, p_obs := y / n]
  
  ld <- list(
    # full dataset
    N = nrow(d_mod), 
    y = d_mod[, y], 
    n = d_mod[, n], 
    trt = d_mod[, trt], 
    reg = d_mod[, reg], 
    loc = d_mod[, loc],
    prior_only = 0
  )
  
  list(
    d_mod = d_mod,
    ld = ld
  )
  
}