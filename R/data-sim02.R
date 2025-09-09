library(data.table)







get_sim02_trial_data <- function(
    l_spec
){
  # cohort
  if(is.null(l_spec$ic)){ ic <- 1 } else { ic <- l_spec$ic }
  if(is.null(l_spec$is)){ is <- 1 } else { is <- l_spec$is }
  if(is.null(l_spec$ie)){ ie <- 1 } else { ie <- l_spec$ie }
  if(is.null(l_spec$t0)){ t0 <- 1 } else { t0 <- l_spec$t0 }
  if(is.null(l_spec$tfu)){ tfu <- 1 } else { tfu <- l_spec$tfu }
  
  d <- data.table(
    ic = ic,
    id = is:ie,
    t0 = t0,
    reg = rbinom(l_spec$N[ic], 1, prob = l_spec$p_reg_alloc) + 1
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
  
  d[, ia := NA_integer_]
  d[, t_anlys := NA_real_]
  
  # Given the simplicity of the model we can specify the linear predictor 
  # directly in terms of risk increments.
  
  d[, tfu := tfu]
  d[, p := l_spec$b_0 + l_spec$b_reg[reg] + l_spec$b_loc[loc] + l_spec$b_trt[trt]]
  d[, y := rbinom(.N, 1, p)]
  
  d
}



get_sim02_stan_data <- function(d_mod){
  
  # convert from binary representation to binomial (successes/trials)
  dd <- copy(d_mod)
  
  dd <- dd[, .(y = sum(y), n = .N, eta = unique(p)), 
                 keyby = .(reg, loc, trt)]
  
  dd[, p_obs := y / n]
  
  ld <- list(
    # full dataset
    N = nrow(dd), 
    y = dd[, y], 
    n = dd[, n], 
    trt = dd[, trt], 
    reg = dd[, reg], 
    loc = dd[, loc],
    prior_only = 0
  )
  
  list(
    dd = dd,
    ld = ld
  )
  
}