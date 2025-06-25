# Experiment with independent set of models with reduced linear predictor.

source("./R/init.R")
source("./R/data-sim01.R")

# Command line arguments list scenario (true dose response),
# the number of simulations to run, the number of cores
# to use and the simulator to use.
args = commandArgs(trailingOnly=TRUE)

# Load cfg based on cmd line args.
if (length(args)<1) {
  log_info("Setting default run method (does nothing)")
  args[1] = "run_none_sim01"
  args[2] = "./sim01/cfg-sim01-v01.yml"
} else {
  log_info("Run method ", args[1])
  log_info("Scenario config ", args[2])
}

# Log setup
f_log_sim <- file.path("./logs", "log-sim.txt")
log_appender(appender_file(f_log_sim))
log_info("*** START UP ***")
# log_threshold(TRACE)

f_cfgsc <- file.path("./etc", args[2])
g_cfgsc <- config::get(file = f_cfgsc)
stopifnot("Config is null" = !is.null(g_cfgsc))

ix <- 1
m1 <- cmdstanr::cmdstan_model("stan/sim01-v01.stan")

output_dir_mcmc <- paste0(getwd(), "/tmp")



# Main trial loop.
run_trial <- function(
    ix,
    l_spec,
    return_posterior = F
){
  
  log_info("Entered  run_trial for trial ", ix)
  
  # Get enrolment times for arbitrary large set of patients
  # Simpler to produce enrol times all in one hit rather than try to do them 
  # incrementally with a non-hom pp.
  
  # events per day
  lambda = 1.0
  # ramp up over 12 months 
  rho = function(t) pmin(t/360, 1)
  
  loc_t0 <- get_sim01_enrol_time(sum(l_spec$N), lambda, rho)
  
  # loop controls
  stop_enrol <- FALSE
  l_spec$ia <- 1 # interim number
  N_analys <- length(l_spec$N)
  
  # posterior summaries
  d_post_smry_1 <- CJ(
    ia = 1:N_analys,
    par = c("p0", "p1", "rd", "rr")
  )
  d_post_smry_1[, mu := NA_real_]
  d_post_smry_1[, med := NA_real_]
  d_post_smry_1[, se := NA_real_]
  d_post_smry_1[, q_025 := NA_real_]
  d_post_smry_1[, q_975 := NA_real_]
  
  
  # decisions 
  # superior, ni, 
  # futile (for superiority - idiotic)
  # futile (for ni - idiotic)
  g_dec_type <- c("sup", 
                  "fut"
  )
  
  # superiority probs
  pr_dec <- array(
    NA,
    # num analysis x num domains x num decision types
    dim = c(N_analys, length(g_dec_type)),
    dimnames = list(
      1:N_analys,  g_dec_type)
  )
  
  decision <- array(
    NA,
    # num analysis x num domains x num decision types
    dim = c(N_analys, length(g_dec_type)),
    dimnames = list(
      1:N_analys,  g_dec_type)
  )
  
  # store all simulated trial pt data
  d_all <- data.table()
  
  if(return_posterior){
    d_post_all <- data.table()
  }
  
  
  ## LOOP -------
  while(!stop_enrol){
    
    log_info("Trial ", ix, " analysis ", l_spec$ia)
    
    # next chunk of data on pts.
    if(l_spec$ia == 1){
      # starting pt index in data
      l_spec$is <- 1
      l_spec$ie <- l_spec$is + l_spec$N[l_spec$ia] - 1
    } else {
      l_spec$is <- l_spec$ie + 1
      l_spec$ie <- l_spec$is + l_spec$N[l_spec$ia] - 1
    }
    
    # id and time
    l_spec$t0 = loc_t0[l_spec$is:l_spec$ie]
    
    # Our analyses only occur on those that have reached 12 months post 
    # randomisation. As such, we are assuming that the analysis takes place
    # 12 months following the last person to be enrolled in the current 
    # analysis set.
    
    d <- get_sim01_trial_data(l_spec)
    
    log_info("Trial ", ix, " new data generated ", l_spec$ia)
    
    # combine the existing and new data
    d_all <- rbind(d_all, d)
    
    # create stan data format based on the relevant subsets of pt
    lsd <- get_sim01_stan_data(d_all)
    
    lsd$ld$pri_a <- l_spec$prior$a
    lsd$ld$pri_b_trt <- l_spec$prior$b_trt
    lsd$ld$pri_b_reg <- l_spec$prior$b_reg
    lsd$ld$pri_b_loc <- l_spec$prior$b_loc
    
    # str(lsd$ld)
    foutname <- paste0(
      format(Sys.time(), format = "%Y%m%d%H%M%S"),
      "-sim-", ix, "-intrm-", l_spec$ia)
    
    
    
    # fit model - does it matter that I continue to fit the model after the
    # decision is made...?
    snk <- capture.output(
      f_1 <- m1$sample(
        lsd$ld, iter_warmup = 1000, iter_sampling = 2000,
        parallel_chains = 1, chains = 1, refresh = 0, show_exceptions = F,
        max_treedepth = 11,
        output_dir = output_dir_mcmc,
        output_basename = foutname
      )
    )
    
    log_info("Trial ", ix, " fitted models ", l_spec$ia)
    
    # extract posterior - marginal probability of outcome by group
    # dair vs rev
    d_post <- data.table(f_1$draws(
      variables = c(
        
        "p_1", "p_2", "rd", "rr"
        
      ),   
      format = "matrix"))
    
    if(return_posterior){
      d_post_all <- rbind(
        d_post_all,
        cbind(ia = l_spec$ia, d_post)
      )
    }
    
    d_post_long <- melt(d_post, measure.vars = names(d_post))
    
    d_post_smry_1[ia == l_spec$ia, 
                  mu := d_post_long[, mean(value), 
                                    keyby = .(variable)]$V1] 
    d_post_smry_1[ia == l_spec$ia, 
                  med := d_post_long[, median(value), 
                                     keyby = .(variable)]$V1]
    d_post_smry_1[ia == l_spec$ia, 
                  se := d_post_long[, sd(value), 
                                    keyby = .(variable)]$V1]
    d_post_smry_1[ia == l_spec$ia, 
                  q_025 := d_post_long[, quantile(value, prob = 0.025), 
                                       keyby = .(variable)]$V1]
    d_post_smry_1[ia == l_spec$ia, 
                  q_975 := d_post_long[, quantile(value, prob = 0.975), 
                                       keyby = .(variable)]$V1]
    
    log_info("Trial ", ix, " extracted posterior ", l_spec$ia)
    
    # superiority is implied by a high probability that the risk diff 
    # greater than zero
    pr_dec[l_spec$ia, "sup"]  <-   d_post[, mean(rd < l_spec$delta$sup)]
    # futility for the superiority decision is implied by a low probability 
    # that the risk diff is greater than some small value (2% difference)
    pr_dec[l_spec$ia, "fut"]  <-   d_post[, mean(rd < l_spec$delta$fut)]
    
    log_info("Trial ", ix, " calculated decision quantities ", l_spec$ia)
    
    decision[l_spec$ia, "sup"] <- pr_dec[l_spec$ia, "sup"] > l_spec$thresh$sup
    decision[l_spec$ia, "fut"] <- pr_dec[l_spec$ia, "fut"] < l_spec$thresh$fut
    
    log_info("Trial ", ix, " compared to thresholds ", l_spec$ia)
    
    # Earlier decisions are retained - once superiority has been decided, 
    # we retain this conclusion irrespective of subsequent post probabilities.
    # The following simply overwrites any decision reversal.
    # This means that there could be inconsistency with a silo pr_sup and the 
    # decision reported.
    decision[1:l_spec$ia, "sup"] <- apply(decision[1:l_spec$ia, "sup", drop = F], 2, function(z){ cumsum(z) > 0 })
    decision[1:l_spec$ia, "fut"] <- apply(decision[1:l_spec$ia, "fut", drop = F], 2, function(z){ cumsum(z) > 0 })
    
    # # superiority decisions apply to domains 1, 3 and 4
    # if(any(decision[l_spec$ia, "sup"])){
    #   # since there are only two treatments per cell, if a superiority decision 
    #   # is made then we have answered all the questions and we can stop 
    #   # enrolling into that cell. if there were more than two treatments then
    #   # we would need to take a different approach.
    #   
    #   dec_sup <- T
    #   
    #   
    # }
    # # stop enrolling if futile wrt superiority decision
    # if(any(decision[l_spec$ia,  "fut"])){
    #   dec_fut <- T
    # }
    
    # have we answered all questions of interest?
    if(
      # if rev (in late acute silo) is superior to dair (or superiority decision
      # decision is futile to pursue)
      (decision[l_spec$ia, "sup"] | decision[l_spec$ia, "fut"]) 
    ){
      log_info("Stop trial all questions addressed ", ix)
      stop_enrol <- T  
    } 
    
    log_info("Trial ", ix, " updated allocation control ", l_spec$ia)
    
    # next interim
    l_spec$ia <- l_spec$ia + 1
    
    if(l_spec$ia > N_analys){
      stop_enrol <- T  
    }
  }
  
  # did we stop (for any reason) prior to the final interim?
  stop_at <- N_analys
  
  # any na's in __any__ of the decision array rows means we stopped early:
  # we can just use sup to test this:
  if(any(is.na(decision[, "sup"]))){
    # interim where the stopping rule was met
    stop_at <- min(which(is.na(decision[, "sup"]))) - 1
    
    if(stop_at < N_analys){
      log_info("Stopped at analysis ", stop_at, " filling all subsequent entries")
      decision[(stop_at+1):N_analys, "sup"] <- decision[rep(stop_at, N_analys-stop_at), "sup"]
      decision[(stop_at+1):N_analys, "fut"] <- decision[rep(stop_at, N_analys-stop_at), "fut"]
    }
  }
  
  l_ret <- list(
    # data collected in the trial
    d_all = d_all[, .(y = sum(y), .N), keyby = .(ia, reg, loc, trt)],
    
    d_post_smry_1 = d_post_smry_1,
    
    pr_dec = pr_dec,
    decision = decision,
    
    stop_at = stop_at
  )
  
  if(return_posterior){
    l_ret$d_post_all <- copy(d_post_all)
  }
  # 
  
  return(l_ret)
}






run_sim01 <- function(){
  
  log_info(paste0(match.call()[[1]]))
  
  if(unname(Sys.info()[1]) == "Darwin"){
    log_info("On mac, reset cores to 5")
    g_cfgsc$mc_cores <- 5
  }
  
  l_spec <- list()
  # N by analysis
  l_spec$N <- g_cfgsc$N_pt
  l_spec$p_reg_alloc <- g_cfgsc$reg_alloc
  l_spec$p_loc_alloc_a <- g_cfgsc$loc_alloc_a
  l_spec$p_loc_alloc_d <- g_cfgsc$loc_alloc_d
  
  # model parameters
  # model params
  l_spec$bmu <- g_cfgsc$bmu
  l_spec$breg <- unlist(g_cfgsc$breg)
  l_spec$bloc <- unlist(g_cfgsc$bloc)
  l_spec$btrt <- unlist(g_cfgsc$btrt)
  
  l_spec$prior <- list()
  # location, scale
  l_spec$prior$a <- unlist(g_cfgsc$pri_a)
  l_spec$prior$b_reg <- unlist(g_cfgsc$pri_b_reg)
  l_spec$prior$b_loc <- unlist(g_cfgsc$pri_b_loc)
  l_spec$prior$b_trt <- unlist(g_cfgsc$pri_b_trt)
  
  l_spec$delta <- list()
  l_spec$delta$sup <- g_cfgsc$dec_delta_sup
  l_spec$delta$fut <- g_cfgsc$dec_delta_fut
  
  # domain specific
  l_spec$thresh <- list()
  l_spec$thresh$sup <- unlist(g_cfgsc$dec_thresh_sup)
  l_spec$thresh$fut <- unlist(g_cfgsc$dec_thresh_fut)
  
  l_spec$nex <- g_cfgsc$nex
  if(l_spec$nex > 0){
    log_info("Creating ", l_spec$nex, " example trials with full posterior")
    l_spec$ex_trial_ix <- sort(sample(1:g_cfgsc$nsim, size = l_spec$nex, replace = F))
  }
  return_posterior <- F
  str(l_spec)
  
  
  
  
  e = NULL
  log_info("Start simulation")
  r <- pbapply::pblapply(
    X=1:g_cfgsc$nsim, cl = g_cfgsc$mc_cores, FUN=function(ix) {
      # X=1:5, mc.cores = g_cfgsc$mc_cores, FUN=function(ix) {
      log_info("Simulation ", ix);
      
      if(ix %in% l_spec$ex_trial_ix){
        return_posterior = T  
      } else {
        return_posterior = F
      }
      
      ll <- tryCatch({
        run_trial(
          ix,
          l_spec,
          return_posterior = return_posterior
        )
      },
      error=function(e) {
        log_info("ERROR in MCLAPPLY LOOP (see terminal output):")
        message(" ERROR in MCLAPPLY LOOP " , e);
        log_info("Traceback (see terminal output):")
        message(traceback())
        stop(paste0("Stopping with error ", e))
      })
      
      # ll$d_all[, sum(N)]
      ll
    })
  
  
  
  log_info("Length of result set ", length(r))
  log_info("Sleep for 5 before processing")
  Sys.sleep(5)
  
  for(i in 1:length(r)){
    log_info("Element at index ",i, " is class ", class(r[[i]]))
    if(any(class(r[[i]]) %like% "try-error")){
      log_info("Element at index ",i, " has content ", r[[i]])  
    }
    log_info("Element at index ",i, " has names ", 
             paste0(names(r[[i]]), collapse = ", "))
  }
  
  
  d_pr_dec <- data.table()
  for(i in 1:length(r)){
    
    log_info("Appending pr_sup for result ", i)
    
    if(is.recursive(r[[i]])){
      d_pr_dec <- rbind(
        d_pr_dec,
        cbind(
          sim = i, ia = as.integer(rownames(r[[i]]$pr_dec)), r[[i]]$pr_dec
        ) 
      )  
    } else {
      log_info("Value for r at this index is not recursive ", i)
      message("r[[i]] ", r[[i]])
      message(traceback())
      stop(paste0("Stopping due to non-recursive element "))
    }
    
  }
  
  # decisions based on probability thresholds
  d_decision <- rbind(
    
    data.table(do.call(rbind, lapply(1:length(r), function(i){ 
      m <- r[[i]]$decision[, "sup"]
      cbind(sim = i, ia = 1:length(m), quant = "sup", value = m) 
    } ))), 
    
    data.table(do.call(rbind, lapply(1:length(r), function(i){ 
      m <- r[[i]]$decision[, "fut"]
      cbind(sim = i, ia = 1:length(m), quant = "fut", value = m) 
    } )))
    
  )
  
  d_decision[, `:=`(sim = as.integer(sim), 
                    ia = as.integer(ia),
                    value = as.logical(value)
  )]
  
  d_post_smry_1 <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_post_smry_1
  } ), idcol = "sim")
  
  
  # data from each simulated trial
  d_all <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_all
  } ), idcol = "sim")
  
  
  
  
  d_post_all <- data.table(do.call(rbind, lapply(1:length(r), function(i){
    # if the sim contains full posterior (for example trial) then return
    if(!is.null(r[[i]]$d_post_all)){
      cbind(sim = i, r[[i]]$d_post_all)
    }
    
  } )))
  
  l <- list(
    cfg = g_cfgsc,
    
    
    
    d_pr_dec = d_pr_dec, 
    
    d_decision = d_decision,
    d_post_smry_1 = d_post_smry_1,
    
    d_all = d_all,
    
    d_post_all = d_post_all
  )
  
  log_info("Command line arguments ", paste(args[2], collapse = ", "))
  
  toks <- unlist(tstrsplit(args[2], "[-.]"))
  log_info("Tokens ", paste(toks, collapse = ", "))
  
  fname <- paste0("data/sim01/sim01-", toks[5],  "-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".qs")
  
  log_info("Saving results file ", fname)
  
  qs::qsave(l, file = fname)
}

run_none_sim01 <- function(){
  log_info("run_none_sim01: Nothing doing here bud.")
}

main_sim01 <- function(){
  funcname <- paste0(args[1], "()")
  log_info("Main, invoking ", funcname)
  eval(parse(text=funcname))
}

main_sim01()


