
source("./R/init.R")
source("./R/data.R")
source("./R/data-sim02.R")

# Command line arguments list scenario (true dose response),
# the number of simulations to run, the number of cores
# to use and the simulator to use.
args = commandArgs(trailingOnly=TRUE)

# Load cfg based on cmd line args.
if (length(args)<1) {
  log_info("Setting default run method (does nothing)")
  args[1] = "run_none_sim02"
  args[2] = "./sim02/cfg-sim02-v01.yml"
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
m1 <- cmdstanr::cmdstan_model("stan/sim02-v01.stan")

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
  
  lambda = l_spec$pt_per_day
  # ramp up over x days 
  rho = function(t) pmin(t/l_spec$ramp_up_days, 1)
  
  loc_t0 <- get_enrol_time(sum(l_spec$N), lambda, rho)
  # follow up time - deterministic based on follow up length
  loc_t <- loc_t0 + l_spec$fu_days
  
  # lambda = 0.75
  # # ramp up over x months
  # rho = function(t) pmin(t/120, 1)
  # 
  # rr <- unlist(pblapply(1:100, cl = 4, FUN=function(ii){
  #   ttt <- get_enrol_time(1000, lambda, rho)
  #   max(ttt)
  # }))
  # mean(rr) / 365
  
  
  # loop controls
  stop_enrol <- FALSE
  l_spec$ic <- 1 # interim number
  N_analys <- length(l_spec$N)
  
  
  # posterior summaries
  g_par_p = c("p_1", "p_2")
  g_par_rd = c("rd_2_1")
  d_post_smry_1 <- CJ(
    ic = 1:N_analys,
    par = factor(c(g_par_p, g_par_rd))
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
  g_rule_type <- c("sup", "fut")
  
  d_pr_dec <- CJ(
    ic = 1:N_analys,
    rule = factor(g_rule_type),
    par = factor(g_par_rd),
    p = NA_real_,
    dec = NA_integer_
  )
  
  # store all simulated trial pt data
  d_all <- data.table()
  
  if(return_posterior){
    d_post_all <- data.table()
    d_freq <- data.table()
  }
  
  
  ## LOOP -------
  while(!stop_enrol){
    
    log_info("Trial ", ix, " cohort ", l_spec$ic)
    
    # next chunk of data on pts.
    if(l_spec$ic == 1){
      # starting pt index in data
      l_spec$is <- 1
      l_spec$ie <- l_spec$is + l_spec$N[l_spec$ic] - 1
    } else {
      l_spec$is <- l_spec$ie + 1
      l_spec$ie <- l_spec$is + l_spec$N[l_spec$ic] - 1
    }
    
    # id and time
    l_spec$t0 = loc_t0[l_spec$is:l_spec$ie]
    l_spec$tfu = loc_t[l_spec$is:l_spec$ie]
    
    # We are assuming that the analysis takes place on pt having reached endpoint
    
    d <- get_sim02_trial_data(l_spec)
    
    log_info("Trial ", ix, " cohort ", l_spec$ic, " generated")
    
    # combine the existing and new data
    d_all <- rbind(d_all, d)
    
    # select subset of pt reaching fu to use in analyses
    if(l_spec$ie == sum(l_spec$N)){
      log_info("Trial ", ix, " final analysis, using all pt")
      t_now <- d_all[, max(tfu)]
      d_mod <- copy(d_all)
      
      d_all[is.na(ia) & id %in% d_mod$id, ia := l_spec$ic]
      d_all[is.na(t_anlys) & id %in% d_mod$id, t_anlys := t_now]
      
    } else {
      t_now <- d_all[, max(t0)]
      d_mod <- d_all[tfu <= t_now]
      
      # set analysis set indicator for when pt entered analyses
      d_all[is.na(ia) & id %in% d_mod$id, ia := l_spec$ic]
      d_all[is.na(t_anlys) & id %in% d_mod$id, t_anlys := t_now]
    }
    
    
    # create stan data format based on the relevant subsets of pt
    lsd <- get_sim02_stan_data(d_mod)
    
    lsd$ld$pri_b_0 <- l_spec$prior$b_0
    lsd$ld$pri_b_trt <- l_spec$prior$b_trt
    lsd$ld$pri_b_reg <- l_spec$prior$b_reg
    lsd$ld$pri_b_loc <- l_spec$prior$b_loc
    
    # str(lsd$ld)
    foutname <- paste0(
      format(Sys.time(), format = "%Y%m%d%H%M%S"),
      "-sim-", ix, "-intrm-", l_spec$ia)
    
    
    snk <- capture.output(
      f_1 <- m1$sample(
        lsd$ld, iter_warmup = 1000, iter_sampling = 2000,
        parallel_chains = 1, chains = 1, refresh = 0, show_exceptions = F,
        max_treedepth = 11,
        output_dir = output_dir_mcmc,
        output_basename = foutname
      )
    )
    
    log_info("Trial ", ix, " fitted models ", l_spec$ic)
    
    # extract posterior - marginal probability of outcome by group
    # dair vs rev
    d_post <- data.table(f_1$draws(
      variables = c(g_par_p, g_par_rd),   
      format = "matrix"))
    
    if(return_posterior){
      d_tmp <- data.table(f_1$draws(
        variables = c(g_par_p, g_par_rd,
                      "b_0", "b_trt", "b_reg", "b_loc") ,   
        format = "matrix"))
      
      d_post_all <- rbind(
        d_post_all,
        cbind(ic = l_spec$ic, d_tmp)
      )
      # colMeans(d_post_all)
      
      d_tmp <- copy(lsd$dd)
      d_tmp[, reg := factor(reg)]
      d_tmp[, loc := factor(loc)]
      d_tmp[, trt := factor(trt)]
      
      f_2 <- glm(cbind(y, n-y) ~ trt + reg + loc, data = d_tmp, family = binomial)
      s <- summary(f_2)
      
      pred_f_2 <- data.table(avg_predictions(f_2, by = "trt"))
      pred_f_2[, s.value := NULL]
      pred_f_2[, df := NULL]
      setnames(pred_f_2, "trt", "par")
      pred_f_2[, par := factor(par, levels = 1:2, labels = c("p_1", "p_2"))]
      
      d_tmp <- data.table(
        par = rownames(s$coef), 
        s$coef, 
        # no profiling
        confint.default(f_2, type = "Wald"))
      colnames(d_tmp) <- c(
        "par", "estimate", "std.error", "statistic", 
        "p.value", "conf.low", "conf.high")
      d_tmp[par == "(Intercept)", par := "b_0"]
      d_tmp[par == "trt2", par := "b_trt[2]"]
      d_tmp[par == "reg2", par := "b_reg[2]"]
      d_tmp[par == "loc2", par := "b_loc[2]"]
      
      d_tmp <- rbind(
        d_tmp, pred_f_2
        )
      
      d_freq <- rbind(
        d_freq,
        cbind(ic = l_spec$ic, d_tmp)
      )
      
    }
    
    d_post_long <- melt(d_post, measure.vars = names(d_post), variable.name = "par")
    
    # sanity check
    # d_fig <- copy(d_post_long)
    # d_fig[par %like% "p", par_type := "prob"]
    # d_fig[par %like% "rd", par_type := "rd"]
    # d_fig_2 <- d_all[, .(value = mean(y)), keyby = trt]
    # setnames(d_fig_2, "trt", "par")
    # p1 <- ggplot(d_fig[par %like% "p", ],
    #              aes(x = value, group = par, col = par)) +
    #   geom_vline(
    #     data = d_all[, .(value = mean(y)), keyby = trt],
    #     aes(xintercept = value), lwd = 0.3
    #   ) +
    #   geom_density()
    # p2 <- ggplot(d_fig[par %in% c("rd_2_1", "rd_3_1"), ],
    #              aes(x = value, group = par, col = par)) +
    #   geom_density()
    # p1/p2
    
    # merge posterior summaries for current interim
    d_post_smry_1[
      d_post_long[, .(ic = l_spec$ic,
                      mu = mean(value),
                      med = median(value),
                      se = sd(value),
                      q_025 = quantile(value, prob = 0.025),
                      q_975 = quantile(value, prob = 0.975)
      ), keyby = par], 
      on = .(ic, par), `:=`(
        mu = i.mu,
        med = i.med,
        se = i.se,
        q_025 = i.q_025,
        q_975 = i.q_975
      )]
    
    log_info("Trial ", ix, " extracted posterior ", l_spec$ic)
    
    # compute and merge the current probability and decision trigger status
    
    # The trial setup is such that if the intervention were beneficial, it would 
    # reduce the occurrence of medically attended LRIs.
    # Therefore, when we compare a beneficial treatment with the soc, the central
    # value of the posterior would be negative and ideally the whole posterior
    # would be below zero.
  
    d_pr_dec[
      
      rbind(
        
        # superiority decision is implied by a high probability 
        # that the risk diff is below zero e.g. Pr(RD < 0) is high
        
        d_post_long[par %in% g_par_rd, .(
          ic = l_spec$ic,
          rule = factor("sup", levels = g_rule_type),
          p = mean(value < l_spec$delta$sup),
          dec = as.integer(mean(value < l_spec$delta$sup) > l_spec$thresh$sup)
        ), keyby = par],
        
        # futility for the superiority decision is implied by a low probability 
        # that the risk diff is lower than some small value 
        # e.g. Pr(RD < -0.02) is low
        
        d_post_long[par %in% g_par_rd, .(
          ic = l_spec$ic,
          rule = factor("fut", levels = g_rule_type),
          p = mean(value < l_spec$delta$fut),
          dec = as.integer(mean(value < l_spec$delta$fut) < l_spec$thresh$fut)
        ), keyby = par]
      ),
      
      on = .(ic, rule, par), `:=`(
        p = i.p, dec = i.dec  
      )
    ]
    
    
    # have we answered all questions of interest?
    d_stop <- d_pr_dec[
      ic <= l_spec$ic & par %in% c("rd_2_1"), 
      .(resolved = as.integer(sum(dec) > 0)), keyby = .(par)]
    
    if(d_stop[par == "rd_2_1", resolved]) {
      log_info("Stop trial all questions addressed ", ix)
      stop_enrol <- T  
    } 
    
    log_info("Trial ", ix, " updated allocation control ", l_spec$ia)
    
    # next interim
    l_spec$ic <- l_spec$ic + 1
    
    if(l_spec$ic > N_analys){
      stop_enrol <- T  
    }
  }
  
  # did we stop (for any reason) prior to the final interim?
  stop_at <- l_spec$ic - 1
  
  l_ret <- list(
    # data collected in the trial
    # d_all = d_all[, .(ty = max(ty), y = sum(y), .N), keyby = .(ia, reg, loc, trt)],
    
    d_all = d_all,
    d_post_smry_1 = d_post_smry_1,
    d_pr_dec = d_pr_dec,
    stop_at = stop_at
  )
  
  if(return_posterior){
    l_ret$d_post_all <- copy(d_post_all)
    l_ret$d_freq <- copy(d_freq)
  }
  # 
  
  return(l_ret)
}






run_sim02 <- function(){
  
  log_info(paste0(match.call()[[1]]))
  
  if(unname(Sys.info()[1]) == "Darwin"){
    log_info("On mac, reset cores to 5")
    g_cfgsc$mc_cores <- 5
  }
  
  l_spec <- list()
  l_spec$desc <- g_cfgsc$desc
  l_spec$nsim <- g_cfgsc$nsim
  l_spec$nex <- g_cfgsc$nex
  
  # N by analysis
  l_spec$N <- g_cfgsc$N_pt
  
  l_spec$pt_per_day <- g_cfgsc$pt_per_day
  l_spec$ramp_up_days <-  g_cfgsc$ramp_up_days
  
  l_spec$fu_days <- g_cfgsc$fu_days
  l_spec$fu_lab <- g_cfgsc$fu_lab
  
  # distribution of demographics (region and locality)
  l_spec$p_reg_alloc <- g_cfgsc$reg_alloc
  l_spec$p_loc_alloc_a <- g_cfgsc$loc_alloc_a
  l_spec$p_loc_alloc_d <- g_cfgsc$loc_alloc_d
  
  # model parameters
  # model params
  l_spec$b_0 <- g_cfgsc$b_0
  l_spec$b_reg <- unlist(g_cfgsc$b_reg)
  l_spec$b_loc <- unlist(g_cfgsc$b_loc)
  l_spec$b_trt <- unlist(g_cfgsc$b_trt)
  
  l_spec$prior <- list()
  # location, scale
  l_spec$prior$b_0 <- unlist(g_cfgsc$pri_b_0)
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
  
  # compute marginal risk by group
  # probability of being in each combination
  # P(A, U) = P(A)P(U | A) etc
  # rand is 1:1 so each group has the same weights
  l_spec$p_alice_urban <- l_spec$p_reg_alloc * (1-l_spec$p_loc_alloc_a)
  l_spec$p_alice_remote <- l_spec$p_reg_alloc * l_spec$p_loc_alloc_a
  l_spec$p_darwin_urban <- (1-l_spec$p_reg_alloc) * (1-l_spec$p_loc_alloc_d)
  l_spec$p_darwin_remote <- (1-l_spec$p_reg_alloc) * l_spec$p_loc_alloc_d
  
  # weight combinations by probability of being in strata 
  l_spec$p_y_ctl <- l_spec$p_alice_urban * (l_spec$b_0 + l_spec$b_reg[1] + l_spec$b_loc[1]) + 
    l_spec$p_alice_remote * (l_spec$b_0 + l_spec$b_reg[1] + l_spec$b_loc[2])  + 
    l_spec$p_darwin_urban * (l_spec$b_0 + l_spec$b_reg[2] + l_spec$b_loc[1])  + 
    l_spec$p_darwin_remote * (l_spec$b_0 + l_spec$b_reg[2] + l_spec$b_loc[2]) 
  
  l_spec$p_y_trt <- l_spec$p_alice_urban * (l_spec$b_0 + l_spec$b_reg[1] + l_spec$b_loc[1] + l_spec$b_trt[2]) + 
    l_spec$p_alice_remote * (l_spec$b_0 + l_spec$b_reg[1] + l_spec$b_loc[2] + l_spec$b_trt[2])  + 
    l_spec$p_darwin_urban * (l_spec$b_0 + l_spec$b_reg[2] + l_spec$b_loc[1] + l_spec$b_trt[2])  + 
    l_spec$p_darwin_remote * (l_spec$b_0 + l_spec$b_reg[2] + l_spec$b_loc[2] + l_spec$b_trt[2]) 
  
  if(l_spec$nex > 0){
    log_info("Creating ", pmin(l_spec$nex, g_cfgsc$nsim), " example trials with full posterior")
    l_spec$ex_trial_ix <- sort(sample(1:g_cfgsc$nsim, size = pmin(l_spec$nex, g_cfgsc$nsim), replace = F))
  }
  return_posterior <- T
  str(l_spec)
  
  e = NULL
  log_info("Start simulation")
  r <- pbapply::pblapply(
    X=1:g_cfgsc$nsim, cl = g_cfgsc$mc_cores, FUN=function(ix) {
      # X=1:5, cl = g_cfgsc$mc_cores, FUN=function(ix) {
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
          sim = i, r[[i]]$d_pr_dec
        ) 
      )   
    } else {
      log_info("Value for r at this index is not recursive ", i)
      message("r[[i]] ", r[[i]])
      message(traceback())
      stop(paste0("Stopping due to non-recursive element "))
    }
    
  }
  
  
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
  
  d_freq <- data.table(do.call(rbind, lapply(1:length(r), function(i){
    # if the sim contains full posterior (for example trial) then return
    if(!is.null(r[[i]]$d_freq)){
      cbind(sim = i, r[[i]]$d_freq)
    }
    
  } )))
  
  l <- list(
    cfg = l_spec,
    d_pr_dec = d_pr_dec, 
    d_post_smry_1 = d_post_smry_1,
    d_all = d_all,
    d_post_all = d_post_all,
    d_freq = d_freq
  )
  
  log_info("Command line arguments ", paste(args[2], collapse = ", "))
  
  toks <- unlist(tstrsplit(args[2], "[-.]"))
  log_info("Tokens ", paste(toks, collapse = ", "))
  
  fname <- paste0("data/sim02/sim02-", toks[5],  "-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".qs")
  
  log_info("Saving results file ", fname)
  
  qs::qsave(l, file = fname)
}

run_none_sim02 <- function(){
  log_info("run_none_sim02: Nothing doing here bud.")
}

main_sim02 <- function(){
  funcname <- paste0(args[1], "()")
  log_info("Main, invoking ", funcname)
  eval(parse(text=funcname))
}

main_sim02()


