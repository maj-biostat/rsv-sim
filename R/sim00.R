
source("./R/init.R")
source("./R/data.R")
# same data generating process
source("./R/data-sim01.R")

# Command line arguments list scenario (true dose response),
# the number of simulations to run, the number of cores
# to use and the simulator to use.
args = commandArgs(trailingOnly=TRUE)

# Load cfg based on cmd line args.
if (length(args)<1) {
  log_info("Setting default run method (does nothing)")
  args[1] = "run_none_sim00"
  args[2] = "./sim00/cfg-sim00-v01.yml"
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

output_dir_mcmc <- paste0(getwd(), "/tmp")

ix <- 1

# Main trial loop.
run_trial00 <- function(
    ix,
    l_spec,
    return_posterior = F
){
  
  log_info("Entered  run_trial for trial ", ix)
  
  # loop controls
  N_total <- sum(l_spec$N)
  
  if(return_posterior){
    d_post_all <- data.table()
  }
  
  l_spec$is <- 1
  l_spec$ie <- N_total
  
  l_spec$N <- N_total
  
  # generate the enrolment times for the entire trial sample one time 
  # and then index as necessary
  l_spec$t0 <- c(0, cumsum(rexp(l_spec$N-1, l_spec$recruit_rate)))
  
  # can use the same data generation process as sim01
  d <- get_sim01_trial_data(l_spec)
  # d[]
  
  # ggplot(d, aes(x = t0/365, y = id)) +
  #   geom_step() +
  #   theme_bw()
  # 
  
  pri_a <- l_spec$prior$p[1]
  pri_b <- l_spec$prior$p[2]
  d_post <- data.table(
    p_1 = rbeta(2000, pri_a + d[trt == 1, sum(y)], pri_b + d[trt == 1, .N] - d[trt == 1, sum(y)]),
    p_2 = rbeta(2000, pri_a + d[trt == 2, sum(y)], pri_b + d[trt == 2, .N] - d[trt == 2, sum(y)])
  )
  d_post[, rd_2_1 := p_2 - p_1]
  d_post <- melt(d_post, measure.vars = names(d_post), variable.name = "par")
  
  # ggplot(d_post, aes(x = value)) +
  #   geom_density() +
  #   ggh4x::facet_wrap2(~par, nrow = 2, axes = "x") +
  #   theme_bw()
  #
  
  
  d_post_smry_1 <- d_post[, .(
    mu = mean(value),
    med = median(value),
    se = sd(value), 
    q_025 = quantile(value, prob = 0.025), 
    q_975 = quantile(value, prob = 0.975)
  ), keyby = .(par)]
  
  d_pr_dec <- rbind(
    d_post[par %like% "rd_", .(
      rule = "sup",
      p = mean(value < l_spec$delta$sup)
    ), keyby = .(par)],
    d_post[par %like% "rd_", .(
      rule = "fut",
      p = mean(value < l_spec$delta$fut)
    ), keyby = .(par)]
  )
  d_pr_dec[rule == "sup", dec := as.integer(p > l_spec$thresh$sup)]
  d_pr_dec[rule == "fut", dec := as.integer(p < l_spec$thresh$fut)]
  
  l_ret <- list(
    # data collected in the trial
    d_all = d[, .(ty = max(ty), y = sum(y), .N), keyby = .(ia, reg, loc, trt)],
    
    d_post_smry_1 = d_post_smry_1,
    
    d_pr_dec = d_pr_dec
  )
  
  if(return_posterior){
    l_ret$d_post_all <- copy(d_post)
  }
  # 
  
  l_ret
  
  
}

run_sim00 <- function(){
  
  log_info(paste0(match.call()[[1]]))
  
  if(unname(Sys.info()[1]) == "Darwin"){
    log_info("On mac, reset cores to 5")
    g_cfgsc$mc_cores <- 5
  }
  
  l_spec <- list()
  l_spec$desc <- g_cfgsc$desc
  l_spec$nsim <- g_cfgsc$nsim
  # N by analysis
  l_spec$N <- g_cfgsc$N_pt
  l_spec$recruit_rate <- g_cfgsc$recruit_rate
  l_spec$fu_days <- g_cfgsc$fu_days
  l_spec$p_reg_alloc <- g_cfgsc$reg_alloc
  l_spec$p_loc_alloc_a <- g_cfgsc$loc_alloc_a
  l_spec$p_loc_alloc_d <- g_cfgsc$loc_alloc_d
  
  # model parameters
  # model params
  l_spec$bmu <- g_cfgsc$bmu
  l_spec$breg <- unlist(g_cfgsc$breg)
  l_spec$bloc <- unlist(g_cfgsc$bloc)
  l_spec$btrt <- unlist(g_cfgsc$btrt)
  
  # priors on beta binomial
  l_spec$prior <- list()
  # location, scale
  l_spec$prior$p <- unlist(g_cfgsc$pri_p)
  
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
  
  
  
  ## LOOP -------
  
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
        run_trial00(
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
  
  d_pr_dec <- rbindlist(lapply(1:length(r), function(i){ 
    r[[i]]$d_pr_dec
  } ), idcol = "sim")
  
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
    d_post_smry_1 = d_post_smry_1,
    d_all = d_all,
    d_post_all = d_post_all
  )
  
  log_info("Command line arguments ", paste(args[2], collapse = ", "))
  
  toks <- unlist(tstrsplit(args[2], "[-.]"))
  log_info("Tokens ", paste(toks, collapse = ", "))
  
  fname <- paste0("data/sim00/sim00-", toks[5],  "-", format(Sys.time(), "%Y%m%d-%H%M%S"), ".qs")
  
  log_info("Saving results file ", fname)
  
  qs::qsave(l, file = fname)
  
  
  
}

run_none_sim00 <- function(){
  log_info("run_none_sim00: Nothing doing here bud.")
}

main_sim00 <- function(){
  funcname <- paste0(args[1], "()")
  log_info("Main, invoking ", funcname)
  eval(parse(text=funcname))
}

main_sim00()

