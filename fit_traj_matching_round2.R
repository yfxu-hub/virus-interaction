# ---------------------------------------------------------------------------------------------------------------------
# Code to fit DETERMINISTIC flu/RSV interaction model
# Round 2: Fit all seasons simultaneously, constraining shared parameters to be the same for all seasons
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Set seed:
set.seed(749501349)

# Load libraries:
library(nloptr)
library(parallel)
library(doMC)
library(tidyverse)
library(pomp)
library(viridis)
library(gridExtra)
library(zoo)
library(readxl)
#library(doParallel)

# Get cluster environmental variables:
jobid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")); print(jobid)
no_jobs <- as.integer(Sys.getenv("NOJOBS")); print(no_jobs)
sobol_size <- as.integer(Sys.getenv("SOBOLSIZE")); print(sobol_size)
which_round <- as.integer(Sys.getenv("WHICHROUND")); print(which_round)
search_type <- as.character(Sys.getenv("SEARCHTYPE")); print(search_type)
int_eff <- as.character(Sys.getenv("INTERACTIONEFFECT")); print(int_eff)
prof_lik <- as.logical(Sys.getenv("PROFLIK")); print(prof_lik)
# prof_val <- as.numeric(as.character(Sys.getenv("PROFVAL"))); print(prof_val)
sens <- as.character(Sys.getenv("SENS")); print(sens)
fit_canada <- as.logical(Sys.getenv("FITCANADA")); print(fit_canada)
run_parallel <- as.logical(Sys.getenv("RUNPARALLEL")); print(run_parallel)

# # Set parameters for local run:
# jobid <-1
# no_jobs <- 1
# sobol_size <- 2
# which_round <- 2
# search_type <- 'round2_CIs'
# int_eff <- 'susc' # 'susc', 'sev', or 'both' - fit impact of interaction on susceptibility or severity, or both?
# prof_lik <- FALSE
# sens <- 'sinusoidal_forcing' # 'main', 'less_circ_h3', 'sinusoidal_forcing', 'no_ah', 'no_int', 'no_rsv_immune', 'h3_covar', 'rhino_covar'
# fit_canada <- FALSE
# run_parallel <- TRUE

# Set parameter values:
vir1_name<-'covid'
vir2_name<-'h1'

# Set parameters for run:
#outbreak<-read_excel("/home/yfxu/virus_interference/data/multi_pathogons_outbreak.xlsx")%>%as.data.frame()
outbreak<-read_excel("/home/yfxu/virus_interference/data/multi_pathogons_outbreak_1.xlsx")%>%as.data.frame()
match_v<-intersect(grep(vir1_name, outbreak$season, value = TRUE),grep(vir2_name, outbreak$season, value = TRUE))
outbreak<-outbreak[which(outbreak$season %in% match_v),]
seasons <- outbreak$season

vir1<-vir1_name
vir2<-vir2_name

#seasons<-c('s22-23','s23-24','s24-25')
#seasons<-c("covid-h1")

debug_bool <- FALSE
# if (sens == 'less_circ_h3') {
#   seasons <- c('s17-18', 's18-19')
# }
# if (fit_canada) {
#   seasons <- c('s10-11', 's11-12', 's12-13', 's13-14')
# }

if (sens == 'sinusoidal_forcing') {
  time_max <- 23.75 # Maximal execution time (in hours)
} else {
  time_max <- 14.75 # Maximal execution time (in hours)
}

Ri_max1 <- 4.0
Ri_max2 <- 3.0
d2_max <- 10.0

if (prof_lik) {
  
  jobid_orig <- ceiling(jobid / no_jobs)
  jobid <- (jobid - 1) %% no_jobs + 1
  
  print(jobid_orig)
  print(jobid)
  
  prof_param <- 'theta_lambda1'
  # prof_param <- 'theta_lambda2'
  
  # prof_val <- seq(0.0, 0.2, by = 0.01)[jobid_orig]
  prof_val <- seq(0, 0.02, by = 0.001)[jobid_orig]
  # prof_val <- seq(0, 1.0, by = 0.05)[jobid_orig]
  # prof_val <- seq(0.6, 0.7, by = 0.005)[jobid_orig]
  print(prof_val)
  
} else {
  jobid_orig <- jobid
}

# # Fit for synthetic data from age-structured model?:
# age_structured <- TRUE

# ---------------------------------------------------------------------------------------------------------------------

# Functions

create_obj_fxn <- function(po, estpars) {
  # Creates the objective function for a given season
  # param po: A pomp model object for a specific season
  # param estpars: A vector listing the parameters to be fit
  # returns: The negative log likelihood
  
  ofun <- traj_objfun(data = po, 
                      est = estpars, 
                      partrans = po@partrans,
                      verbose=TRUE)
  
  return(ofun)
}



calculate_global_loglik <- function(trans_vals) {
  # Calculates the log-likelihood for each season, and combines to yield a global log-likelihood
  # param trans_vals: Unnamed vector of transformed parameters; fxn only works if this is the only input?
  # returns: The global, negative log-likelihood
  
  # Add names to vector:
  if(is.null(names(trans_vals))) names(trans_vals) <- x0_trans_names
  
  # Split params into shared and unit params:
  unit_in <- as.data.frame(sapply(paste0('^', seasons, '_'), grep, x = names(trans_vals)))
  if (ncol(unit_in) == 1) {
    unit_in <- as.data.frame(matrix(data = c(unlist(unit_in)), nrow = 1))
  }
  names(unit_in) <- seasons
  shared_params <- trans_vals[-unique(unlist(unit_in))]
  
  unit_params <- list()
  for (i in 1:length(seasons)) {
    unit <- trans_vals[unit_in[, i]]
    if (length(unit) > 0) {
      names(unit) <- str_split(names(unit), '_', simplify = TRUE)[, 2]
      unit_params[[i]] <- unit
    }
  }
  
  # Get -ll for each season:
  units_ll <- rep(NA, length(seasons))
  
  for (i in 1:length(seasons)) {
    if (length(unit_params) > 0) {
      params_temp <- c(shared_params, unit_params[[i]])
    } else {
      params_temp <- c(shared_params)
    }
    
    units_ll[i] <- obj_fun_list[[i]](params_temp)
  }
  
  # Calculate global -ll:
  glob_ll <- sum(units_ll)
  return(glob_ll)
}

transform_params <- function(orig_vals, po, seas, params_all, params_shared) {
  # Transforms parameters as needed
  # param orig_vals: Untransformed parameter values
  # param po: pomp object with correct parameter transformations
  # param seas: Vector containing all seasons of interest
  # param params_all: Names of all parameters to be estimated
  # param params_shared: Names of the shared parameters to be estimated
  # returns: Vector of transformed parameter values
  
  names(orig_vals) <- params_all
  
  orig_vals_shared <- orig_vals[which(params_shared %in% params_all)]
  
  coef(po, params_shared) <- orig_vals_shared
  trans_vals_shared <- coef(po, params_shared, transform = TRUE)
  trans_vals <- trans_vals_shared
  
  po_save <- po
  
  for (i in 1:length(seas)) {
    po <- po_save
    
    orig_vals_unit <- orig_vals[grep(paste0('^', seas[i], '_'), params_all, value = TRUE)]
    names(orig_vals_unit) <- gsub(paste0(seas[i], '_'), '', names(orig_vals_unit))
    
    coef(po, names(orig_vals_unit)) <- orig_vals_unit
    trans_vals_unit <- coef(po, names(orig_vals_unit), transform = TRUE)
    names(trans_vals_unit) <- paste0(seas[i], '_', names(orig_vals_unit))
    trans_vals <- c(trans_vals, trans_vals_unit)
  }
  
  return(trans_vals)
}

back_transform_params <- function(trans_vals, po, seas, params_all, params_shared) {
  # Un-transforms parameters as needed
  # param orig_vals: Transformed parameter values
  # param po: pomp object with correct parameter 
  # param seas: Vector containing all seasons of interest
  # param params_all: Names of all parameters to be estimated
  # param params_shared: Names of the shared parameters to be estimated
  # returns: Vector of un-transformed parameter values
  
  names(trans_vals) <- params_all
  
  trans_vals_shared <- trans_vals[which(params_shared %in% params_all)]
  
  coef(po, params_shared, transform = TRUE) <- trans_vals_shared
  orig_vals_shared <- coef(po, params_shared)
  orig_vals <- orig_vals_shared
  
  po_save <- po
  
  for (i in 1:length(seas)) {
    po <- po_save
    
    trans_vals_unit <- trans_vals[grep(paste0('^', seas[i], '_'), params_all, value = TRUE)]
    names(trans_vals_unit) <- gsub(paste0(seas[i], '_'), '', names(trans_vals_unit))
    
    coef(po, names(trans_vals_unit), transform = TRUE) <- trans_vals_unit
    orig_vals_unit <- coef(po, names(trans_vals_unit))
    
    names(orig_vals_unit) <- paste0(seas[i], '_', names(trans_vals_unit))
    orig_vals <- c(orig_vals, orig_vals_unit)
  }
  
  return(orig_vals)
}

# ---------------------------------------------------------------------------------------------------------------------

# Fit using trajectory matching

# Loop through years and construct pomp models:
po_list <- vector('list', length(seasons))
for (yr_index in 1:length(seasons)) {
  yr <- seasons[yr_index]
  print(yr)
  
  # Load data and create pomp object:
  source('/home/yfxu/virus_interference/code/resp_interaction_model.R')
  
  # Check whether any data for given season:
  if (exists('resp_mod')) {
    
    # If doing profile likelihood, set interaction parameter:
    if (prof_lik) {
      coef(resp_mod, prof_param) <- prof_val
    }
    
    # Add pomp object to list:
    po_list[[yr_index]] <- resp_mod
    
  }
  
  # Remove pomp object before repeating loop:
  rm(resp_mod)
  
}

# Check that there are no empty elements:
expect_true(all(!lapply(po_list, length) == 0))

# Choose parameters to estimate:
if (int_eff == 'susc') {
  if (prof_lik) {
    
    if (sens == 'sinusoidal_forcing') {
      
      shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                          'alpha', 'phi', 'b1', 'b2', 'phi1', 'phi2')
      
    } else {
      
      shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                          'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
      
    }
    
    shared_estpars <- shared_estpars[shared_estpars != prof_param]
  } else {
    
    if (sens == 'sinusoidal_forcing') {
      shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                          'alpha', 'phi', 'b1', 'b2', 'phi1', 'phi2')
      
    } else if (sens == 'h3_covar') {
      shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                          'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2',
                          'beta_h3')
    } else if (sens == 'rhino_covar') {
      shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                          'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2',
                          'beta_rhino')
    } else {
      shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                          'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
    }
    
  }
} else if (int_eff == 'both') {
  shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'theta_rho1', 'theta_rho2',
                      'delta1', 'd2', 'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
} else {
  stop('Unrecognized int_eff value.')
}

unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')
if (sens == 'no_rsv_immune') {
  unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10')
}

unit_sp_estpars <- c()
for (i in 1:length(seasons)) {
  unit_sp_estpars <- c(unit_sp_estpars, paste(seasons[i], unit_estpars, sep = '_'))
}
rm(i)

true_estpars <- c(shared_estpars, unit_estpars)
estpars <- c(shared_estpars, unit_sp_estpars)

# Set upper/lower values for global params:
start_range <- data.frame(rho1 = c(0, 0.9999), #c(0, 1.0)
                          rho2 = c(0, 0.9999),
                          theta_lambda1 = c(0, 0.9999),
                          theta_lambda2 = c(0, 0.9999),
                          theta_rho1 = c(0, 5.0),
                          theta_rho2 = c(0, 5.0),
                          delta1 = c(7 / 60, 7),
                          d2 = c(0, 10),
                          alpha = c(0, 0.5),
                          phi = c(0, 52.25),
                          eta_temp1 = c(-0.5, 0.5),
                          eta_temp2 = c(-0.5, 0.5),
                          eta_ah1 = c(-0.5, 0.5),
                          eta_ah2 = c(-0.5, 0.5),
                          b1 = c(0.05, 0.2),
                          b2 = c(0.05, 0.2),
                          phi1 = c(0, 52.25),
                          phi2 = c(0, 52.25),
                          beta_h3 = c(0, 5.0),
                          beta_rhino = c(0, 5.0))

# Set upper/lower values for unit params (broad):
unit_start_range <- data.frame(Ri1 = c(1.0, Ri_max1),
                               Ri2 = c(2.0, Ri_max2),
                               I10 = c(0, 1e-3),
                               I20 = c(0, 1e-3),
                               R10 = c(0, 0.3),
                               R20 = c(0, 0.3),
                               R120 = c(0, 0.3))
if (sens == 'no_rsv_immune') {
  unit_start_range <- data.frame(Ri1 = c(1.0, Ri_max1),
                                 Ri2 = c(2.0, Ri_max2),
                                 I10 = c(0, 1e-3),
                                 I20 = c(0, 1e-3),
                                 R10 = c(0, 0.3))
}

#Get 95% CI from round 1 for unit params:
#tj_res_list <- read_rds('/lustre1/g/sph_timtsang/yfxu/virus_interference/src_simulation/multipathogons/fitbyoutbreak/traj_match_round1_byvirseas_TOP.rds')
tj_res_list <- read_rds(paste0('/lustre1/g/sph_timtsang/yfxu/virus_interference/src_simulation/multipathogons/self_define_outbreak/Sinusoidal_Forcing/round1/',vir1,"_",vir2,'_traj_match_round1_byvirseas_TOP.rds'))
ci_list <- vector('list', length(seasons))

# tj_res_list<-list(s23=res_r1_fir,`s23-24`=res_r1_sec)
# seasons<- c('s23', 's23-24')
# 
# ci_list <- vector('list', length(seasons))

for (i in 1:length(ci_list)) {
  yr <- seasons[i]

  # if(yr =="s23"){
  #   tj_res_temp <- tj_res_list[[which(!str_detect(names(tj_res_list), "s23-24"))]] %>%
  #     select(all_of(unit_estpars)) %>% na.omit()
  # }else{
  #   tj_res_temp <- tj_res_list[[which(str_detect(names(tj_res_list), yr))]] %>%
  #     select(all_of(unit_estpars)) %>% na.omit()
  # }
  tj_res_temp <- tj_res_list[[which(str_detect(names(tj_res_list), yr))]] %>%
        select(all_of(unit_estpars)) %>% na.omit()

  ci_temp <- as.data.frame(rbind(summarise(tj_res_temp, across(.cols = everything(), min)),
                                 summarise(tj_res_temp, across(.cols = everything(), max))))
  
  if (sens != 'no_rsv_immune') {
    
    # Check that initial conditions can't sum to >1:
    sums <- ci_temp %>% select(I10:R120) %>% rowSums()
    
    if (sums[1] > 1.0) {
      print('Lower bounds sum to more than 1!')
      stop()
    }
    
    if (sums[2] > 1.0) {
      
      # Reduce the upper bounds of R10/R20/R120 proportionally:
      orig_upper_bounds <- ci_temp[2, ] %>% select(R10:R120)
      red_needed <- sums[2] - 0.9999999
      new_upper_bounds <- orig_upper_bounds - (red_needed * (orig_upper_bounds / sum(orig_upper_bounds)))
      
      # Ensure upper bounds still greater than lower:
      orig_lower_bounds <- ci_temp[1, ] %>% select(R10:R120)
      expect_true(all(new_upper_bounds > orig_lower_bounds))
      
      # Ensure upper bounds now sum to 1 or less:
      ci_temp[2, c('R10', 'R20', 'R120')] <- new_upper_bounds
      expect_lt(ci_temp[2, ] %>% select(I10:R120) %>% sum(), 1.0)
      
    }
    
  }
  
  # Store in list:
  ci_list[[i]] <- ci_temp
}
rm(i)

# Get data frame of all ranges:
if (search_type == 'round2_CIs') {
  
  start_range_thetarho <- start_range %>% select(theta_rho1:theta_rho2)
  #/home/yfxu/virus_interference/results/round2_CIs/from_2_1
  file_name<-paste("round2CI_fitforselfdefineoutbreak_SF_",vir1_name,"_",vir2_name,".rds",sep='')
  #write_rds(ci_start, file = paste0(new_dir, file_name))
  start_range <- read_rds(paste0('/home/yfxu/virus_interference/results/round2_CIs/from_2_', which_round - 1, '/',file_name))
  
  if (int_eff == 'both') {
    start_range <- start_range %>%
      bind_cols(start_range_thetarho) %>%
      select(all_of(estpars))
  }
  
  rm(start_range_thetarho)
  
  if (prof_lik) {
    start_range <- start_range[, estpars]
  }
  
} else {
  
  for (i in 1:length(seasons)) {
    
    if (search_type == 'broad') {
      start_range_temp <- unit_start_range
    } else if (search_type == 'round1_CIs') {
      start_range_temp <- ci_list[[i]]
    } else if (search_type == 'round2_CIs') {
      print('ERROR: Round2 CIs used.')
    } else {
      stop('Unrecognized search type!')
    }
    
    names(start_range_temp) <- paste0(seasons[i], '_', names(unit_start_range))
    start_range <- start_range %>%
      bind_cols(start_range_temp)
    rm(start_range_temp)
    
  }
  rm(i)
  
  start_range <- start_range[, estpars]
  
}

# Get starting values for each parameter:
start_values <- sobol_design(lower = setNames(as.numeric(start_range[1, ]), names(start_range[1, ])),
                             upper = setNames(as.numeric(start_range[2, ]), names(start_range[2, ])),
                             nseq = sobol_size)

if (search_type == 'round2_CIs') {
  
  start_values <- start_values %>%
    mutate(phi = if_else(phi > 52.25, phi - 52.25, phi))
  
  if ('phi1' %in% names(start_values)) {
    start_values <- start_values %>%
      mutate(phi1 = if_else(phi1 > 52.25, phi1 - 52.25, phi1),
             phi2 = if_else(phi2 > 52.25, phi2 - 52.25, phi2))
  }
}

# Force no interaction?:
if (sens == 'no_int') {
  
  estpars <- estpars[!(estpars %in% c('theta_lambda1', 'theta_lambda2', 'delta1', 'd2'))]
  true_estpars <- true_estpars[!(true_estpars %in% c('theta_lambda1', 'theta_lambda2', 'delta1', 'd2'))]
  shared_estpars <- shared_estpars[!(shared_estpars %in% c('theta_lambda1', 'theta_lambda2', 'delta1', 'd2'))]
  
  start_values <- start_values[!(names(start_values) %in% c('theta_lambda1', 'theta_lambda2', 'delta1', 'd2'))]
  
}

# Remove eta_ah1, eta_ah2 from consideration?:
if (sens == 'no_ah') {
  
  estpars <- estpars[!(estpars %in% c('eta_ah1', 'eta_ah2'))]
  true_estpars <- true_estpars[!(true_estpars %in% c('eta_ah1', 'eta_ah2'))]
  shared_estpars <- shared_estpars[!(shared_estpars %in% c('eta_ah1', 'eta_ah2'))]
  
  start_values <- start_values[!(names(start_values) %in% c('eta_ah1', 'eta_ah2'))]
  
}

# Check that starting values and estpars are correct:
print(start_range)
print(summary(start_values))
print(estpars)

# Get list of season-specific objective functions:
obj_fun_list <- lapply(po_list, function(ix) {
  create_obj_fxn(ix, estpars = true_estpars)
}) # equivalent to Libbie's GlobalOfun fxn

# Set maximal execution time for each estimation:
if (run_parallel) {
  
  if (sens == 'sinusoidal_forcing') {
    nmins_exec <- time_max * 60 / 2
  } else {
    nmins_exec <- time_max * 60 / 4
  }
  
} else {
  nmins_exec <- time_max * 60 / (sobol_size / no_jobs)
}

print(sprintf("Max estimation time=%.1f min", nmins_exec))

# Get unique identifiers:
sub_start <- (1 + (jobid - 1) * sobol_size / no_jobs) : (jobid * sobol_size / no_jobs)
print(sub_start)

# Fit:
if (run_parallel) {
  
  # Set up parallelization:
  print(detectCores())
  
  n_cores <-length(sub_start)
  
  #use_cluster <- makeCluster(n_cores, outfile = '')
  # use_cluster <- makeSOCKcluster(n_cores, outfile = '')
  # use_cluster <- makeMPIcluster(n_cores, outfile = '')
  # use_cluster <- makePSOCKcluster(n_cores, outfile = '')
  # use_cluster <- makeForkCluster(n_cores, outfile = '')
  # print(use_cluster)
  
  registerDoMC(n_cores)
  # registerDoSNOW(cl = use_cluster)
  #registerDoParallel(cl = use_cluster)
  
  print(getDoParRegistered())
  print(getDoParWorkers())
  
  # Transform start values:
  start_values_tran <- t(
    apply(start_values, 1, function(ix) {
      transform_params(ix, po_list[[1]], seasons, estpars, shared_estpars)
    }, simplify = TRUE)
  )
  x0_trans_names <- colnames(start_values_tran)
  print(x0_trans_names)
  
  # Fit:
  tic <- Sys.time()
  m <- foreach(i = sub_start, .packages = c('tidyverse', 'testthat', 'pomp', 'nloptr')) %dopar% {
    
    x0_trans <- start_values_tran[i, ]
    
    return(
      try(
        nloptr(x0 = x0_trans,
               eval_f = calculate_global_loglik,
               opts = list(algorithm = "NLOPT_LN_SBPLX",
                           maxtime = 60 * nmins_exec,
                           maxeval = -1, # Negative value: criterion is disabled
                           xtol_rel = -1, # Default value: 1e-4
                           print_level = 0))
      )
    )
    
  }
  toc <- Sys.time()
  etime <- toc - tic
  units(etime) <- 'hours'
  print(etime)
  
  # TEMPORARY - SAVE JUST IN CASE PROCESSING CODE RUNS OUT OF TIME:
  saveRDS(m, file = sprintf('/lustre1/g/sph_timtsang/yfxu/virus_interference/src_simulation/multipathogons/self_define_outbreak/Sinusoidal_Forcing/%s_%s_round2/res_TEMP_%s_%s_%d_%d_%d_PARALLEL.rds',
                            vir1,vir2,vir1,
                            int_eff,
                            jobid_orig,
                            sub_start[1],
                            which_round
                            )
  )
  #m<-readRDS('/lustre1/g/sph_timtsang/yfxu/virus_interference/src_simulation/results_round2/h1b/res_TEMP_flu_h1_plus_b_susc_1_1_PARALLEL.rds')
  # Process results:
  m <- lapply(m, function(ix) {
    #ix<-m[[1]]
    if (!inherits(ix, 'try-error')) {
      
      x0_fit <- ix$solution
      names(x0_fit) <- x0_trans_names
      x0_fit_untrans <- back_transform_params(x0_fit, po_list[[1]], seasons, estpars, shared_estpars)
      
      out <- list(estpars = x0_fit_untrans,
                  ll = -ix$objective,
                  conv = ix$status,
                  message = ix$message,
                  niter = ix$iterations)
      # saveRDS(out, file = sprintf('/lustre1/g/sph_timtsang/yfxu/virus_interference/src_simulation/results_round2/h1b/res_out_%s_%s_%d_%d_PARALLEL.rds',
      #                             vir1,
      #                             int_eff,
      #                             jobid_orig,
      #                             which_round
      #                             #sub_start[1]
      #                             )
      #)
    } else {
      out <- 'error'
    }
    
    return(out)
    
  })
  
  # Write to file:
  if (prof_lik) {
    
    saveRDS(m, file = sprintf('/lustre1/g/sph_timtsang/yfxu/virus_interference/src_simulation/multipathogons/self_define_outbreak/Sinusoidal_Forcing/%s_%s_round2/res_%s_%s_%d_%d_%d_%.3f_PARALLEL.rds',
                              vir1,vir2,vir1,
                              int_eff,
                              jobid_orig,
                              sub_start[1],
                              which_round,
                              prof_val)
    )
    
  } else {
    
    saveRDS(m, file = sprintf('/lustre1/g/sph_timtsang/yfxu/virus_interference/src_simulation/multipathogons/self_define_outbreak/Sinusoidal_Forcing/%s_%s_round2/res_%s_%s_%d_%d_%d_PARALLEL.rds',
                              vir1,vir2,vir1,
                              int_eff,
                              jobid_orig,
                              sub_start[1],
                              which_round
                              )
    )
   
  }
  
} else {
  
  for (i in seq_along(sub_start)) {
    
    print(paste0('Estimation: ', sub_start[i]))
    
    # Get start values:
    x0 <- as.numeric(start_values[sub_start[i], ])
    x0_trans <- transform_params(x0, po_list[[1]], seasons, estpars, shared_estpars)
    x0_trans_names <- names(x0_trans)
    
    # Check that parameter transformations correct:
    x0_orig <- back_transform_params(x0_trans, po_list[[1]], seasons, estpars, shared_estpars)
    expect_equal(x0, unname(x0_orig))
    rm(x0_orig)
    
    # Calculate initial log-likelihood:
    print(-1 * calculate_global_loglik(x0_trans))
    
    # Fit models:
    tic <- Sys.time()
    m <- try(
      nloptr(x0 = x0_trans, 
             eval_f = calculate_global_loglik,
             opts = list(algorithm = "NLOPT_LN_SBPLX",
                         maxtime = 60 * nmins_exec,
                         maxeval = -1, # Negative value: criterion is disabled
                         xtol_rel = -1, # Default value: 1e-4
                         print_level = 0))
    )
    toc <- Sys.time()
    etime <- toc - tic
    units(etime) <- 'hours'
    print(etime)
    
    # If estimation is successful, save results:
    if (!inherits(m, 'try-error')) {
      x0_fit <- m$solution
      names(x0_fit) <- x0_trans_names
      x0_fit_untrans <- back_transform_params(x0_fit, po_list[[1]], seasons, estpars, shared_estpars)
      
      out <- list(estpars = x0_fit_untrans,
                  ll = -m$objective,
                  conv = m$status,
                  message = m$message,
                  niter = m$iterations,
                  etime = as.numeric(etime))
      
      # Write to file:
      if (prof_lik) {
        saveRDS(out, file = sprintf('/lustre1/g/sph_timtsang/yfxu/virus_interference/src_simulation/multipathogons/self_define_outbreak/Sinusoidal_Forcing/%s_%s_round2/res_%s_%s_%d_%d_%.3f.rds',
                                    vir1,vir2,vir1,
                                    int_eff,
                                    jobid_orig,
                                    #which_round,
                                    sub_start[i],
                                    prof_val)
        )
      } else {
        saveRDS(out, file = sprintf('/lustre1/g/sph_timtsang/yfxu/virus_interference/src_simulation/multipathogons/self_define_outbreak/Sinusoidal_Forcing/%s_%s_round2/res_%s_%s_%d_%d_%d.rds',
                                    vir1,vir2,vir1,
                                    int_eff,
                                    jobid_orig,
                                    sub_start[i],
                                    which_round)
        )
      }
      
      # Print results:
      print(out$ll)
      print(out$estpars, digits = 2)
      print(out$conv)
      print(out$message)
      
    }
  }
  
  rm(i)
}

# Clean up:
# rm(list = ls())

print('Done!')
save(ci_list,file=sprintf('/lustre1/g/sph_timtsang/yfxu/virus_interference/src_simulation/multipathogons/self_define_outbreak/Sinusoidal_Forcing/%s_%s_round2/ci_list_%s_%s.rds',
                          vir1,vir2,vir1,
                          which_round)
     )
