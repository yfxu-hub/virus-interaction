# ---------------------------------------------------------------------------------------------------------------------
# Code to fit DETERMINISTIC flu/RSV interaction model
# Round 1: Fit each season individually (w/ or w/o interaction) to obtain good start values for Ri/initial conditions
# ---------------------------------------------------------------------------------------------------------------------
rm(list = ls())
library(tidyverse)
library(pomp)
library(viridis)
library(gridExtra)
library(zoo)
library(readxl)
# Setup

# Set seed:
set.seed(749501349)

# Get cluster environmental variables:
jobid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")); print(jobid)
no_jobs <- as.integer(Sys.getenv("NOJOBS")); print(no_jobs)
sobol_size <- as.integer(Sys.getenv("SOBOLSIZE")); print(sobol_size)
search_type <- as.character(Sys.getenv("SEARCHTYPE")); print(search_type)
sens <- as.character(Sys.getenv("SENS")); print(sens)
fit_canada <- FALSE


vir1_name<-'covid'
vir2_name<-'rsv'

if(vir2_name == 'rsv'){
  rho1_ori <- 1.5
}else{
  rho1_ori <- 0.5
}
# Set parameters for run:
#outbreak<-read_excel("/home/yfxu/virus_interference/data/multi_pathogons_outbreak.xlsx")%>%as.data.frame()

# outbreak<-read_excel("/home/yfxu/virus_interference/data/multi_pathogons_outbreak_allperiod.xlsx")%>%as.data.frame()
outbreak<-read_excel("/home/yfxu/virus_interference/data/multi_pathogons_outbreak_1.xlsx")%>%as.data.frame()
match_v<-intersect(grep(vir1_name, outbreak$season, value = TRUE),grep(vir2_name, outbreak$season, value = TRUE))
outbreak<-outbreak[which(outbreak$season %in% match_v),]

season<-outbreak$season

#season<-c('s22-23','s23-24','s24-25')

jobid <- (jobid - 1) %% no_jobs + 1; print(jobid)

vir1<-vir1_name
vir2 <- vir2_name

time_max <- 20 # 5.0 # Maximal execution time (in hours)
debug_bool <- FALSE

Ri_max1 <- 4.0
Ri_max2 <- 3.0
d2_max <- 10.0

# # Fit for synthetic data from age-structured model?:
# age_structured <- TRUE

# ---------------------------------------------------------------------------------------------------------------------

# Run trajectory matching
# Note: Separate runs for each season

# Load data and create pomp object:
for(yr in season){
  #yr<-season[2]
source('/home/yfxu/virus_interference/code/resp_interaction_model.R')

  # Check that sufficient epidemic activity:
  if (exists('resp_mod')) {
    
    # Set start ranges for estimated parameters:
    if (sens == 'no_rsv_immune') {
      estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'rho1', 'rho2')
    } else {
      estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120', 'rho1', 'rho2')
    }
    
    start_range <- data.frame(Ri1 = c(1.0, Ri_max1),
                              Ri2 = c(1.0, Ri_max2),
                              I10 = c(0.000001, 0.5),
                              I20 = c(0.000001, 0.5),
                              R10 = c(0, 0.5),
                              R20 = c(0, 0.5),
                              R120 = c(0, 0.5),
                              rho1 = c(0.001, rho1_ori),
                              rho2 = c(0.001, 0.5),
                              theta_lambda1 = c(0, 1.0),
                              theta_lambda2 = c(0, 1.0),
                              theta_rho1 = c(0, 1.0),
                              theta_rho2 = c(0, 1.0),
                              delta1 = c(7 / 60, 7),
                              d2 = c(7 / 60, 7),
                              alpha = c(0, 0.5),
                              phi = c(0, 52.25),
                              eta_temp1 = c(-0.5, 0.5),
                              eta_temp2 = c(-0.5, 0.5),
                              eta_ah1 = c(-0.5, 0.5),
                              eta_ah2 = c(-0.5, 0.5))
    start_range <- start_range[, estpars]
    
    if (search_type == 'broad') {
      start_values <- sobol_design(lower = setNames(as.numeric(start_range[1, ]), names(start_range[1, ])),
                                   upper = setNames(as.numeric(start_range[2, ]), names(start_range[2, ])),
                                   nseq = sobol_size)
    }
    
    print(estpars)
    print(start_range)
    print(summary(start_values))
    
    # Create objective function for call to nloptr:
    obj_fun <- traj_objfun(data = resp_mod,
                           est = estpars,
                           partrans = resp_mod@partrans,
                           verbose = TRUE)
    
    # Set maximal execution time for each estimation:
    nmins_exec <- time_max * 60 / (sobol_size / no_jobs)
    print(sprintf("Max estimation time=%.1f min", nmins_exec))
    
    # Get unique identifiers:
    sub_start <- (1 + (jobid - 1) * sobol_size / no_jobs) : (jobid * sobol_size / no_jobs)
    
    # Loop through start values and perform trajectory matching:
    for (i in seq_along(sub_start)) {
      
      print(paste0('Estimation: ', sub_start[i]))
      
      # Get param start values:
      x0 <- as.numeric(start_values[sub_start[i], ])
      coef(resp_mod, estpars) <- x0
      x0_trans <- coef(resp_mod, estpars, transform = TRUE)
      print(-1 * obj_fun(x0_trans))
      
      # Run trajectory matching using subplex algorithm:
      # http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms
      tic <- Sys.time()
      m <- try(
        nloptr(x0 = unname(x0_trans),
               eval_f = obj_fun,
               opts = list(algorithm = 'NLOPT_LN_SBPLX',
                           maxtime = 60.0 * nmins_exec,
                           maxeval = -1, # disabled
                           xtol_rel = -1, # disabled; default: 1e-4
                           print_level = 0))
      )
      toc <- Sys.time()
      etime <- toc - tic
      units(etime) <- 'mins'
      print(etime)
      
      # If estimation is successful, save results:
      if (!inherits(m, 'try-error')) {
        coef(resp_mod, estpars, transform = TRUE) <- m$solution
        
        # Collect all results:
        out <- list(allpars = coef(resp_mod),
                    estpars = coef(resp_mod, estpars),
                    ll = -m$objective,
                    conv = m$status,
                    message = m$message,
                    niter = m$iterations,
                    etime = as.numeric(etime))
        
        # Write to file:
        saveRDS(out,
                file = sprintf('/lustre1/g/sph_timtsang/yfxu/virus_interference/src_simulation/multipathogons/self_define_outbreak/Sinusoidal_Forcing/round1/res_%s_%s_%s_%d_%d.rds',
                               vir1, vir2,
                               as.character(yr),
                               jobid,
                               sub_start[i])
        )
        
        # Print results:
        print(out$ll)
        print(out$estpars, digits = 2)
        print(out$conv)
        print(out$message)
      }
      
    }
    
  } else {
    print('No data for given season.')
  }
  
}

print('Done.')
# # Process "round 1" results
# files <- list.files(
#   path       = "/lustre1/g/sph_timtsang/yfxu/virus_interference/src_simulation/results_round1_fit_by_array/fluab", 
#   pattern    = "\\.rds$", 
#   full.names = TRUE
# )
# files_fir<-files[!grepl("s23-24", files)]
# files_sec<-files[grepl("s23-24", files)]
# 
# res_list1 <- lapply(files_fir, readRDS)
# names(res_list1) <- basename(files_fir)
# save(res_list1,file="/lustre1/g/sph_timtsang/yfxu/virus_interference/src_simulation/results_round1_fit_by_array/res_list_firstseason.RData")
# 
# res_list2 <- lapply(files_sec, readRDS)
# names(res_list2) <- basename(files_sec)
# save(res_list2,file="/lustre1/g/sph_timtsang/yfxu/virus_interference/src_simulation/results_round1_fit_by_array/res_list_secondseason.RData")
# 
# res_list <- lapply(files_fir, readRDS)
# names(res_list) <- basename(files)
# # Compile:
# res_r1 <- lapply(res_list, getElement, 'estpars') %>%
#   bind_rows() %>%
#   bind_cols('loglik' = lapply(res_list, getElement, 'll') %>%
#               unlist()) %>%
#   bind_cols('message' = lapply(res_list, getElement, 'message') %>%
#               unlist()) %>%
#   mutate(year = 's23-24') %>%
#   select(year, Ri1:message)
# #rm(res_list)
# 
# # Remove where no convergence occurs:
# res_r1 <- res_r1 %>%
#   filter(!str_detect(message, 'maxtime')) %>%
#   select(-message)
# 
# # Check that fit parameters do not produce simulations where state variables go negative:
# p_mat <- parmat(coef(resp_mod), nrep = nrow(res_r1))
# 
# for (param in estpars) {
#   p_mat[param, ] <- res_r1 %>% pull(param)
# }
# 
# expect_equal(trajectory(resp_mod, params = p_mat, format = 'data.frame') %>%
#                pivot_longer(X_SS:H2, names_to = 'state') %>%
#                filter(value < 0) %>%
#                nrow(), 0)
# 
# #rm(p_mat, param)
# 
# # Arrange results by log-likelihood and keep only best estimates:
# res_r1 <- res_r1 %>%
#   arrange(desc(loglik))
# 
# no_best <- nrow(subset(res_r1, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = length(estpars))))
# res_r1 <- res_r1[1:no_best, ]
# save(res_r1,file="/lustre1/g/sph_timtsang/yfxu/virus_interference/src_simulation/results_round1_fit_by_array/firstseasonres_r1_inter499.RData")
# 
# # here all 10 fall within an acceptable range of the max log likelihood; in the real analysis,
# # this is done with fits from all 500 starting parameter sets
# 
