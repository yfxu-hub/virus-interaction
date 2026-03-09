# ---------------------------------------------------------------------------------------------------------------------
# Run model to assess how timing and coverage of vaccination impact RSV outbreak dynamics, plus sensitivity analyses
# ---------------------------------------------------------------------------------------------------------------------

# Setup

# Load libraries:
library(doMC)
library(tidyverse)
library(pomp)
library(viridis)
library(gridExtra)
library(zoo)
library(readxl)
library(nloptr)
library(parallel)

# Set vaccination coverage levels and time points:
vacc_cov_vec <- round(seq(0.05, 1.0, by = 0.05), digits = 2) # seq(0.1, 1.0, by = 0.1)
# vacc_time_vec <- round(seq(0, 52, by = 1)) # seq(0, 52, by = 2)
setwd("/home/yfxu/virus_interference")

# Set parameters for run:
vir1_name<-'h3'
vir2_name <- 'covid'
which_round<-4
k=1

vir1<-vir1_name
vir2<-vir2_name

Ri_max1 <- 3.0
Ri_max2 <- 4.0
d2_max <- 10.0

debug_bool <- FALSE

# Choose season and vaccine coverage level:
jobid <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID")); print(jobid)
# jobid <- 18

p_vacc <- vacc_cov_vec[(jobid - 1) %% 20 + 1]; print(p_vacc)

# Which assumptions about vaccine efficacy/duration are made?
sens_sim <- as.character(Sys.getenv("SENS")); print(sens_sim)
sens_sim<-'main'

outbreak<-read_excel("/home/yfxu/virus_interference/data/multi_pathogons_outbreak_1.xlsx")%>%as.data.frame()
match_v<-intersect(grep(vir1_name, outbreak$season, value = TRUE),grep(vir2_name, outbreak$season, value = TRUE))
outbreak<-outbreak[which(outbreak$season %in% match_v),]
seasons <- outbreak$season

# Set vaccination efficacy against covid:
# if (sens_sim == 'vacceff60') {
#   vacc_eff <- 0.6
# } else if (sens_sim == 'vacceff95') {
#   vacc_eff <- 0.95
# } else {
#   vacc_eff <- 0.50
# }
vacc_eff <-0.20
print(vacc_eff)

# Set interaction parameter names:
int_params <- c('theta_lambda1', 'theta_lambda2', 'delta1', 'd2')

# Read in MLEs:
# mle_hk <- read_rds('results/MLEs_hk.rds')[1, ]
file_name<-paste("MLE_",vir2_name,"_",vir1_name,"_round",which_round,".rds",sep='')
mle_hk <- read_rds(paste("/home/yfxu/virus_interference/results/",file_name,sep=''))[1, ]
mle_hk_ori<-mle_hk

nm <- colnames(mle_hk)
protected <- grepl(paste0("^",vir2,"-",vir1,"-"), nm)
basenames <- c("rho","theta_lambda","eta_temp","eta_ah")
# c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2','alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
# seasonal<-c('R10','R20','R120','Ri1','Ri2','I10','I20')
for(b in basenames) {
  name1 <- paste0(b, "1")
  name2 <- paste0(b, "2")
  if(name1 %in% nm && name2 %in% nm) {
    if(grepl("^covid-h3-", name1) || grepl("^covid-h3-", name2)) next
    tmp <- mle_hk[name1]
    mle_hk[name1] <- mle_hk[name2]
    mle_hk[name2] <- tmp
  }
}
mle_hk[which(colnames(mle_hk) == "delta1")]<-mle_hk_ori[which(colnames(mle_hk_ori) == "delta1")]*mle_hk_ori[which(colnames(mle_hk_ori) == "d2")]
mle_hk[which(colnames(mle_hk) == "d2")]<-1/mle_hk_ori[which(colnames(mle_hk_ori) == "d2")]

for(yr in seasons){
  vir1_ori_pa<-c('R10','Ri1','I10')
  vir2_ori_pa<-c('R20','Ri2','I20')
  
  mle_hk[which(colnames(mle_hk) %in% paste0(yr,"_",vir1_ori_pa))]<-mle_hk_ori[which(colnames(mle_hk) %in%  paste0(yr,"_",vir2_ori_pa))]
  mle_hk[which(colnames(mle_hk) %in% paste0(yr,"_",vir2_ori_pa))]<-mle_hk_ori[which(colnames(mle_hk) %in%  paste0(yr,"_",vir1_ori_pa))]
}
  
check<-rbind(mle_hk,mle_hk_ori) 
  
# ---------------------------------------------------------------------------------------------------------------------

# Get values for impact of NATURAL influenza infection:
int_param_vals <- mle_hk %>%
  select(all_of(int_params))

if (sens_sim == 'fit_can') {
  int_param_vals <- mle_can %>%
    select(all_of(int_params))
}

# ---------------------------------------------------------------------------------------------------------------------

# Run main simulation study code (Hong Kong / "subtropical" scenario)

# Set parameters for run:
sens <- 'main'
fit_canada <- FALSE

# Get year:
# seasons <- c('s13-14', 's14-15', 's15-16', 's16-17', 's17-18', 's18-19')

for(se in 1:length(seasons)){
  # yr <- seasons[ceiling(jobid / 20)]; print(yr)
  # se<-1
  yr<-seasons[se]
  
  # Specify shared parameters:
  shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                      'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
  
  # Get MLEs for each season:
  mle <- mle_hk
  
  # Perform model checks:
  source('/home/yfxu/virus_interference/code/vaccination_simulation_study/resp_interaction_model_VACC.R')
  
  # Set desired model parameter values:
  model_params <- mle %>%
    dplyr::select(all_of(shared_estpars), contains(yr)) %>%
    rename_with(~str_remove(.x, paste0(yr, '_')), contains(yr)) %>%
    unlist()
  model_params[int_params] <- int_param_vals %>% unlist()
  
  if (sens_sim == 'vacc_can') {
    temp_eff <- mle_can %>% select('theta_lambda1') %>% unlist()
    model_params <- c(model_params, unname(temp_eff), unname(model_params['delta1']), vacc_eff)
    rm(temp_eff)
  } else if (sens_sim == 'deltavacc1month') {
    model_params <- c(model_params, unname(model_params['theta_lambda1']), 7 / 30, vacc_eff)
  } else if (sens_sim == 'deltavacc6months') {
    model_params <- c(model_params, unname(model_params['theta_lambda1']), 7 / 182, vacc_eff)
  } else {
    # model_params <- c(model_params, unname(k*model_params['theta_lambda1']), k*unname(model_params['delta1']), vacc_eff)
    model_params <- c(model_params, k, unname(model_params['delta1']), vacc_eff)
  }
  
  names(model_params)[names(model_params) == ''] <- c('theta_lambda_vacc', 'delta_vacc', 'vacc_eff')
  print(model_params)
  
  resp_mod <- create_SITRxSITR_mod_VACC(dat = dat_pomp,
                                        Ri1_max = Ri_max1,
                                        Ri2_max = Ri_max2,
                                        d2_max = d2_max,
                                        t0_eff = 0,
                                        debug_bool = debug_bool,
                                        sens = sens,
                                        test_diff = FALSE)
  
  resp_mod <- set_model_parameters(resp_mod, model_params, vaccinate = TRUE)
  
  model_params <- parmat(params = coef(resp_mod), nrep = 2)
  
  # Also run where 1) RSV has no impact on flu, and 2) no flu is circulating:
  # model_params['theta_lambda2', ] <- 1.0
  # model_params['I10', ] <- 0
  
  # Set vaccination coverage:
  model_params['theta_lambda_vacc', 1] <- 0
  model_params['delta_vacc', 1] <- 0
  model_params['p_vacc', ] <- c(0, p_vacc)
  model_params['vacc_eff', ] <- c(0, vacc_eff)
  
  # Initiate results data frame:
  res <- NULL
  
  # Loop through vaccination time points:
  vacc_time_vec <- round(seq(0, nrow(dat_pomp), by = 1))
  for (t_vacc in vacc_time_vec) {
    # t_vacc<-5
    # Run deterministic model:
    sim_temp <- run_simulation_with_vaccination(dat_pomp, t_vacc, model_params, Ri_max1, Ri_max2, d2_max, debug_bool, sens, test_diff = FALSE) %>%
      dplyr::select(time, H1:H2, .id) %>%
      mutate(vacc_cov = p_vacc,
             vacc_time = t_vacc,
             season = yr)
    
    res <- res %>% bind_rows(sim_temp)
    
  }
  
  # Write simulation to file:
  write_rds(res, paste0('/home/yfxu/virus_interference/vaccination_simulation_study/simulations/',vir1,'-',vir2,'/sim_determ_', sens_sim, '_','vacc_eff',vacc_eff,'k_',k,'_', yr, '_', p_vacc * 100, 'perc_SUBTROPICAL.rds'))
  
}

# # Check that, if p_vacc = 0 (no vaccination), all vaccine timepoints yield same results:
# res_comp1 <- res %>% filter(.id == 1, vacc_time == min(vacc_time))
# for (t in unique(res$vacc_time)[-which.min(res$vacc_time)]) {
#   res_comp2 <- res %>% filter(.id == 1 & vacc_time == t)
#   print(all.equal(res_comp1$H1, res_comp2$H1))
#   print(all.equal(res_comp1$H2, res_comp2$H2))
# }
# rm(t, res_comp1, res_comp2)

# ---------------------------------------------------------------------------------------------------------------------

# Run main simulation study code (Canada / "temperate" scenario)

# Set parameters for run:
# vir1 <- 'flu'
# sens <- 'sinusoidal_forcing'
# fit_canada <- TRUE

# # Get year:
# seasons <- c('s10-11', 's11-12', 's12-13', 's13-14')
# yr <- seasons[ceiling(jobid / 20)]; print(yr)
# 
# # Specify shared parameters:
# shared_estpars <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
#                     'alpha', 'phi', 'b1', 'b2', 'phi1', 'phi2')

# Check whether all seasons already run:
# if (!is.na(yr)) {
#   
#   # Get MLEs for each season:
#   mle <- mle_can
#   
#   # Perform model checks:
#   source('src/vaccination_simulation_study/resp_interaction_model_VACC.R')
#   
#   # Set desired model parameter values:
#   model_params <- mle %>%
#     dplyr::select(all_of(shared_estpars), contains(yr)) %>%
#     rename_with(~str_remove(.x, paste0(yr, '_')), contains(yr)) %>%
#     unlist()
#   model_params[int_params] <- int_param_vals %>% unlist()
#   
#   if (sens_sim == 'vacc_can') {
#     temp_eff <- mle_can %>% select('theta_lambda1') %>% unlist()
#     model_params <- c(model_params, unname(temp_eff), unname(model_params['delta1']), vacc_eff)
#     rm(temp_eff)
#   } else if (sens_sim == 'deltavacc1month') {
#     model_params <- c(model_params, unname(model_params['theta_lambda1']), 7 / 30, vacc_eff)
#   } else if (sens_sim == 'deltavacc6months') {
#     model_params <- c(model_params, unname(model_params['theta_lambda1']), 7 / 182, vacc_eff)
#   } else {
#     model_params <- c(model_params, unname(model_params['theta_lambda1']), unname(model_params['delta1']), vacc_eff)
#   }
#   
#   names(model_params)[names(model_params) == ''] <- c('theta_lambda_vacc', 'delta_vacc', 'vacc_eff')
#   print(model_params)
#   
#   resp_mod <- create_SITRxSITR_mod_VACC(dat = dat_pomp,
#                                         Ri1_max = Ri_max1,
#                                         Ri2_max = Ri_max2,
#                                         d2_max = d2_max,
#                                         t0_eff = 0,
#                                         debug_bool = debug_bool,
#                                         sens = sens,
#                                         test_diff = TRUE)
#   
#   resp_mod <- set_model_parameters(resp_mod, model_params, vaccinate = TRUE)
#   
#   model_params <- parmat(params = coef(resp_mod), nrep = 2)
#   
#   # Also run where 1) RSV has no impact on flu, and 2) no flu is circulating:
#   # model_params['theta_lambda2', ] <- 1.0
#   # model_params['I10', ] <- 0
#   
#   # Set vaccination coverage:
#   model_params['p_vacc', ] <- c(0, p_vacc)
#   
#   # Initiate results data frame:
#   res <- NULL
#   
#   # Loop through vaccination time points:
#   for (t_vacc in vacc_time_vec) {
#     
#     # Run deterministic model:
#     sim_temp <- run_simulation_with_vaccination(dat_pomp, t_vacc, model_params, Ri_max1, Ri_max2, d2_max, debug_bool, sens, test_diff = TRUE) %>%
#       dplyr::select(time, H1:H2, .id) %>%
#       mutate(vacc_cov = p_vacc,
#              vacc_time = t_vacc,
#              season = yr)
#     
#     res <- res %>% bind_rows(sim_temp)
#     
#   }
#   
#   # Write simulation to file:
#   write_rds(res, paste0('results/vaccination_simulation_study/simulations/sim_determ_', sens_sim, '_', yr, '_', p_vacc * 100, 'perc_TEMPERATE.rds'))
#   
# }

# Print message when completed:
print('Done!')
