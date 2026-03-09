# ---------------------------------------------------------------------------------------------------------------------
# Set up model of flu/RSV transmission
# Note: Adapted from code of Dr. Matthieu Domenech de Celles
# ---------------------------------------------------------------------------------------------------------------------

# Load packages
library(tidyverse)
library(pomp)
library(viridis)
library(gridExtra)

# Load required functions:
source('/home/yfxu/virus_interference/code/functions_flu_RSV.R')
source('/home/yfxu/virus_interference/code/test_code.R')
source('/home/yfxu/virus_interference/code/vaccination_simulation_study/functions_simulation_study.R')

# Set test_diff_use
if (fit_canada) {
  test_diff_use <- TRUE
} else {
  test_diff_use <- FALSE
}

# # Set additional variables locally (for testing):
# vir1 <- 'flu_h1_plus_b'
# vir2 <- 'rsv'
# yr <- 's13-14'
# # yr <- 's15-16'
# # yr <- 's16-17'
# # yr <- 's17-18'
# # yr <- 's18-19'
# Ri_max1 <- 2.0
# Ri_max2 <- 3.0
# d2_max <- 10.0
# debug_bool <- TRUE

# ---------------------------------------------------------------------------------------------------------------------

# Load and format data

# Read in data and choose season:
# if (fit_canada) {
#   
#   can_dat <- read_csv('data/formatted/dat_canada.csv')
#   dat_pomp <- can_dat %>%
#     filter(season == yr)
#   
# } else {
#   
#   hk_dat <- read_rds('data/formatted/dat_hk_byOutbreak.rds')
#   dat_pomp <- hk_dat[[paste(str_sub(vir1, 5, str_length(vir1)), vir2, sep = '_')]] %>%
#     filter(season == yr)
#   
# }
df1 <- read_csv("/home/yfxu/virus_interference/data/multi_virus_comb_case_ili.csv")%>%as.data.frame()

colnames(df1)[which(colnames(df1) == "Severe acute respiratory syndrome coronavirus 2")]<-'covid'
colnames(df1)[which(colnames(df1) == "H1")]<-'h1'
colnames(df1)[which(colnames(df1) == "H3")]<-'h3'
colnames(df1)[which(colnames(df1) == "Type B")]<-'flub'
colnames(df1)[which(colnames(df1) == "RSV")]<-'rsv'

df1<-df1[5:nrow(df1),]

# df1<-df1[5:NROW(df1),]
df1[] <- lapply(df1, function(x) {
  if (is.numeric(x)) na.approx(x, na.rm = FALSE) else x
})

#impute the missing data
df1$n_T<-df1$`No. of specimen tested`
df1$n_P1<-df1[,vir1_name]
df1$n_P2<-df1[,vir2_name]
df1$i_ILI<-df1$PMP/1000
df1$pop<-rep(7531800,nrow(df1))
df1$temp<-as.vector(scale(df1$avg_temperature))
df1$ah<-as.vector(scale(df1$avg_abs_humidity))

df1$rhino_inc<-rep(0,nrow(df1))
df1$h3_inc<-rep(0,nrow(df1))

for(vir in c(vir1_name,vir2_name)){
  if(vir == "covid"){
    gamma_v <- 14
  }
  if(vir == "h1"|vir == "h3"|vir == "flub"){
    gamma_v <- 5
  }
  if(vir == 'rsv'){
    gamma_v <- 10
  }
  if(vir == vir1_name){
    gamma_vir1<-gamma_v
  }
  if(vir == vir2_name){
    gamma_vir2<-gamma_v
  }
}

match_v<-intersect(grep(vir1_name, outbreak$season, value = TRUE),grep(vir2_name, outbreak$season, value = TRUE))
outbreak<-outbreak[which(outbreak$season %in% match_v),]

df1$season<-rep(NA,length=nrow(df1))
for(r in 1:nrow(outbreak)){
  df1$season[which(df1$start_date >= as.Date(outbreak[r,"start_date"]) & df1$start_date <= as.Date(outbreak[r,"end_date"]))]<-outbreak[r,1]
}
df1<-df1 %>% group_by(season) %>% mutate(time = row_number())

hk_dat<-df1

dat_pomp<-hk_dat %>%
  filter(season == yr)
nrow_check <- nrow(dat_pomp)
nrow_check
# # Get climate data:
# if (fit_canada) {
#   
#   dat_pomp <- dat_pomp %>%
#     mutate(temp = 0,
#            ah = 0)
#   expect_true(nrow(dat_pomp) == nrow_check)
#   
# } else {
#   
#   dat_clim <- read_csv('data/formatted/clim_dat_hk_NORM.csv')
#   
#   dat_pomp <- dat_pomp %>%
#     inner_join(dat_clim,
#                by = c('Year' = 'year',
#                       'Week' = 'week')) %>%
#     select(time:pop, temp, ah, rh)
#   expect_true(nrow(dat_pomp) == nrow_check)
#   
#   rm(dat_clim)
#   
# }

# If no data for this season, skip:
if (nrow(dat_pomp) > 0) {
  
  # Format data:
  # if (!fit_canada) {
  #   dat_pomp <- dat_pomp %>%
  #     rename('i_ILI' = 'PMP') %>% #or GOPC
  #     mutate(i_ILI = i_ILI / 1000)
  #   # https://www.censtatd.gov.hk/en/
  #   # https://www.censtatd.gov.hk/en/web_table.html?id=1A#
  # }
  
  # Plot data:
  if (debug_bool) {
    
    # Plot ILI incidence:
    p1 <- ggplot(data = dat_pomp, aes(x = time, y = i_ILI)) + geom_line() +
      labs(x = 'Time (Weeks)', y = 'ILI Incidence Rate (per 1000 Consultations)') +
      theme_classic()
    # print(p1)
    
    if (fit_canada) {
      p2 <- ggplot(data = dat_pomp %>%
                     pivot_longer(n_P1:n_P2, names_to = 'virus', values_to = 'n_pos') %>%
                     pivot_longer(n_T1:n_T2, names_to = 'virus1', values_to = 'n_T') %>%
                     mutate(virus = str_remove(virus, 'n_P'),
                            virus1 = str_remove(virus1, 'n_T')) %>%
                     filter(virus == virus1),
                   aes(x = time, y = n_pos / n_T, color = virus)) +
        geom_line() + labs(x = 'Time (Weeks)', y = 'Positivity Fraction') +
        theme_classic()
    } else {
      p2 <- ggplot(data = dat_pomp %>%
                     pivot_longer(n_P1:n_P2, names_to = 'virus', values_to = 'n_pos'),
                   aes(x = time, y = n_pos / n_T, color = virus)) +
        geom_line() + labs(x = 'Time (Weeks)', y = 'Positivity Fraction') +
        theme_classic()
    }
    print(p2)
    
  }
  
  # ---------------------------------------------------------------------------------------------------------------------
  
  # Create pomp model and run basic model checks
  
  # Create model:
  if (fit_canada) {
    
    resp_mod <- create_SITRxSITR_mod_VACC(dat = dat_pomp,
                                          Ri1_max = Ri_max1,
                                          Ri2_max = Ri_max2,
                                          d2_max = d2_max,
                                          t0_eff = 0,
                                          debug_bool = debug_bool,
                                          sens = sens,
                                          test_diff = TRUE)
    
  } else {
    
    resp_mod <- create_SITRxSITR_mod_VACC(dat = dat_pomp,
                                          Ri1_max = Ri_max1,
                                          Ri2_max = Ri_max2,
                                          d2_max = d2_max,
                                          t0_eff = 0,
                                          debug_bool = debug_bool,
                                          sens = sens,
                                          test_diff = FALSE)
    
  }
  
  # If parameter information available, update model parameters:
  if (exists('mle')) {
    
    model_params <- mle %>%
      dplyr::select(all_of(shared_estpars), contains(yr)) %>%
      rename_with(~str_remove(.x, paste0(yr, '_')), contains(yr)) %>%
      unlist()
    
    resp_mod <- set_model_parameters(resp_mod, model_params, vaccinate = FALSE)
    
  }
  
  # Check transformations:
  check_transformations(resp_mod)
  
  # Check parameters:
  check_params(resp_mod)
  
  # Check initial conditions:
  expect_true(all.equal(sum(rinit(resp_mod)), as.numeric(coef(resp_mod, 'N'))))
  
  # Check constant population size:
  check_correct_N_CONST(resp_mod, unique(dat_pomp$pop), n_sim = 0)
  
  # Run deterministic simulation:
  sim_determ <- trajectory(object = resp_mod, format = 'data.frame') %>%
    dplyr::select(H1:.id) %>%
    pivot_longer(H1:H2, names_to = 'Vir', values_to = 'Inc')
  p3 <- ggplot(data = sim_determ, aes(x = time, y = Inc, group = Vir, col = Vir)) +
    geom_line() + geom_point() +
    labs(x = 'Time (Weeks)', y = 'Incidence', col = 'Virus') +
    theme_classic()
  if (debug_bool) print(p3)
  
  # Check that measurement density model works:
  ll <- logLik(traj_objfun(data = resp_mod))
  if (debug_bool) print(ll)
  
  # Check that dynamics are independent when there is no interaction:
  if (!exists('mle')) {
    p4 <- check_independent_dynamics(resp_mod)
    if (debug_bool) print(p4)
  }
  
  # Check that all vaccination compartments are always 0:
  sim_determ <- trajectory(object = resp_mod, format = 'data.frame') %>%
    dplyr::select(V_SS1:V_RS)
  expect_true(all(sim_determ == 0))
  
  # ---------------------------------------------------------------------------------------------------------------------
  
  # Add vaccination and continue checks
  
  # First, vaccinate at the beginning of the season:
  t_vacc <-0
  
  if (exists('model_params')) {
    
    model_params <- c(model_params, 0.1, 0.8, 0.2)
    names(model_params)[names(model_params) == ''] <- c('theta_lambda_vacc', 'vacc_eff', 'p_vacc')
    
  } else {
    
    model_params <- c(0.1, 0.8, 0.2)
    names(model_params) <- c('theta_lambda_vacc', 'vacc_eff', 'p_vacc')
    
  }
  
  sim_determ <- run_simulation_with_vaccination(dat_pomp, t_vacc, model_params, Ri_max1, Ri_max2, d2_max, debug_bool, sens, test_diff = test_diff_use)
  
  # Check that correct proportion of people are vaccinated:
  if (debug_bool) {
    print(paste0('VSS1 should be: ', round(unique(dat_pomp$pop) * model_params['p_vacc'] * (1.0 - sum(coef(resp_mod, c('I10', 'I20', 'R10', 'R20', 'R120')))))))
    print(paste0('VSR should be: ', round(unique(dat_pomp$pop) * model_params['p_vacc'] * coef(resp_mod, 'R20'))))
    print(paste0('VRS should be: ', round(unique(dat_pomp$pop) * model_params['p_vacc'] * coef(resp_mod, 'R10'))))
    print(paste0('VSI should be: ', round(unique(dat_pomp$pop) * model_params['p_vacc'] * coef(resp_mod, 'I20'))))
  }
  
  # Check that population size is constant:
  check_correct_N_CONST_VACC(sim_determ, unique(dat_pomp$pop))
  
  # Check that no compartments go negative:
  expect_true(all(sim_determ %>% dplyr::select(X_SS:H2) >= 0))
  
  # # Check that error is returned if vaccination rate is too high:
  # This will be relevant if we decide to vaccinate as a proportion of the population, rather than as a proportion of X_SS+X_SR+X_RS
  
  # Now vaccinate in the middle of the season, and recheck:
  t_vacc = round(nrow(dat_pomp)/2)
  sim_determ <- run_simulation_with_vaccination(dat_pomp, t_vacc, model_params, Ri_max1, Ri_max2, d2_max, debug_bool, sens, test_diff = test_diff_use)
  
  # Check that population size is constant:
  check_correct_N_CONST_VACC(sim_determ, unique(dat_pomp$pop))
  
  # Check that no compartments go negative:
  expect_true(all(sim_determ %>% dplyr::select(X_SS:H2) >= 0))
  # print(all(sim_determ %>% dplyr::select(X_SS:H2) >= 0))
  
  # Check for t_vacc == 1:
  t_vacc = 1
  sim_determ <- run_simulation_with_vaccination(dat_pomp, t_vacc, model_params, Ri_max1, Ri_max2, d2_max, debug_bool, sens, test_diff = test_diff_use)
  
  # Check that population size is constant:
  check_correct_N_CONST_VACC(sim_determ, unique(dat_pomp$pop))
  
  # Check that no compartments go negative:
  expect_true(all(sim_determ %>% dplyr::select(X_SS:H2) >= 0))
  # print(all(sim_determ %>% dplyr::select(X_SS:H2) >= 0))
  
  # Check that dynamics remain independent, so long as vaccine has no effect on either virus:
  if (!exists('mle')) {
    
    t_vacc <- 0
    model_params <- parmat(params = coef(resp_mod), nrep = 3)
    model_params['p_vacc', ] <- 0.2
    model_params['I10', ] <- c(1e-5, 1e-5, 0)
    model_params['I20', ] <- c(1e-5, 0, 1e-5)
    p5 <- check_independent_dynamics_VACC(dat_pomp, t_vacc, model_params, Ri_max1, Ri_max2, d2_max, debug_bool)
    if (debug_bool) print(p5)
    
    t_vacc <- 10
    p6 <- check_independent_dynamics_VACC(dat_pomp, t_vacc, model_params, Ri_max1, Ri_max2, d2_max, debug_bool)
    if (debug_bool) print(p6)
    
  }
  
  # # Vaccine affects flu, but not RSV:
  # model_params['vacc_eff', ] <- 0.8
  # model_params['theta_lambda_vacc', ] <- 1.0
  # check_independent_dynamics_VACC(dat_pomp, t_vacc, model_params, Ri_max1, Ri_max2, d2_max, debug_bool)
  # t_vacc <- 0
  # check_independent_dynamics_VACC(dat_pomp, t_vacc, model_params, Ri_max1, Ri_max2, d2_max, debug_bool)
  # 
  # # Vaccine affects RSV, but not flu:
  # model_params['vacc_eff', ] <- 0
  # model_params['theta_lambda_vacc', ] <- 0.1
  # check_independent_dynamics_VACC(dat_pomp, t_vacc, model_params, Ri_max1, Ri_max2, d2_max, debug_bool)
  # 
  # # Vaccine affects both:
  # model_params['vacc_eff', ] <- 0.8
  # model_params['theta_lambda_vacc', ] <- 0.1
  # check_independent_dynamics_VACC(dat_pomp, t_vacc, model_params, Ri_max1, Ri_max2, d2_max, debug_bool)
  # 
  # model_params['delta_vacc', ] <- 7/60
  # check_independent_dynamics_VACC(dat_pomp, t_vacc, model_params, Ri_max1, Ri_max2, d2_max, debug_bool)
  # 
  # model_params['delta_vacc', ] <- 7 / 5
  # model_params['theta_lambda1', ] <- 0.1
  # check_independent_dynamics_VACC(dat_pomp, t_vacc, model_params, Ri_max1, Ri_max2, d2_max, debug_bool)
  # 
  # model_params['delta_vacc', ] <- 7 / 60
  # model_params['delta1', ] <- 7 / 60
  # check_independent_dynamics_VACC(dat_pomp, t_vacc, model_params, Ri_max1, Ri_max2, d2_max, debug_bool)
  
  # Check that only flu impacted if vaccine only affects flu:
  t_vacc <- 0
  model_params <- parmat(params = coef(resp_mod), nrep = 2)
  model_params['vacc_eff', ] <- 0.8
  model_params['p_vacc', ] <- c(0.0, 0.2)
  model_params['theta_lambda1', ] <- 1.0
  model_params['theta_lambda2', ] <- 1.0
  
  p7 <- check_single_virus_impact(dat_pomp, t_vacc, model_params, Ri_max1, Ri_max2, d2_max, debug_bool, sens, test_diff = test_diff_use)
  if (debug_bool) print(p7)
  
  t_vacc <-5
  p8 <- check_single_virus_impact(dat_pomp, t_vacc, model_params, Ri_max1, Ri_max2, d2_max, debug_bool, sens, test_diff = test_diff_use)
  if (debug_bool) print(p8)
  
  # Check that only RSV impacted if vaccine only affects RSV:
  t_vacc <- 0
  model_params <- parmat(params = coef(resp_mod), nrep = 2)
  model_params['theta_lambda_vacc', ] <- 0.1
  model_params['p_vacc', ] <- c(0.0, 0.2)
  model_params['theta_lambda1', ] <- 1.0
  model_params['theta_lambda2', ] <- 1.0
  
  p9 <- check_single_virus_impact(dat_pomp, t_vacc, model_params, Ri_max1, Ri_max2, d2_max, debug_bool, sens, test_diff = test_diff_use)
  if (debug_bool) print(p9)
  
  t_vacc <- 10
  p10 <- check_single_virus_impact(dat_pomp, t_vacc, model_params, Ri_max1, Ri_max2, d2_max, debug_bool, sens, test_diff = test_diff_use)
  if (debug_bool) print(p10)
  
  # Clean up:
  rm(sim_determ, p3, p7, p8, p9, p10, ll, resp_mod, t_vacc, model_params)
  
  if (!exists('mle')) {
    rm(p4, p5, p6)
  }
  
}

# ---------------------------------------------------------------------------------------------------------------------
