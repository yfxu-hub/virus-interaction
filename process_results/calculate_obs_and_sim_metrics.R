# ---------------------------------------------------------------------------------------------------------------------
# Calculate observed attack rates, peak timings, and outbreak concentration for flu and RSV outbreaks, and compare
# to metrics calculated for simulations when parameters are set to best-fit values
# ---------------------------------------------------------------------------------------------------------------------
rm(list = ls())
# Setup

# Load libraries:
library(tidyverse)
library(testthat)

# Load necessary functions:
source('C:/Users/yfxu/SPH Dropbox/Yanfang Xu/Rt_project/virus_interference/code/function/functions_evalutate_res.R')

# ---------------------------------------------------------------------------------------------------------------------

# Calculate observed metrics

# Read in data:
hk_dat <- read_csv("C:/Users/yfxu/SPH Dropbox/Yanfang Xu/Rt_project/virus_interference/code/data/multi_virus_comb_case_ili.csv")#$h1_plus_b_rsv
#can_dat <- read_csv('data/formatted/dat_canada.csv')

# Format:
hk_dat$pop<-rep(7531800,nrow(hk_dat))
hk_dat$n_T<-hk_dat$`No. of specimen tested`
hk_dat<-hk_dat %>%
  mutate(n_T1 = n_T,
         n_T2 = n_T,
         loc = 'hk')
hk_dat$season<-vector(length=nrow(hk_dat))
hk_dat$season[which(hk_dat$year == 2023 & hk_dat$Week <= 45)]<-'s23'
hk_dat$season[-which(hk_dat$year == 2023 & hk_dat$Week <= 45)]<-'s23-24'
hk_dat$n_P2<-hk_dat$`Severe acute respiratory syndrome coronavirus 2`
vir1_name = 'allflu'
if(vir1_name == 'allflu'){
  hk_dat$n_P1<-hk_dat$`Type A`+hk_dat$`Type B`
  hk_dat$flu_h1_plus_b<-hk_dat$`Type A`+hk_dat$`Type B`
}
if(vir1_name == 'h1b'){
  f1$n_P1<-hk_dat$H1+hk_dat$`Type B`
  hk_dat$flu_h1_plus_b<-hk_dat$H1+hk_dat$`Type B`
}
seasons<-c("s23","s23-24")
# hk_dat <- hk_dat %>%
#   select(time, season, pop, n_T:n_P2) %>%
#   mutate(n_T1 = n_T,
#          n_T2 = n_T) %>%
#   select(time:pop, n_T1:n_T2, n_P1:n_P2) %>%
#   mutate(loc = 'hk')

# can_dat <- can_dat %>%
#   select(time, season, pop, n_T1:n_T2, n_P1, n_P2) %>%
#   mutate(loc = 'can')

# Loop through all seasons and calculate metrics (HK):
metrics_hk <- vector('list', length = length(unique(hk_dat$season)))

for (i in 1:length(unique(hk_dat$season))) {
  
  # Get relevant season's data:
  dat_temp <- hk_dat %>% filter(season == unique(hk_dat$season)[i])
  
  # Calculate peak timing, peak intensity, and attack rate:
  obs_metrics <- calculate_metrics(dat_temp)# %>%
  # mutate(pt_diff = pt1 - pt2)
  
  # Calculate peak timing based on proportion (rather than number) positive:
  obs_metrics <- obs_metrics %>% bind_cols(dat_temp %>%
                                             mutate(prop1 = n_P1 / n_T1, prop2 = n_P2 / n_T2) %>%
                                             summarise(pt_prop1 = which.max(prop1), pt_prop2 = which.max(prop2))) %>%
    mutate(pt_diff = pt_prop1 - pt_prop2)
  
  # Calculate attack rate as total proportion positive:
  obs_metrics <- obs_metrics %>% bind_cols(dat_temp %>%
                                             summarise(ar_prop1 = sum(n_P1, na.rm = TRUE) / sum(n_T1, na.rm = TRUE) * 100,
                                                       ar_prop2 = sum(n_P2, na.rm = TRUE) / sum(n_T2, na.rm = TRUE) * 100))
  
  # Calculate duration and concentration (number of weeks containing 75% of reported cases):
  obs_metrics <- obs_metrics %>% bind_cols(calculate_duration_and_concentration(dat_temp))
  
  # Save in list:
  metrics_hk[[i]] <- obs_metrics %>% mutate(season = unique(hk_dat$season)[i])
  
}

metrics_hk <- bind_rows(metrics_hk)
#rm(i, hk_dat, dat_temp, obs_metrics)

# # Loop through all seasons and calculate metrics (Canada):
# metrics_can <- vector('list', length = length(unique(can_dat$season)))
# 
# for (i in 1:length(unique(can_dat$season))) {
#   
#   # Get relevant season's data:
#   dat_temp <- can_dat %>% filter(season == unique(can_dat$season)[i])
#   
#   # Calculate peak timing, peak intensity, and attack rate:
#   obs_metrics <- calculate_metrics(dat_temp)
#   
#   # Calculate peak timing based on proportion (rather than number) positive:
#   obs_metrics <- obs_metrics %>% bind_cols(dat_temp %>%
#                                              mutate(prop1 = n_P1 / n_T1, prop2 = n_P2 / n_T2) %>%
#                                              summarise(pt_prop1 = which.max(prop1), pt_prop2 = which.max(prop2))) %>%
#     mutate(pt_diff = pt_prop1 - pt_prop2)
#   
#   # Calculate attack rate as total proportion positive:
#   obs_metrics <- obs_metrics %>% bind_cols(dat_temp %>%
#                                              summarise(ar_prop1 = sum(n_P1, na.rm = TRUE) / sum(n_T1, na.rm = TRUE) * 100,
#                                                        ar_prop2 = sum(n_P2, na.rm = TRUE) / sum(n_T2, na.rm = TRUE) * 100))
#   
#   # Calculate duration and concentration (number of weeks containing 75% of reported cases):
#   obs_metrics <- obs_metrics %>% bind_cols(calculate_duration_and_concentration(dat_temp))
#   
#   # Save in list:
#   metrics_can[[i]] <- obs_metrics %>% mutate(season = unique(can_dat$season)[i])
#   
# }
# 
# metrics_can <- bind_rows(metrics_can)
# rm(i, can_dat, dat_temp, obs_metrics)

# Print:
metrics_hk %>%
  select(season, contains('1'), pt_diff) %>%
  mutate(pt1 = if_else(season == 's16-17', pt1 + 45 - 53, pt1 + 45 - 52),
         pt_prop1 = if_else(season == 's16-17', pt_prop1 + 45 - 53, pt_prop1 + 45 - 52)) %>%
  pivot_longer(-season) %>%
  group_by(name) %>%
  summarise(mean = mean(value),
            median = median(value),
            min = min(value),
            max = max(value)) %>%
  print()

metrics_hk %>%
  select(season, contains('2')) %>%
  mutate(pt2 = if_else(season == 's16-17', pt2 + 45 - 53, pt2 + 45 - 52),
         pt_prop2 = if_else(season == 's16-17', pt_prop2 + 45 - 53, pt_prop2 + 45 - 52)) %>%
  pivot_longer(-season) %>%
  group_by(name) %>%
  summarise(mean = mean(value),
            median = median(value),
            min = min(value),
            max = max(value)) %>%
  print()

# metrics_can %>%
#   select(season, contains('1'), pt_diff) %>%
#   mutate(pt1 = pt1 + 34 - 52,
#          pt_prop1 = pt_prop1 + 34 - 52) %>%
#   pivot_longer(-season) %>%
#   group_by(name) %>%
#   summarise(mean = mean(value),
#             median = median(value),
#             min = min(value),
#             max = max(value)) %>%
#   print()
# 
# metrics_can %>%
#   select(season, contains('2')) %>%
#   mutate(pt2 = pt2 + 34 - 52,
#          pt_prop2 = pt_prop2 + 34 - 52) %>%
#   pivot_longer(-season) %>%
#   group_by(name) %>%
#   summarise(mean = mean(value),
#             median = median(value),
#             min = min(value),
#             max = max(value)) %>%
#   print()

# Plot "realistic" values:
p1 <- metrics_hk %>% mutate(loc = 'Hong Kong') %>%
  #bind_rows(metrics_can %>% mutate(loc = 'Canada')) %>%
  pivot_longer(-c(season:loc), names_to = 'metric', values_to = 'value') %>%
  mutate(vir = if_else(str_detect(metric, '2'), 'COVID', 'Influenza'),
         metric = str_remove(metric, '1'),
         metric = str_remove(metric, '2')) %>%
  filter(metric %in% c('ar_prop', 'dur', 'pt_prop', 'pt_diff')) %>%
  ggplot(aes(x = loc, y = value)) +
  geom_violin(fill = 'gray90') +
  facet_grid(metric ~ vir, scales = 'free_y') +
  theme_classic() +
  labs(x = 'Virus', y = 'Value')
print(p1)

# ---------------------------------------------------------------------------------------------------------------------

# Calculate simulated metrics and compare to observed

# Get names of fitted parameters:
shared_estpars_hk <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                       'alpha', 'phi', 'eta_temp1', 'eta_temp2', 'eta_ah1', 'eta_ah2')
shared_estpars_can <- c('rho1', 'rho2', 'theta_lambda1', 'theta_lambda2', 'delta1', 'd2',
                        'alpha', 'phi', 'b1', 'b2', 'phi1', 'phi2')
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')

true_estpars_hk <- c(shared_estpars_hk, unit_estpars)
true_estpars_can <- c(shared_estpars_can, unit_estpars)

# Set parameter values necessary for loading models:
prof_lik <- FALSE

# Read in MLEs:
params <- read.csv("C:/Users/yfxu/SPH Dropbox/Yanfang Xu/Rt_project/virus_interference/progress_update/20250806/params_estimate_CI_0806(flu_usingalltypeAandtypeB).csv")%>%as.data.frame()
mle<-params$mean
names(mle)<-params$X
mle<-as.data.frame(t(mle))
#mle_hk <- read_rds('results/MLEs_flu_h1_plus_b.rds')
#mle_can <- read_rds('results/round2_fit/sens/canada/MLEs_flu.rds')
mle_hk<-mle %>% as.data.frame()
# Simulate several outbreaks for each season (HK):
fit_canada <- FALSE
vir1 <- 'flu_h1_plus_b'
true_estpars <- true_estpars_hk
seasons<-c("s23","s23-24")
source('C:/Users/yfxu/SPH Dropbox/Yanfang Xu/Rt_project/virus_interference/code/function/setup_global_likelilhood.R')
#hk_dat <- read_rds('data/formatted/dat_hk_byOutbreak.rds')$'h1_plus_b_rsv'
res_list <- vector('list', length = length(seasons))
for (i in 1:length(seasons)) {
  
  # Set seed:
  set.seed(1078543)
  
  # Run simulations:
  sim_temp <- run_sim(po_list[[i]], seasons[i], mle_hk, shared_estpars_hk, unit_estpars, model_type = 'stochastic', n_sim = 100)
  
  # Calculate metrics (simulated):
  sim_metrics <- calculate_metrics(sim_temp)
  
  # Calculate metrics (observed):
  obs_metrics <- calculate_metrics(hk_dat %>% filter(season == seasons[i]))
  
  # Also calculate correlation coefficients for full outbreaks:
  sim_corr <- sim_temp %>%
    full_join(hk_dat %>% filter(season == seasons[i]),
              by = 'time') %>%
    group_by(.id) %>%
    summarise(corr1 = cor(n_P1.x, n_P1.y, use = 'pairwise.complete.obs'),
              corr2 = cor(n_P2.x, n_P2.y, use = 'pairwise.complete.obs'))
  
  # Determine which simulations within 2wk/25% of observed values:
  sim_metrics <- check_accuracy_metrics(sim_metrics, obs_metrics, seasons[i], pt_acc = 2, pi_acc = 25)
  
  # Store results:
  res_list[[i]] <- sim_metrics
  
}

res_hk <- do.call('rbind', res_list) %>%
  as_tibble() %>%
  mutate(loc = 'hk')


#plot simulate data
seasons<-c("s23","s23-24")
res_list <- vector('list', length = length(seasons))
sim_data<-NULL
for (i in 1:length(seasons)) {
  #i<-1
  # Set seed:
  set.seed(1078543)
  # Run simulations:
  sim_temp <- run_sim(po_list[[i]], seasons[i], mle_hk, shared_estpars_hk, unit_estpars, model_type = 'stochastic', n_sim = 100)
  sim_data<-rbind(sim_data,sim_temp)
}
sim_data<-sim_data %>%
  group_by(season,time)%>%
  summarise(mean_n_P1 = mean(n_P1),
            mean_n_P2 = mean(n_P2))
hk_data1<-hk_dat[,c("Week","start_date","n_P1","n_P2")]
hk_data1<-cbind(hk_data1,sim_data)
hk_data1$start_date<-as.Date(hk_data1$start_date)

df_long_flu <- pivot_longer(hk_data1, cols = c("n_P1", "mean_n_P1"),
                        names_to = "type", values_to = "y")
flu_p<-ggplot() +
  geom_point(data = subset(df_long_flu, type == "n_P1"),
             aes(x = start_date, y = y, color = type), size = 1.5) +
  geom_line(data = subset(df_long_flu,  type == "mean_n_P1"),
            aes(x = start_date, y = y, color = type), size = 1.2) +
  scale_color_manual(values = c("n_P1" = "black", "mean_n_P1" = "blue"),
                     labels = c("n_P1" = "real data", "mean_n_P1" = "simulation")) +
  #scale_shape_manual(values = c("n_P1" = 16, "mean_n_P1" = NA)) +
  #scale_linetype_manual(values = c("n_P1" = "blank", "mean_n_P1" = "solid")) +
  labs(title="Influenza",
       x="Date",
       y="Case number")+
  theme_minimal()


df_long_covid <- pivot_longer(hk_data1, cols = c("n_P2", "mean_n_P2"),
                            names_to = "type", values_to = "y")
covid_p<-ggplot() +
  geom_point(data = subset(df_long_covid, type == "n_P2"),
             aes(x = start_date, y = y, color = type), size = 1.5) +
  geom_line(data = subset(df_long_covid,  type == "mean_n_P2"),
            aes(x = start_date, y = y, color = type), size = 1.2) +
  scale_color_manual(values = c("n_P2" = "black", "mean_n_P2" = "blue"),
                     labels = c("n_P2" = "real data", "mean_n_P2" = "simulation")) +
  #scale_shape_manual(values = c("n_P1" = 16, "mean_n_P1" = NA)) +
  #scale_linetype_manual(values = c("n_P1" = "blank", "mean_n_P1" = "solid")) +
  labs(title="COVID",
       x="Date",
       y="Case number")+
  theme_minimal()
covid_p

# # Simulate several outbreaks for each season (Canada):
# fit_canada <- TRUE
# vir1 <- 'flu'
# true_estpars <- true_estpars_can
# 
# source('src/functions/setup_global_likelilhood.R')
# #can_dat <- read_csv('data/formatted/dat_canada.csv')
# 
# res_list <- vector('list', length = length(seasons))
# for (i in 1:length(seasons)) {
#   
#   # Set seed:
#   set.seed(1078543)
#   
#   # Run simulations:
#   sim_temp <- run_sim(po_list[[i]], seasons[i], mle_can, shared_estpars_can, unit_estpars, model_type = 'stochastic', n_sim = 100)
#   
#   p_temp <- ggplot() +
#     geom_line(data = sim_temp, aes(x = time, y = n_P2, group = .id), col = 'gray80', alpha = 0.5) +
#     geom_point(data = can_dat %>% filter(season == seasons[i]), aes(x = time, y = n_P2)) +
#     theme_classic()
#   print(p_temp)
#   
#   # Calculate metrics (simulated):
#   sim_metrics <- calculate_metrics(sim_temp)
#   
#   # Calculate metrics (observed):
#   obs_metrics <- calculate_metrics(can_dat %>% filter(season == seasons[i]))
#   
#   # Also calculate correlation coefficients for full outbreaks:
#   sim_corr <- sim_temp %>%
#     full_join(can_dat %>% filter(season == seasons[i]),
#               by = 'time') %>%
#     group_by(.id) %>%
#     summarise(corr1 = cor(n_P1.x, n_P1.y, use = 'pairwise.complete.obs'),
#               corr2 = cor(n_P2.x, n_P2.y, use = 'pairwise.complete.obs'))
#   
#   # Determine which simulations within 2wk/25% of observed values:
#   sim_metrics <- check_accuracy_metrics(sim_metrics, obs_metrics, seasons[i], pt_acc = 2, pi_acc = 25)
#   
#   # Store results:
#   res_list[[i]] <- sim_metrics
#   
# }
# 
# res_can <- do.call('rbind', res_list) %>%
#   as_tibble() %>%
#   mutate(loc = 'canada')

# # Combine results from all locations:
# res <- bind_rows(res_hk, res_can)
# print(res)

# Clean up
#rm(list = ls())
