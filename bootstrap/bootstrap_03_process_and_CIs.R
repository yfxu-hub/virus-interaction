# ---------------------------------------------------------------------------------------------------------------------
# Process results of parametric bootstrapping and output 95% CIs and plots
# ---------------------------------------------------------------------------------------------------------------------
rm(list = ls())
# Load libraries:
library(tidyverse)
library(testthat)
library(rethinking)
library(gt)
library(stringr)
library(readxl)
library(dplyr)

vir1<-'h1'
vir2<-'rsv'
which_round<-4
# Results from Canada?:
fit_canada <- FALSE
#=====================load boostrapt===========================
# Get names of all results files:
#file_list <- list.files(path = 'results/bootstrapping/flu_H1_plus_B/', full.names = TRUE)
file_list <- list.files(path = 'C:/Users/yfxu/SPH Dropbox/Yanfang Xu/Rt_project/virus_interference/code/data/boostrapt/self_define_outbreak/round4/',
                        pattern=paste0('^',"res_",vir1,"_",vir2),full.names = TRUE)
# file_list <- file_list[!grepl("res_out_flu_h1_plus_b_susc", basename(file_list))]
# file_list<-"C:/Users/yfxu/SPH Dropbox/Yanfang Xu/Rt_project/virus_interference/code/data/results_flu_A_B/result_bootstrap/round3/res_flu_h1_plus_b_susc_1_3_PARALLEL_10times.rds"

# Are results from a sensitivity analysis?:
  sens <- 'main'


# Ensure no results missing:
# if (str_detect(file_list[[1]], 'PARALLEL')) {
#   
#   expect_true(length(file_list) == 500)
#   
#   run_list <- str_split(file_list, '_') %>%
#     purrr::map(~ .x[!is.na(as.numeric(.x))]) %>%
#     unlist() %>%
#     as.numeric()
#   print(which(!(1:500 %in% run_list)))
#   
# } else {
#   
#   expect_true(length(file_list) == 500 * 10)
#   
# }

# Read in all results:
res_full = list()
for (i in seq_along(file_list)) {
  res_full[[i]] <- read_rds(file_list[[i]])
}
rm(i)

# Get parameter estimates and log-likelihoods:
if (str_detect(file_list[[1]], 'PARALLEL')) {
  
  res_df <- lapply(res_full, function(ix) {
    lapply(ix, getElement, 'estpars') %>%
      bind_rows() %>%
      bind_cols('loglik' = lapply(ix, getElement, 'll') %>%
                  unlist()) %>%
      # mutate(dataset = str_split(file_list, '_') %>% purrr::map(~ .x[which(!is.na(as.numeric(str_split(file_list, '_')[[1]])))]) %>% unlist()) %>%
      bind_cols('conv' = lapply(ix, getElement, 'message') %>%
                  unlist())
  })
  
  res_df <- lapply(1:length(res_df), function(ix) {
    #ix<-1
    res_df[[ix]]$dataset <- as.numeric(str_split(file_list[ix], '_')[[1]])[which(!is.na(as.numeric(str_split(file_list[ix], '_')[[1]])))[2]]
    res_df[[ix]]
  })
  
  res_df <- res_df %>%
    bind_rows()
  
  #expect_true(nrow(res_df) == length(file_list) * 10)
  
} else {
  
  res_df <- lapply(res_full, getElement, 'estpars') %>%
    bind_rows() %>%
    bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
                unlist()) %>%
    mutate(dataset = str_split(file_list, '_') %>% purrr::map(~ .x[which(!is.na(as.numeric(str_split(file_list, '_')[[1]])))]) %>% unlist()) %>%
    bind_cols('conv' = lapply(res_full, getElement, 'message') %>%
                unlist())
  
  expect_true(nrow(res_df) == length(file_list))
  
}

expect_true(all(is.finite(res_df$loglik)))

# Remove estimates that did not converge:
table(res_df$conv)

res_df <- res_df %>%
  filter(str_detect(conv, 'XTOL_REACHED')) %>%
  dplyr::select(-conv)
res_df1<-res_df %>% arrange(desc(loglik))
# Check whether at least one estimate remains for all synthetic datasets:
#expect_true(length(unique(res_df$dataset)) == 500)

# Keep only top fit for each synthetic dataset:
res_df <- res_df %>%
  group_by(dataset) %>%
  filter(loglik == max(loglik)) %>%
  ungroup()
#expect_true(nrow(res_df) == 500)

# Are all top fits within a reasonable range of log-likelihood values?
hist(res_df$loglik, breaks = 50)

# Get season-specific and shared parameter names:
unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')
shared_estpars <- res_df %>%
  dplyr::select(-contains(unit_estpars), -loglik, -dataset) %>%
  names()
true_estpars <- c(shared_estpars, unit_estpars)

# # Check that no states go below zero for any of the top-fit parameter sets:
# prof_lik <- FALSE
# 
# if (fit_canada) {
#   vir1 <- 'flu'
# } else {
#   vir1 <- 'flu_h1_plus_b'
# }
# 
# source('src/functions/setup_global_likelilhood.R')
# 
# traj_list <- lapply(1:length(seasons), function(ix) {
#   pars_temp <- res_df %>%
#     select(all_of(shared_estpars), contains(seasons[ix]))
#   names(pars_temp) <- true_estpars
# 
#   p_mat <- parmat(coef(po_list[[ix]]), nrep = nrow(pars_temp))
#   for (param in names(pars_temp)) {
#     p_mat[param, ] <- pars_temp %>% pull(param)
#   }
# 
#   trajectory(object = po_list[[ix]],
#              params = p_mat,
#              format = 'data.frame') %>%
#     select(!(H1_tot:H2_tot)) %>%
#     pivot_longer(X_SS:H2,
#                  names_to = 'state')
# })
# 
# expect_false(any(lapply(traj_list, function(ix) {
#   any(ix[, 'value'] < 0)
# }) %>%
#   unlist()))

# Calculate composite parameters:
if (fit_canada) {
  
  res_df_unit <- res_df %>%
    dplyr::select(contains('I') | contains('R') | dataset) %>%
    dplyr::select(-c(phi, rho1, rho2, phi1, phi2))
  
} else {
  
  res_df_unit <- res_df %>%
    dplyr::select(contains('I') | contains('R') | dataset) %>%
    dplyr::select(-c(phi, rho1, rho2))
  
  if (sens == 'sinusoidal_forcing') {
    res_df_unit <- res_df_unit %>%
      dplyr::select(-c(phi1, phi2))
  }
  
  if (sens == 'rhino_covar') {
    res_df_unit <- res_df_unit %>%
      dplyr::select(-beta_rhino)
  }
  
  if (sens == 'susc_plus_sev') {
    res_df_unit <- res_df_unit %>%
      dplyr::select(-c(theta_rho1, theta_rho2))
  }
  
}
#res_df_unit1<-res_df_unit
res_df_unit <- res_df_unit %>%
  pivot_longer(-c(loglik, dataset)) %>%
  # mutate(season = str_sub(name, 1, 6),
  #        name = str_sub(name, 8)) %>%
  separate(name, into = c("season", "name"), sep = "_") %>%
  pivot_wider(names_from = name,
              values_from = value) %>%
  mutate(`R10 + R120` = R10 + R120,
         `R20 + R120` = R20 + R120,
         R01 = Ri1 / (1 - (I10 + R10 + R120)),
         R02 = Ri2 / (1 - (I20 + R20 + R120))) %>%
  pivot_longer(Ri1:R02,
               names_to = 'parameter',
               values_to = 'value') %>%
  mutate(parameter = paste(season, parameter, sep = '_')) %>%
  dplyr::select(parameter:value)

res_df_shared <- res_df %>%
  dplyr::select(all_of(shared_estpars), loglik, dataset)

if (sens != 'no_int') {
  
  res_df_shared <- res_df_shared %>%
    mutate(delta2 = d2 * delta1)
  
}

res_df_shared <- res_df_shared %>%
  pivot_longer(-c(loglik, dataset),
               names_to = 'parameter',
               values_to = 'value') %>%
  dplyr::select(parameter:value)

if (sens != 'no_int') {
  
  res_df <- bind_rows(res_df_shared, res_df_unit) %>%
    mutate(parameter = factor(parameter, levels = c(shared_estpars, 'delta2', unique(res_df_unit$parameter))))
  
} else {
  
  res_df <- bind_rows(res_df_shared, res_df_unit) %>%
    mutate(parameter = factor(parameter, levels = c(shared_estpars, unique(res_df_unit$parameter))))
  
}

# Calculate 95% confidence intervals for each parameter:
# ci_res <- res_df %>%
#   group_by(parameter) %>%
#   # summarise(lower = min(value),
#   #           upper = max(value))
#   summarise(lower = HPDI(value, p = 0.95)[1],
#             upper = HPDI(value, p = 0.95)[2])
ci_res <- res_df %>%
  group_by(parameter) %>%
  summarise(lower = quantile(value, p = 0.025),
            upper = quantile(value, p = 0.975))


#==============================Read in MLEs and add to data frame:=======================
if (fit_canada) {
  
  mle <- read_rds('results/round2_fit/sens/canada/MLEs_flu.rds')[1, ]
  
  mle_unit <- mle %>%
    dplyr::select(contains('I') | contains('R')) %>%
    dplyr::select(-c(phi, rho1, rho2, phi1, phi2))
  
} else {
  
mle_ori <- read_rds(paste0('C:/Users/yfxu/SPH Dropbox/Yanfang Xu/Rt_project/virus_interference/code/data/self_define_outbreak/','MLE_',vir1,"_",vir2,'_round',which_round,'.rds'))
mle<-mle_ori[1,]

mle_unit <- mle[1,] %>%
    dplyr::select(contains('I') | contains('R')) %>%
    dplyr::select(-c(phi, rho1, rho2))

}

mle_unit <- mle_unit %>%
  pivot_longer(everything(), names_to = "name", values_to = "value") %>%
  mutate(
    season = str_remove(name, "_.*$"),
    name   = if_else(str_detect(name, "_"),
                     str_replace(name, "^[^_]*_", ""),
                     NA_character_)
  )%>%
  pivot_wider(names_from = name,
              values_from = value) %>%
  mutate(`R10 + R120` = R10 + R120,
         `R20 + R120` = R20 + R120,
         R01 = Ri1 / (1 - (I10 + R10 + R120)),
         R02 = Ri2 / (1 - (I20 + R20 + R120))) %>%
  pivot_longer(Ri1:R02,
               names_to = 'parameter',
               values_to = 'mle') %>%
  mutate(parameter = paste(season, parameter, sep = '_')) %>%
  dplyr::select(parameter:mle)

mle_shared <- mle %>%
  dplyr::select(all_of(shared_estpars)) 

if (sens != 'no_int') {
  
  mle_shared <- mle_shared %>%
    mutate(delta2 = d2 * delta1)
  
}

mle_shared <- mle_shared %>%
  pivot_longer(everything(),
               names_to = 'parameter',
               values_to = 'mle')

mle <- bind_rows(mle_shared, mle_unit)

ci_res <- ci_res %>%
  left_join(mle, by = 'parameter') 
 # select(parameter, mle, lower:upper)

# Generate tables of results

# mle_lo <-pivot_longer(mle_ori, cols = 1:ncol(mle_ori), names_to = "column", values_to = "param")
# mle_ci<-mle_lo %>% 
#   dplyr::group_by(column) %>%
#   # summarise(mle_lower = quantile(param, p = 0.025),
#   #           mle_upper = quantile(param, p = 0.975))
#   summarise(mle_lower = min(param, p = 0.025),
#             mle_upper = max(param, p = 0.975))
# no_best <- nrow(subset(mle_ori, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = df_use)))

# mle_mins <- sapply(mle_ori, min, na.rm = TRUE)
# mle_maxs <- sapply(mle_ori, max, na.rm = TRUE)
# ci95 <- function(x) {
#   n <- length(x)
#   mean_x <- mean(x, na.rm = TRUE)
#   se <- sd(x, na.rm = TRUE) / sqrt(n)
#   lower <- mean_x - 1.96 * se
#   upper <- mean_x + 1.96 * se
#   return(c(mean = mean_x, lower = lower, upper = upper))
# }
# ci95 <- function(x) {
#   x <- x[!is.na(x)] 
#   if (length(unique(x[!is.na(x)])) != 1) {
#   n <- length(x)
#   mean_x <- mean(x, na.rm = TRUE)
#   lower <- t.test(x)$conf.int[1]
#   upper <- t.test(x)$conf.int[2]
#   return(c(mean = mean_x, lower = lower, upper = upper))}
#   if (length(unique(x[!is.na(x)])) == 1) {
#     mean_x <- mean(x, na.rm = TRUE)
#     return(c(mean = mean_x, lower = NA, upper = NA))
#   }
# }
library(boot)
set.seed(1000)
ci95_boot <- function(x, R=1000) {
  # x<-mle_ori$d2
  boot_mean <- function(data, idx) mean(data[idx])
  b <- boot(x, boot_mean, R=R)
  cib<-boot.ci(b, type="perc")
  lower <- cib[["percent"]][4]
  upper <- cib[["percent"]][5]
  return(c(lower = lower, upper = upper))
}
if(vir2 =='rsv'){
  mle_mins <- t(sapply(mle_ori[,2:ncol(mle_ori)], ci95_boot))
  mle_mins<-rbind(rep(1,ncol(mle_mins)),mle_mins)
  rownames(mle_mins)[1]<-'rho1'
}else{
  mle_mins <- t(sapply(mle_ori, ci95_boot))
}

# mle_maxs <- sapply(mle_ori, ci95)

# mle_ci<-data.frame(parameter=colnames(mle_ori),mle_lower=mle_mins, mle_upper=mle_maxs) 
mle_ci<-data.frame(parameter=colnames(mle_ori),mle_lower=mle_mins[,"lower"], mle_upper=mle_mins[,"upper"]) 

# Check whether MLEs fall within CIs:
ci_res %>% filter(mle <= lower)
ci_res %>% filter(mle >= upper)

ci_res2<-left_join(ci_res,mle_ci,by=c('parameter'))

# ci_res2[1:12,c('lower','upper')]<-ci_res2[1:12,c('mle_lower','mle_upper')]
# ci_res<-ci_res2[,1:4]

ci_res2[,c('lower','upper')]<-ci_res2[,c('mle_lower','mle_upper')]
ci_res<-ci_res2[,1:4]

ci_res %>% filter(mle <= lower)
ci_res %>% filter(mle >= upper)
ci_res<-na.omit(ci_res)

ci_res1<-ci_res
ci_res1$parameter<-gsub("covid", "COVID", ci_res1$parameter, perl = TRUE)
ci_res1$parameter<-gsub("h1", "H1N1", ci_res1$parameter, perl = TRUE)
ci_res1$parameter<-gsub("h3", "H3N2", ci_res1$parameter, perl = TRUE)
ci_res1$parameter<-gsub("rsv", "RSV", ci_res1$parameter, perl = TRUE)

res_table <- ci_res1 %>%
  gt() %>%
  tab_header(title = 'Best-Fit Parameter Values') %>%
  fmt_number(columns = c(mle, lower, upper), decimals = 3, suffixing = TRUE)
print(res_table)

gtsave(res_table, filename = paste0('C:/Users/yfxu/SPH Dropbox/Yanfang Xu/Rt_project/virus_interference/progress_update/20251107/',
"MLE_with_mle_ci_",vir1,"_",vir2,'.pdf'))


#=====generate table with Description======
codebook<-read_excel("C:/Users/yfxu/SPH Dropbox/Yanfang Xu/Rt_project/virus_interference/code/data/parameters_codebook.xlsx")%>% as.data.frame()
colnames(codebook)[1]<-'name'
colnames(codebook)[2]<-'meaning'
share_code<-ci_res[which(ci_res$parameter %in% shared_estpars),]
share_code<-left_join(share_code,codebook,by=c("parameter"="name"))
ci_res %>% filter(mle <= lower)
ci_res %>% filter(mle >= upper)

share_code<-share_code %>%
  mutate(lower = ifelse(lower<=mle,lower,mle),
         upper = ifelse(upper>=mle, upper, mle))

write.csv(share_code,file=paste0('C:/Users/yfxu/SPH Dropbox/Yanfang Xu/Rt_project/virus_interference/progress_update/20251107/Table_',vir1,'_',vir2,'_withmleci.csv'),row.names = FALSE)



#======Plot range of fit values for each parameter:=======
res_df_unit <- res_df_unit %>%
  # mutate(season = str_sub(parameter, 1, 6),
  #        parameter = str_sub(parameter, 8)) %>%
  mutate(
    season = str_remove(parameter, "_.*$"),
    parameter   = if_else(str_detect(parameter, "_"),
                     str_replace(parameter, "^[^_]*_", ""),
                     NA_character_)
  )%>%
  mutate(parameter = factor(parameter, levels = c('Ri1', 'Ri2', 'R01', 'R02', 'I10', 'I20', 'R10', 'R20', 'R120', 'R10 + R120', 'R20 + R120')))%>%
  mutate(season=toupper(season))
res_df_unit$season<-gsub("H1(?!N1)", "H1N1", res_df_unit$season, perl = TRUE)
res_df_unit$season<-gsub("H3(?!N2)", "H3N2", res_df_unit$season, perl = TRUE)

p1 <- ggplot(data = res_df %>% filter(parameter %in% c(shared_estpars, 'delta2')),
             aes(x = value, y = after_stat(count)/nrow(res_df))) +
  geom_freqpoly(bins = 40) + facet_wrap(~ parameter, scales = 'free') +
  theme_classic() + labs(x = 'Parameter Value', y = 'Proportion of Fits', title = 'Shared Parameters')
p2 <- ggplot(data = res_df_unit, aes(x = value, y = after_stat(count)/nrow(res_df), col = season)) +
  geom_freqpoly(bins = 100) + facet_wrap(~ parameter, scales = 'free') +
  theme_classic() + scale_color_brewer(palette = 'Set1') +
  labs(x = 'Parameter Value', y = 'Proportion of Fits', title = 'Season-Specific Parameters')

print(p1)
print(p2)

ggsave(file=paste0("C:/Users/yfxu/SPH Dropbox/Yanfang Xu/Rt_project/virus_interference/progress_update/20251017/",vir1,'_',vir2,'_','seasonal_specifit_parameters_distribution.pdf'), 
       plot =p2, width =13, height =5, units = "in")

# Clean up:
#rm(list = ls())
