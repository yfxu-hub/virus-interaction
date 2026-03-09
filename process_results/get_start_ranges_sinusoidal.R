# ---------------------------------------------------------------------------------------------------------------------
# Get good-quality start values for the amplitude of sinusoidal forcing from fit climate forcing parameters
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(tidyverse)
library(viridis)
library(gridExtra)

# Read in MLEs:
mle <- read_rds('results/MLEs_flu_h1_plus_b.rds')[1, ]

# Read in data:
hk_dat <- read_rds('data/formatted/dat_hk_byOutbreak.rds')
dat_clim <- read_csv('data/formatted/clim_dat_hk_NORM.csv')

# Set duration of infection:
gamma1 <- 7/5
gamma2 <- 7/10

# Get results by season:
mle <- mle %>%
  select(eta_temp1:eta_ah2,
         contains('R10'),
         contains('R20'),
         contains('R120'),
         contains('Ri')) %>%
  pivot_longer(-c(eta_temp1:eta_ah2),
               names_to = 'parameter',
               values_to = 'mle') %>%
  mutate(season = str_sub(parameter, 2, 6),
         parameter = str_sub(parameter, 8)) %>%
  pivot_wider(names_from = parameter,
              values_from = mle) %>%
  select(eta_temp1:eta_ah2, Ri1:Ri2, R10:R120, season)

# Get vector of seasons:
seasons <- unique(mle$season)

# Calculate beta over time based on fit parameter values:
beta_t <- vector('list', length = length(seasons))
for (yr_index in 1:length(seasons)) {
  
  seas <- seasons[yr_index]
  
  dat_temp <- hk_dat[['h1_plus_b_rsv']] %>%
    filter(season == paste0('s', seas)) %>%
    inner_join(dat_clim,
               by = c('Year' = 'year',
                      'Week' = 'week')) %>%
    select(time, temp, ah)
  
  mle_temp <- mle %>%
    filter(season == seas)
  
  beta1_temp <- unlist(mle_temp['Ri1']) / (1.0 - (unlist(mle_temp['R10']) + unlist(mle_temp['R120']))) * exp(unlist(mle_temp['eta_ah1']) * dat_temp$ah + unlist(mle_temp['eta_temp1']) * dat_temp$temp) * gamma1
  beta2_temp <- unlist(mle_temp['Ri2']) / (1.0 - (unlist(mle_temp['R20']) + unlist(mle_temp['R120']))) * exp(unlist(mle_temp['eta_ah2']) * dat_temp$ah + unlist(mle_temp['eta_temp2']) * dat_temp$temp) * gamma2
  
  force1_temp <- exp(unlist(mle_temp['eta_ah1']) * dat_temp$ah + unlist(mle_temp['eta_temp1']) * dat_temp$temp)
  force2_temp <- exp(unlist(mle_temp['eta_ah2']) * dat_temp$ah + unlist(mle_temp['eta_temp2']) * dat_temp$temp)
  
  beta1_temp <- bind_cols(1:max(dat_temp$time), beta1_temp, force1_temp, seas)
  beta2_temp <- bind_cols(1:max(dat_temp$time), beta2_temp, force2_temp, seas)
  
  names(beta1_temp) <- c('time', 'beta1', 'force1', 'season')
  names(beta2_temp) <- c('time', 'beta2', 'force2', 'season')
  
  beta_temp <- beta1_temp %>%
    inner_join(beta2_temp,
               by = c('time', 'season')) %>%
    select(time, beta1, beta2, force1, force2, season)
  
  beta_t[[yr_index]] <- beta_temp
  
}
rm(yr_index, seas, dat_temp, mle_temp, beta1_temp, beta2_temp, force1_temp, force2_temp, beta_temp, hk_dat, dat_clim)

# Combine results into tibble:
beta_t <- bind_rows(beta_t)

# Fit sine wave to each seasonal beta and get range of amplitudes:
amp1_vec = amp2_vec = c()
for (seas in seasons) {
  
  beta_temp <- beta_t %>%
    filter(season == seas) %>%
    arrange(time)
  force1_temp <- beta_temp$force1
  force2_temp <- beta_temp$force2
  
  m1 <- lm(force1 ~ sin(2 * pi * time / 52.25) + cos(2 * pi * time / 52.25), data = beta_temp)
  m2 <- lm(force2 ~ sin(2 * pi * time / 52.25) + cos(2 * pi * time / 52.25), data = beta_temp)
  
  amp1 <- sqrt(unname(m1$coefficients[2]) ** 2 + unname(m1$coefficients[3]) ** 2)
  amp2 <- sqrt(unname(m2$coefficients[2]) ** 2 + unname(m2$coefficients[3]) ** 2)
  
  # phi1 <- atan2(unname(m1$coefficients[3]), unname(m1$coefficients[2])) * 52.25 / (2 * pi)
  # phi2 <- atan2(unname(m2$coefficients[3]), unname(m2$coefficients[2])) * 52.25 / (2 * pi)
  # 
  # par(mfrow = c(1, 2))
  # plot(beta_temp$time, beta_temp$force1, pch = 20)
  # lines(beta_temp$time, m1$coefficients[1] + amp1 * sin((2 * pi / 52.25) * (beta_temp$time + phi1)))
  # 
  # plot(beta_temp$time, beta_temp$force2, pch = 20)
  # lines(beta_temp$time, m2$coefficients[1] + amp2 * sin((2 * pi / 52.25) * (beta_temp$time + phi2)))
  
  amp1_vec <- c(amp1_vec, amp1)
  amp2_vec <- c(amp2_vec, amp2)
}

print(summary(amp1_vec))
print(summary(amp2_vec))

# Clean up:
rm(list = ls())
