# ---------------------------------------------------------------------------------------------------------------------
# Set up model of flu/RSV transmission
# Note: Adapted from code of Dr. Matthieu Domenech de Celles
# ---------------------------------------------------------------------------------------------------------------------

# Load packages:
library(tidyverse)
library(pomp)
library(viridis)
library(gridExtra)
library(zoo)
library(readxl)

# Load required functions:
source('/home/yfxu/virus_interference/code/functions_flu_RSV.R')
source('/home/yfxu/virus_interference/code/test_code.R')

# Fit to observed data, not synthetic data from age-structured model:
if (!exists('age_structured')) {
  age_structured <- FALSE
}

# ---------------------------------------------------------------------------------------------------------------------

# Load and format data

# Read in data:
# hk_dat <- read_rds('data/formatted/dat_hk_byOutbreak.rds')
# can_dat <- read_csv('data/formatted/dat_canada.csv')
# us_dat <- read_rds('data/formatted/dat_us_byRegion.rds')

# Load libraries:
library(nloptr)

#load data
#outbreak<-read_excel("/home/yfxu/virus_interference/data/multi_pathogons_outbreak.xlsx")%>%as.data.frame()
outbreak<-read_excel("/home/yfxu/virus_interference/data/multi_pathogons_outbreak_1_with_flursv.xlsx")%>%as.data.frame()

df1 <- read_csv("/home/yfxu/virus_interference/data/multi_virus_comb_case_ili.csv")%>%as.data.frame()

colnames(df1)[which(colnames(df1) == "Severe acute respiratory syndrome coronavirus 2")]<-'covid'
colnames(df1)[which(colnames(df1) == "H1")]<-'h1'
colnames(df1)[which(colnames(df1) == "H3")]<-'h3'
colnames(df1)[which(colnames(df1) == "Type B")]<-'flub'
colnames(df1)[which(colnames(df1) == "RSV")]<-'rsv'

df1<-df1[5:nrow(df1),]

# df1<-df1[5:nrow(df1),]
df1[] <- lapply(df1, function(x) {
  if (is.numeric(x)) na.approx(x, na.rm = FALSE) else x
})

#impute the missing data
# valid <- !is.na(df1$`Severe acute respiratory syndrome coronavirus 2`)
# spline_result <- spline(x = df1$start_date[valid],
#                         y = df1$`Severe acute respiratory syndrome coronavirus 2`[valid],
#                         xout = df1$start_date)
# df1$`Severe acute respiratory syndrome coronavirus 2` <- spline_result$y
df1$n_T<-df1$`No. of specimen tested`
df1$n_P1<-df1[,vir1_name]
df1$n_P2<-df1[,vir2_name]
df1$i_ILI<-df1$PMP/1000
df1$pop<-rep(7531800,nrow(df1))
df1$temp<-as.vector(scale(df1$avg_temperature))
df1$ah<-as.vector(scale(df1$avg_abs_humidity))
#df1$flu_h1_plus_b<-df1$H3+df1$`Type B`
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


#==========set season=========
#spefici outbreak
match_v<-intersect(grep(vir1_name, outbreak$season, value = TRUE),grep(vir2_name, outbreak$season, value = TRUE))
outbreak<-outbreak[which(outbreak$season %in% match_v),]

df1$season<-rep(NA,length=nrow(df1))
for(r in 1:nrow(outbreak)){
  df1$season[which(df1$start_date >= as.Date(outbreak[r,"start_date"]) & df1$start_date <= as.Date(outbreak[r,"end_date"]))]<-outbreak[r,1]
}
df1<-df1 %>% group_by(season) %>% mutate(time = row_number())

#set the season sefine by 45 week
# df1 <- df1 %>%
#   mutate(
#     iso_year = isoyear(start_date),              
#     iso_week = isoweek(start_date),             
#     season_start = if_else(
#       iso_week >= 45,
#       iso_year,
#       iso_year - 1
#     ),                                     
#     season = paste0(
#       "s",
#       str_sub(season_start, 3, 4),         
#       "-",
#       str_sub(season_start + 1, 3, 4)      
#     )
#   ) %>%
#   select(-iso_year, -iso_week, -season_start)

# df1$season<-rep(paste0(vir1,"-",vir2),length=nrow(df1))
# 
# df1<-na.omit(df1)
# df1<-df1 %>% group_by(season) %>% mutate(time = row_number())

hk_dat<-df1

dat_pomp<-hk_dat %>%
  filter(season == yr)
nrow_check <- nrow(dat_pomp)



# # Get H3 incidence data (as covariate):
#   dat_h3 <- hk_dat %>%
#    # rename('i_ILI' = 'GOPC') %>%
#     mutate(i_ILI = i_ILI / 1000,
#            h3_inc_raw = (n_P1 / n_T) * i_ILI,
#            h3_inc = scale(h3_inc_raw)[, 1]) %>%
#     select('time', 'year', 'Week', 'season','h3_inc')
#     #select(time:Year, Week:season, h3_inc)
#   expect_equal(mean(dat_h3$h3_inc, na.rm = TRUE), 0)
#   #expect_equal(sd(dat_h3$h3_inc, na.rm = TRUE), 1)
# 
#   dat_h3 <- dat_h3 %>%
#     mutate(h3_inc_lag1 = lag(h3_inc, 1),
#            h3_inc_lag2 = lag(h3_inc, 2),
#            h3_inc_lagd = lag(h3_inc, 15))
# 
#   dat_pomp <- dat_pomp %>%
#     inner_join(dat_h3,
#                by = c('time', 'year', 'Week', 'season')) %>%
#     #select(time:year, Week:season, n_T:i_ILI, pop:h3_inc_lagd) %>%
#     mutate(h3_inc = if_else(is.na(h3_inc), 0, h3_inc),
#            h3_inc_lag1 = if_else(is.na(h3_inc_lag1), 0, h3_inc_lag1),
#            h3_inc_lag2 = if_else(is.na(h3_inc_lag2), 0, h3_inc_lag2),
#            h3_inc_lagd = if_else(is.na(h3_inc_lagd), 0, h3_inc_lagd))
#   expect_true(nrow(dat_pomp) == nrow_check)
#   rm(dat_h3)
# 
#   dat_pomp <- dat_pomp %>%
#     select(-h3_inc_lag1, -h3_inc_lag2)
#  

# If no data for this season, skip:
if (nrow(dat_pomp) > 0) {

  # Plot data:
  if (debug_bool) {
    # Plot ILI incidence:
    p1 <- ggplot(data = dat_pomp, aes(x = time, y = i_ILI)) + geom_line() +
      labs(x = 'Time (Weeks)', y = 'ILI Incidence Rate (per 1000 Consultations)') +
      theme_classic()
    print(p1)

    if (!fit_canada) {

      p2 <- ggplot(data = dat_pomp %>%
                     pivot_longer(n_P1:n_P2, names_to = 'virus', values_to = 'n_pos'),
                   aes(x = time, y = n_pos / n_T, color = virus)) +
        geom_line() + labs(x = 'Time (Weeks)', y = 'Positivity Fraction') +
        theme_classic()

    } else {

      p2 <- ggplot(data = dat_pomp %>%
                     mutate(n_P1 = n_P1 / n_T1, n_P2 = n_P2 / n_T2) %>%
                     pivot_longer(n_P1:n_P2, names_to = 'virus', values_to = 'perc_pos'),
                   aes(x = time, y = perc_pos, color = virus)) +
        geom_line() + labs(x = 'Time (Weeks)', y = 'Positivity Fraction') +
        theme_classic()

    }
    print(p2)
  }

  # ---------------------------------------------------------------------------------------------------------------------

  # Create pomp model and run basic model checks

  # Create model:
  if (fit_canada) {

    resp_mod <- create_SITRxSITR_mod(dat = dat_pomp,
                                     Ri1_max = Ri_max1,
                                     Ri2_max = Ri_max2,
                                     d2_max = d2_max,
                                     debug_bool = debug_bool,
                                     sens = sens,
                                     loc = 'canada')

  } else {

    resp_mod <- create_SITRxSITR_mod(dat = dat_pomp,
                                     Ri1_max = Ri_max1,
                                     Ri2_max = Ri_max2,
                                     d2_max = d2_max,
                                     debug_bool = debug_bool,
                                     sens = sens,
                                     loc = 'hk')

  }

  # Check transformations:
  check_transformations(resp_mod)

  # Check parameters:
  check_params(resp_mod)

  # Check initial conditions:
  expect_true(all.equal(sum(rinit(resp_mod)), as.numeric(coef(resp_mod, 'N'))))

  # Check constant population size:
  check_correct_N_CONST(resp_mod, unique(dat_pomp$pop))

  # Run deterministic simulation:
  sim_determ <- trajectory(object = resp_mod, format = 'data.frame') %>%
    dplyr::select(H1:.id) %>%
    pivot_longer(H1:H2, names_to = 'Vir', values_to = 'Inc')
  p3 <- ggplot(data = sim_determ, aes(x = time, y = Inc, group = Vir, col = Vir)) +
    geom_line() + geom_point() +
    labs(x = 'Time (Weeks)', y = 'Incidence', col = 'Virus') +
    theme_classic()
  if (debug_bool) print(p3)

  # Run stochastic simulation and check that obs never more than n_samp:
  p4 <- check_obs_lessthan_samples(resp_mod, test_diff = fit_canada)
  if (debug_bool) print(p4)

  # Check that measurement density model works:
  ll <- logLik(traj_objfun(data = resp_mod))
  if (debug_bool) print(ll)

  # Check that dynamics are independent when there is no interaction:
  p5 <- check_independent_dynamics(resp_mod)
  if (debug_bool) print(p5)

  # Quick exploration of how interactions impact dynamics:
  p6 <- quick_explore_interaction(resp_mod, c(0, 0.2, 0.4, 0.6, 0.8, 1.0), n_sim = 5)
  if (debug_bool) do.call('grid.arrange', c(p6, ncol = 2))

  # Clean up:
  rm(sim_determ, p3, p4, p5, p6, ll)

}

