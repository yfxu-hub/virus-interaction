# ---------------------------------------------------------------------------------------------------------------------
# Get start ranges for from round2 trajectory matching results (for profile likelihoods, bootstraps, or to rerun round2)
# ---------------------------------------------------------------------------------------------------------------------
rm(list = ls())
# Load libraries:
library(tidyverse)
library(testthat)
library(readxl)
# Set directory where results from round2 fits are stored:
#res_dir <- 'results/round2_fit/round2_3_fluH1_plus_B/'
#res_dir<-"/lustre1/g/sph_timtsang/yfxu/virus_interference/src_simulation/results_round2/fluab/round2_fit"
# res_dir<-"/lustre1/g/sph_timtsang/yfxu/virus_interference/src_simulation/multipathogons/fitbyoutbreak_round2/covid_h1/"
setwd("/home/yfxu/virus_interference")
# Which round of fits?:
#which_round <- str_split(res_dir, '_')[[1]][which(!is.na(as.numeric(str_split(res_dir, '_')[[1]])))]
which_round<-3

vir1_name<-'h1'
vir2_name<-'rsv'

res_dir<-sprintf('/lustre1/g/sph_timtsang/yfxu/virus_interference/src_simulation/multipathogons/self_define_outbreak/%s_%s_round2/',
vir1_name,vir2_name)

# res_dir<-sprintf('/lustre1/g/sph_timtsang/yfxu/virus_interference/src_simulation/multipathogons/self_define_outbreak/Sinusoidal_Forcing/%s_%s_round2/',
                 # vir1_name,vir2_name)
#res_dir<-'/lustre1/g/sph_timtsang/yfxu/virus_interference/src_simulation/multipathogons/season_defined_by_week45/round2/covid_h1/'

#outbreak<-read_excel("/home/yfxu/virus_interference/data/multi_pathogons_outbreak.xlsx")%>%as.data.frame()
#outbreak<-read_excel("/home/yfxu/virus_interference/data/multi_pathogons_outbreak_allperiod.xlsx")%>%as.data.frame()
outbreak<-read_excel("/home/yfxu/virus_interference/data/multi_pathogons_outbreak_1_with_flursv.xlsx")%>%as.data.frame()
match_v<-intersect(grep(vir1_name, outbreak$season, value = TRUE),grep(vir2_name, outbreak$season, value = TRUE))
outbreak<-outbreak[which(outbreak$season %in% match_v),]
seasons <- outbreak$season

vir1<-vir1_name
vir2<-vir2_name

#seasons<-c('s22-23','s23-24','s24-25')
#seasons<-c(paste0(vir1,'-',vir2))

# Are results from a sensitivity analysis?:
if (str_detect(res_dir, 'sinusoidal')) {
  sens <- 'sinusoidal_forcing'
} else if (str_detect(res_dir, 'h3_covar')) {
  sens <- 'h3_covar'
} else if (str_detect(res_dir, 'less_circ_h3')) {
  sens <- 'less_circ_h3'
} else if (str_detect(res_dir, 'no_rsv_immune')) {
  sens <- 'no_rsv_immune'
} else if (str_detect(res_dir, 'no_ah')) {
  sens <- 'no_ah'
} else if (str_detect(res_dir, 'rhino_covar')) {
  sens <- 'rhino_covar'
} else if (str_detect(res_dir, 'no_int')) {
  sens <- 'no_int'
} else if (str_detect(res_dir, 'susc_plus_sev')) {
  sens <- 'susc_plus_sev'
} else {
  sens <- 'main'
}

# Check whether Canada or Hong Kong:
if (str_detect(res_dir, 'canada')) {
  fit_canada <- TRUE
  sens <- 'sinusoidal_forcing'
} else {
  fit_canada <- FALSE
}

# Check that directory for storing results exists, and create if not:
if (!dir.exists('results/')) {
  dir.create('results/')
}
if (!dir.exists('results/round2_CIs/')) {
  dir.create('results/round2_CIs/')
}

if (sens == 'main') {
  
  new_dir <- paste0('results/round2_CIs/from_2_', which_round, '/')
  if (!dir.exists(new_dir)) {
    dir.create(new_dir)
  }
  
  if (str_detect(res_dir, 'age_structured')) {
    
    if(!dir.exists('results/round2_CIs/sens/')) {
      dir.create('results/round2_CIs/sens/')
    }
    if(!dir.exists('results/round2_CIs/sens/age_structured/')) {
      dir.create('results/round2_CIs/sens/age_structured/')
    }
    
    new_dir <- paste0('results/round2_CIs/sens/age_structured/from_2_', which_round, '/')
    if (!dir.exists(new_dir)) {
      dir.create(new_dir)
    }
  }
  
} else {
  
  if(!dir.exists('results/round2_CIs/sens/')) {
    dir.create('results/round2_CIs/sens/')
  }
  if(!dir.exists(paste0('results/round2_CIs/sens/', sens, '/'))) {
    dir.create(paste0('results/round2_CIs/sens/', sens, '/'))
  }
  
  new_dir <- paste0('results/round2_CIs/sens/', sens, '/from_2_', which_round, '/')
  if (!dir.exists(new_dir)) {
    dir.create(new_dir)
  }
  
  if (fit_canada) {
    new_dir <- 'results/round2_CIs/sens/canada/'
    if (!dir.exists(new_dir)) {
      dir.create(new_dir)
    }
    
    new_dir <- paste0('results/round2_CIs/sens/canada/from_2_', which_round, '/')
    if (!dir.exists(new_dir)) {
      dir.create(new_dir)
    }
  }
  
}

# Function to read in and format results:
load_and_format_mega_results <- function(filepath) {
  #filepath<-res_dir
  # Get list of results files:
  res_files <- list.files(path = filepath, full.names = TRUE)
  # res_files<-sprintf('/lustre1/g/sph_timtsang/yfxu/virus_interference/src_simulation/results_round2/fluab/round2_fit/res_flu_h1_plus_b_susc_1_%d_PARALLEL.rds',
  #                    which_round)
  #selected_files <- res_files[grepl(paste(which_round,"_PARALLEL\\.rds$",sep=''), res_files)]
  selected_files <- res_files[grepl(paste("res_",vir1_name,sep=''), res_files)]
  selected_files <- selected_files[grepl(paste(which_round,"_PARALLEL\\.rds$",sep=''), selected_files)]
  res_files<-selected_files
  
  # Read in results:
  res_full = list()
  for (i in seq_along(res_files)) {
    res_full[[i]] <- read_rds(res_files[[i]])
  }
  
  # Get parameter estimates and log-likelihoods:
  if (str_detect(res_files[1], 'PARALLEL')) {
    
    res_full <- do.call('c', res_full)
    num_errors <- length(which(res_full == 'error'))
    
    if (num_errors > 0) {
      res_full <- res_full[-which(res_full == 'error')]
    }
    
  }
  
  pars_df <- lapply(res_full, getElement, 'estpars') %>%
    bind_rows() %>%
    bind_cols('loglik' = lapply(res_full, getElement, 'll') %>%
                unlist()) %>%
    bind_cols('message' = lapply(res_full, getElement, 'message') %>%
                unlist())
  
  # if (str_detect(res_files[1], 'PARALLEL')) {
  #   expect_true(nrow(pars_df) == (length(res_files) * 25) - num_errors)
  # } else {
  #   expect_true(nrow(pars_df) == length(res_files))
  # }
  expect_true(all(is.finite(pars_df$loglik)))
  
  # Keep only top results:
  pars_df <- pars_df %>%
    arrange(desc(loglik))
  
  df_use <- pars_df %>% select(-c(loglik, message)) %>% names() %>% length()
  no_best <- nrow(subset(pars_df, 2 * (max(loglik) - loglik) <= qchisq(p = 0.95, df = df_use)))
  print(table(pars_df$message))
  
  # If only one parameter set is in this range, MLE has not yet been reached; take top 5% of fits instead:
  if (no_best == 1) {
    is_mle <- FALSE
    print('MLE not reached!')
    no_best <- 25
  } else {
    is_mle <- TRUE
  }
  print(no_best)
  
  # Get tibble of top fits:
  pars_top <- pars_df[1:no_best, ]
  print(summary(pars_top$loglik))
  
  # Remove where no convergence occurs:
  pars_top <- pars_top %>%
    filter(!str_detect(message, 'maxtime')) %>%
    select(-message)
  
  # If none remaining, print warning:
  if(nrow(pars_top) == 0) {
    print('No convergence among best-fit runs!')
    
    no_best <- 25
    pars_top <- pars_df[1:no_best, ]
    print(summary(pars_top$loglik))
    
    is_mle <- FALSE
    
    pars_top <- pars_top %>%
      filter(!str_detect(message, 'maxtime')) %>%
      select(-message)
  }
  print(no_best)
  print(summary(pars_top$loglik))
  
  # Store original results before removing unrealistic parameter values:
  pars_top_orig <- pars_top
  
  # # Remove where d2 > 10 and theta_lambda2 != 1.0:
  # pars_top <- pars_top %>%
  #   filter(!(d2 > 10.0 & theta_lambda2 < 0.99))
  
  # Set unrealistic values to NA:
  if (!str_detect(filepath, 'no_int')) {
    pars_top$delta1[pars_top$delta1 > 7.0] <- NA
    pars_top$d2[pars_top$d2 > 10.0] <- NA
  }
  pars_top$rho1[pars_top$rho1 == 1.0] <- NA
  pars_top$rho2[pars_top$rho2 == 1.0] <- NA
  
  # Since phi=0 is equivalent to phi=52.25, don't use full range; transform so that we can select from only best-supported range:
  if ('phi' %in% names(pars_top)) {
    
    par(mfrow = c(2, 1))
    #hist(pars_top$phi, breaks = 50)
    pars_top <- pars_top %>%
      mutate(phi = if_else(phi < 1, phi + 52.25, phi))
    #hist(pars_top$phi, breaks = 50)
    
  }
  
  # If using sinusoidal forcing, do the same for phi1 and phi2:
  if ('phi1' %in% names(pars_top)) {
    
    par(mfrow = c(2, 2))
    # hist(pars_top$phi1, breaks = 50)
    # hist(pars_top$phi2, breaks = 50)
    
    pars_top <- pars_top %>%
      mutate(phi1 = if_else(phi1 < 1, phi1 + 52.25, phi1),
             phi2 = if_else(phi2 < 1, phi2 + 52.25, phi2))
    
    # hist(pars_top$phi1, breaks = 50)
    # hist(pars_top$phi2, breaks = 50)
    
  }
  
  # Return formatted results:
  return(list(pars_top, pars_top_orig, is_mle))
  
}

# Read in results:
res <- load_and_format_mega_results(res_dir)

# Check that best-fit parameter values do not lead trajectories to drop below 0:
res_orig <- res[[2]]

# if (fit_canada) {
#   vir1 <- 'flu'
# } else {
#   vir1 <- 'flu_h1_plus_b'
# }
prof_lik <- FALSE

unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10', 'R20', 'R120')
if (sens == 'no_rsv_immune') {
  unit_estpars <- c('Ri1', 'Ri2', 'I10', 'I20', 'R10')
}

shared_estpars <- res_orig %>% select(!contains(unit_estpars) & !'loglik') %>% names()
true_estpars <- c(shared_estpars, unit_estpars)


source('/home/yfxu/virus_interference/code/setup_global_likelilhood.R')

traj_list <- lapply(1:length(seasons), function(ix) {
  ix<-1
  pars_temp <- res_orig %>%
    select(all_of(shared_estpars), contains(paste(seasons[ix],"_",sep='')))
  names(pars_temp) <- true_estpars
  
  p_mat <- parmat(coef(po_list[[ix]]), nrep = nrow(pars_temp))
  for (param in names(pars_temp)) {
    p_mat[param, ] <- pars_temp %>% pull(param)
  }
  
  trajectory(object = po_list[[ix]],
             params = p_mat,
             format = 'data.frame') %>%
    select(!(H1_tot:H2_tot)) %>%
    pivot_longer(X_SS:H2,
                 names_to = 'state')
})

expect_false(
  any(
    lapply(traj_list, function(ix) any(ix[, "value"] < 0, na.rm = TRUE)) %>%
      unlist(),
    na.rm = TRUE
  )
)

expect_false(any(lapply(traj_list, function(ix) {
  any(ix[, 'value'] < 0)
}) %>%
  unlist()))

# Are results the MLE?
is_mle <- res[[3]]

res_dir_comp <- str_split(res_dir, '_')[[1]]
res_dir_comp[which(!is.na(as.numeric(res_dir_comp)))] <- as.character(as.numeric(which_round) - 1)
res_dir_prev <- paste(res_dir_comp, collapse = '_')

if (sens == 'susc_plus_sev') {
  res_dir_prev <- paste0('results/round2_fit/round2_', as.numeric(which_round) - 1, '_fluH1_plus_B/')
}

is_mle_prev <- try(
  load_and_format_mega_results(res_dir_prev)[[3]]
)

if (inherits(is_mle_prev, 'try-error')) {
  is_mle_prev <- FALSE
}

# Save MLEs:
if (is_mle & is_mle_prev) {
  
  res_orig <- res_orig %>%select(-loglik)
  file_name<-paste("MLE_",vir1_name,"_",vir2_name,"_round",which_round,".rds",sep='')
  
  if (str_detect(res_dir, 'sens')) {
    
    if (str_detect(res_dir, 'canada')) {
      write_rds(res_orig, file = 'results/round2_fit/sens/canada/MLEs_flu.rds')
    } else if (str_detect(res_dir, '/us/')) {
      write_rds(res_orig, file = 'results/round2_fit/sens/us/MLEs_flu.rds')
    } else {
      write_rds(res_orig, file = paste0(paste(str_split(res_dir, '/')[[1]][1:(length(str_split(res_dir, '/')[[1]]) - 2)], collapse = '/'), '/MLEs_flu_h1_plus_b.rds')
                )
    }
    
  } else if (str_detect(res_dir, 'age_structured')) {
    write_rds(res_orig, file = 'results/age_structured_SA/MLEs_flu_h1_plus_b.rds')
  } else {
    write_rds(res_orig, #file = 'results/MLEs_flu_h1_plus_b.rds'
              file=paste("results/",file_name,sep=''))
  }
  
}

# Get results for determining start ranges:
res <- res[[1]] %>%
  select(-loglik)

# Get minimum and maximum start values:
ci_start <- as.data.frame(rbind(summarise(res, across(.cols = everything(), \(x) min(x, na.rm = TRUE))),
                                summarise(res, across(.cols = everything(), \(x) max(x, na.rm = TRUE)))))

# Any parameters where minimum and maximum are equal?:
no_range <- c()
for (i in 1:ncol(ci_start)) {
  if (identical(ci_start[1, i], ci_start[2, i])) {
    no_range <- c(no_range, i)
  }
}
if (length(no_range) > 0) {
  print('Warning: Some parameters have same min and max values!')
  print(names(ci_start)[no_range])
}

# Possible that d2 ranges are missing if all top fits were > 10; if so, replace:
if (any(ci_start == Inf)) {
  if (any(ci_start$d2 == Inf)) {
    ci_start$d2 <- c(0, 10)
  } else if (any(ci_start$delta1 == Inf)) {
    ci_start$delta1 <- c(7 / 60, 7)
  } else {
    print('Range is NA for some other parameter.')
  }
}

# Check that sums of initial conditions can't sum to >1:
init_cond_estpars <- c('I10', 'I20', 'R10', 'R20', 'R120')

sums <- ci_start %>%
  mutate(minmax = c('min', 'max')) %>%
  select(contains(init_cond_estpars), minmax) %>%
  pivot_longer(-minmax) %>%
  mutate(season = sapply(strsplit(name, "_"), `[`, 1))%>%
  #mutate(season = str_sub(name, 1, 10)) %>%
  group_by(season, minmax) %>%
  summarise(sum = sum(value))

if (any(sums %>% filter(minmax == 'min') %>% pull(sum) > 1.0)) {
  print('Lower bounds sum to more than 1!')
}


if (any(sums %>% filter(minmax == 'max') %>% pull(sum) > 1.0)) {
  seasons <- sums %>%
    filter(minmax == 'max',
           sum > 1.0) %>%
    pull(season)
  
  for (yr in seasons) {
    print(yr)
    # Reduce upper bounds proportionally:
    orig_upper_bounds <- ci_start[2, ] %>%
      # select(contains(c("I10","I20",'R10', 'R20', 'R120'))) %>%
      select(contains(c('R10', 'R20', 'R120'))) %>%
      select(contains(yr))
    red_needed <- sums %>%
      filter(season == yr,
             minmax == 'max') %>%
      pull(sum) - 0.9999999
    new_upper_bounds <- orig_upper_bounds - (red_needed * (orig_upper_bounds / sum(orig_upper_bounds)))
    #red_needed <- min(
     # sums %>% filter(season == yr, minmax == 'max') %>% pull(sum) - 0.9999999,
     # sum(orig_upper_bounds) - sum(orig_lower_bounds) - 1e-6
    #)
    
    # Ensure upper bounds still greater than lower:
    orig_lower_bounds <- ci_start[1, ] %>%
      # select(contains(c("I10","I20",'R10', 'R20', 'R120'))) %>%
      select(contains(c('R10', 'R20', 'R120'))) %>%
      select(contains(yr))
    
    if (!all(new_upper_bounds > orig_lower_bounds)) {
      valid_idx <- which(orig_lower_bounds < new_upper_bounds)
      red_needed_sub <- -sum(new_upper_bounds[valid_idx]) + sum(orig_lower_bounds[valid_idx]) 
      red_needed_sub <- max(red_needed_sub, 0)
      
      new_upper_bounds_try <- orig_upper_bounds
      new_upper_bounds_try[-which(orig_lower_bounds >= new_upper_bounds)] <- orig_upper_bounds[-which(orig_lower_bounds >= new_upper_bounds)] - 
        (red_needed * (orig_upper_bounds[-which(orig_lower_bounds >= new_upper_bounds)] / sum(orig_upper_bounds[-which(orig_lower_bounds >= new_upper_bounds)])))
      
      #new_upper_bounds[-which(orig_lower_bounds >= new_upper_bounds)] <- new_upper_bounds_try[new_upper_bounds > orig_lower_bounds]
      new_upper_bounds<- new_upper_bounds_try
     
      #rm(new_upper_bounds_try)
    }
    #new_upper_bounds <- orig_upper_bounds - (red_needed * (orig_upper_bounds / sum(orig_upper_bounds)))
    #new_upper_bounds[new_upper_bounds < orig_lower_bounds] <- orig_lower_bounds[new_upper_bounds < orig_lower_bounds] + 1e-6
    expect_true(all(new_upper_bounds > 0))
    expect_true(all(new_upper_bounds > orig_lower_bounds))
    
    # Check that upper bounds now sum to 1 or less:
    ci_start[2, which(str_detect(names(ci_start), yr) &
                        (#str_detect(names(ci_start), 'I10') |
                           #str_detect(names(ci_start), 'I20') |
                           str_detect(names(ci_start), 'R10') |
                           str_detect(names(ci_start), 'R20') |
                           str_detect(names(ci_start), 'R120')))] <- new_upper_bounds
    
    expect_lt(ci_start[2, ] %>% select(contains(yr)) %>% select(contains(init_cond_estpars)) %>% sum(), 1.0)
    
  }
  rm(yr)
}

# Write start ranges to file:
file_name<-paste("round2CI_fitforselfdefineoutbreak_",vir1_name,"_",vir2_name,".rds",sep='')
#write_rds(ci_start, file = paste0(new_dir, 'round2CI_startvals_frompre1round.rds'))
write_rds(ci_start, file = paste0(new_dir, file_name))
# Clean up:
#rm(list = ls())

#check the mle and ci

for(i in 1:nrow(res_orig)){
  ci_mle<-rbind(ci_start,res_orig[i,])
  df <- tibble(
    column = colnames(ci_mle),
    lower = as.numeric(ci_mle[1, ]),
    upper = as.numeric(ci_mle[2, ]),
    v = as.numeric(ci_mle[3,])
  ) %>%
    mutate(ok = v >= lower & v <= upper)
  
  if(length(df$ok[which(df$ok ==TRUE)])==length(df$ok)){
    print(i)
  }
}





