# ---------------------------------------------------------------------------------------------------------------------
# Functions to help evaluate model results
# ---------------------------------------------------------------------------------------------------------------------

run_sim <- function(pomp_object, seas, mle, shared_estpars, unit_estpars, model_type, return_obs = FALSE, return_data = FALSE, n_sim = 10, analysis = 'basic') {
  # Function to generate deterministic or stochastic simulations from the model at the MLE
  # params pomp_object: The pomp model object used to run the model
  # params seas: The season to be simulated
  # params mle: The maximum likelihood estimates of all parameter values
  # params shared_estpars: A vector containing the names of all shared parameters
  # params unit_estpars: A vector containing the names of all season-specific parameters
  # params model_type: Should simulations be taken from the deterministic or the stochastic model?
  # params return_obs: Should average number of observed cases also be included? (model_type "deterministic" only)
  # params n_sim: Number of stochastic simulations to be run (model type "stochastic" only)
  # params analysis: Should any special analyses be performed? (options: 'paf')
  # returns: Tibble of simulated values for flu and RSV in a given season
  
  # If mle is a list, get alternative interaction parameters:
  if (!is_tibble(mle)) {
    alt_int_params <- mle[[2]]
    mle <- mle[[1]]
  }
  
  # Get all estpars:
  true_estpars <- c(shared_estpars, unit_estpars)
  
  # Get parameter values:
  pars <- mle[1, ] %>%
    select(all_of(shared_estpars),
           contains(seas))
  names(pars)[(length(names(pars)) - 6):length(names(pars))] <- unit_estpars
  
  if (exists('alt_int_params')) {
    pars[c('theta_lambda1', 'theta_lambda2', 'delta1', 'd2')] <- alt_int_params[1, ] %>%
      select(all_of(c('theta_lambda1', 'theta_lambda2', 'delta1', 'd2')))
  }
  
  coef(pomp_object, true_estpars) <- pars
  
  # Simple simulation or additional analysis?:
  if (analysis == 'paf') {
    
    if (model_type == 'deterministic') {
      
      if (!return_obs) {
        
        # Generate matrix of model parameters for fitting:
        param_mat_temp <- parmat(coef(pomp_object), nrep = 5)
        
        # Remove impact of flu from some parameter sets:
        param_mat_temp['theta_lambda1', 2] <- 1.0
        param_mat_temp['I10', 3] <- 0
        
        # Remove impact of RSV from some parameter sets:
        param_mat_temp['theta_lambda2', 4] <- 1.0
        param_mat_temp['I20', 5] <- 0
        
        # Simulate using deterministic model:
        out_temp <- trajectory(pomp_object, params = param_mat_temp, format = 'data.frame') %>%
          mutate(season = seas) %>%
          select(time:season, H1:H2)
        
        # Check that removal of interaction works as expected:
        expect_true(all.equal(out_temp %>% filter(.id == 2) %>% pull(H2), out_temp %>% filter(.id == 3) %>% pull(H2)))
        expect_true(all.equal(out_temp %>% filter(.id == 4) %>% pull(H1), out_temp %>% filter(.id == 5) %>% pull(H1)))
        
        # Remove unneeded simulations:
        out_temp <- out_temp %>%
          filter(.id %in% c(1, 3, 5))
        
      }
    }
    
  } else {
    
    # Run simple simulation
    
    # Deterministic or stochastic?:
    if (model_type == 'deterministic') {
      
      # Get trajectory at MLE:
      out_temp <- trajectory(pomp_object, format = 'data.frame') %>%
        mutate(season = seas) %>%
        select(time, season, H1:H2)
      
      if ('n_T1' %in% rownames(pomp_object@covar@table)) {
        
        out_temp <- out_temp %>%
          cbind(pomp_object@covar@table[1, ]) %>%
          cbind(pomp_object@covar@table[2, ]) %>%
          cbind(pomp_object@covar@table[3, ]) %>%
          as_tibble()
        names(out_temp)[5:7] <- c('i_ILI', 'n_T1', 'n_T2')
        
      } else {
        
        out_temp <- out_temp %>%
          cbind(pomp_object@covar@table[1, ]) %>%
          cbind(pomp_object@covar@table[2, ]) %>%
          as_tibble()
        names(out_temp)[5:6] <- c('i_ILI', 'n_T')
        
      }
      
      if (return_obs) {
        
        # Calculate mean number of observed cases:
        rho1 <- as.numeric(pars['rho1'])
        rho2 <- as.numeric(pars['rho2'])
        
        if ('alpha' %in% names(pars)) {
          
          alpha <- as.numeric(pars['alpha'])
          phi <- as.numeric(pars['phi'])
          
          rho1_w <- rho1 * (1.0 + alpha * cos(((2 * pi) / 52.25) * (out_temp$time - phi))) * out_temp$H1 / out_temp$i_ILI
          rho2_w <- rho2 * (1.0 + alpha * cos(((2 * pi) / 52.25) * (out_temp$time - phi))) * out_temp$H2 / out_temp$i_ILI
          
        } else {
          
          rho1_w <- rho1 * out_temp$H1 / out_temp$i_ILI
          rho2_w <- rho2 * out_temp$H2 / out_temp$i_ILI
          
        }
        
        rho1_w[rho1_w > 1.0 & !is.na(rho1_w)] <- 1.0
        rho2_w[rho2_w > 1.0 & !is.na(rho2_w)] <- 1.0
        
        expect_equal(nrow(out_temp), length(rho1_w))
        expect_equal(nrow(out_temp), length(rho2_w))
        
        out_temp$rho1_w <- rho1_w
        out_temp$rho2_w <- rho2_w
        
        if ('n_T1' %in% rownames(pomp_object@covar@table)) {
          
          out_temp <- out_temp %>%
            mutate(obs1 = rho1_w * n_T1,
                   obs2 = rho2_w * n_T2)
          
        } else {
          
          out_temp <- out_temp %>%
            mutate(obs1 = rho1_w * n_T,
                   obs2 = rho2_w * n_T)
          
        }
        
        if (analysis != 'iqr') {
          
          out_temp <- out_temp %>%
            select(time:H2, obs1:obs2)
          
        }
        
      }
      
      if (return_data) {
        
        out_temp <- out_temp %>%
          cbind(t(pomp_object@data))
        
      }
      
    } else if (model_type == 'stochastic') {
      
      out_temp <- simulate(pomp_object, nsim = n_sim, format = 'data.frame') %>%
        mutate(season = seas) %>%
        select(season, time:.id, n_P1:n_P2) %>%
        as_tibble()
      
      if (return_data) {
        
        out_temp <- out_temp %>%
          arrange(.id, time) %>%
          cbind(t(pomp_object@data))
        names(out_temp)[6:7] <- c('obs1', 'obs2')
        
        out_temp <- out_temp %>%
          as_tibble() %>%
          arrange(time, .id)
        
      }
      
    } else {
      
      stop('Model type neither deterministic nor stochastic!')
      
    }
    
  }
  
  # Return simulations:
  return(out_temp)
  
}


calculate_metrics <- function(dat) {
  # Function to calculate various outbreak metrics
  # params dat: tibble of the number of cases of both viruses over time
  # returns: Tibble containing the peak timing, peak intensity, and attack rates for both viruses
  
  if ('.id' %in% names(dat)) {
    
    out_metrics <- dat %>%
      group_by(.id) %>%
      summarise(pt1 = which.max(n_P1),
                pt2 = which.max(n_P2),
                pi1 = max(n_P1, na.rm = TRUE),
                pi2 = max(n_P2, na.rm = TRUE),
                ar1 = sum(n_P1, na.rm = TRUE),
                ar2 = sum(n_P2, na.rm = TRUE))
    
  } else {
    
    out_metrics <- dat %>%
      summarise(pt1 = which.max(n_P1),
                pt2 = which.max(n_P2),
                pi1 = max(n_P1, na.rm = TRUE),
                pi2 = max(n_P2, na.rm = TRUE),
                ar1 = sum(n_P1, na.rm = TRUE),
                ar2 = sum(n_P2, na.rm = TRUE))
    
  }
  
  return(out_metrics)
  
}


calculate_duration_and_concentration <- function(dat) {
  # Function to calculate outbreak duration (defined as number of weeks containing 75% of reported cases), and to
  # check whether these weeks are consecutive
  # params dat: tibble of the number of cases of both viruses over time for a single season
  # returns: Tibble containing duration and consecutive-ness of both viruses
  
  # List of results:
  res_list <- vector('list', length = 2)
  
  # Loop through viruses:
  for (i in 1:2) {
    
    vir <- c('n_P1', 'n_P2')[i]
    
    case_counts_temp <- dat %>% pull(vir)
    case_counts_temp[is.na(case_counts_temp)] <- 0
    
    target_sum <- sum(case_counts_temp)
    sum_cases <- 0
    which_weeks <- c()
    
    while(sum_cases < 0.75 * target_sum) {
      
      which_max <- which.max(case_counts_temp)
      max_val <- case_counts_temp[which_max]
      
      sum_cases <- sum_cases + max_val
      which_weeks <- c(which_weeks, which_max)
      
      case_counts_temp[which_max] <- 0
      
    }
    
    res_list[[i]] <- c(vir, length(which_weeks), min(which_weeks) + length(which_weeks) - 1 == max(which_weeks))
    
  }
  
  out <- as_tibble(t(c(dur1 = as.numeric(res_list[[1]][2]), dur2 = as.numeric(res_list[[2]][2]),
                       cons1 = as.numeric(as.logical(res_list[[1]][3])), cons2 = as.numeric(as.logical(res_list[[2]][3])))))
  return(out)
  
}


check_accuracy_metrics <- function(sim_metrics, obs_metrics, seas, pt_acc, pi_acc) {
  # Function to determine what proportion of simulations have peak timing/intensity and
  # attack rates within a given error of observed outbreaks
  # param sim_metrics: Tibble of metrics calculated from simulated outbreaks
  # param obs_metrics: Tibble of metrics calculated from observed outbreaks
  # param seas: Season for which simulations were performed
  # param pt_acc: Peak timing is considered accurate if within this many weeks from observed (numeric)
  # param pi_acc: Peak intensity and attack rates are considered accurate if within this percent of observed (numeric)
  # returns: Tibble containing the percent of simulation accurate within given errors
  
  # Process pi_acc:
  pi_acc <- 1 + (pi_acc / 100)
  
  # For each simulation, calculate whether metrics are within given range:
  sim_metrics <- sim_metrics %>%
    mutate(pt1 = (abs(pt1 - obs_metrics$pt1) <= pt_acc),
           pt2 = (abs(pt2 - obs_metrics$pt2) <= pt_acc),
           pi1 = (abs(pi1 - obs_metrics$pi1) + obs_metrics$pi1) <= pi_acc * obs_metrics$pi1,
           pi2 = (abs(pi2 - obs_metrics$pi2) + obs_metrics$pi2) <= pi_acc * obs_metrics$pi2,
           ar1 = (abs(ar1 - obs_metrics$ar1) + obs_metrics$ar1) <= pi_acc * obs_metrics$ar1,
           ar2 = (abs(ar2 - obs_metrics$ar2) + obs_metrics$ar2) <= pi_acc * obs_metrics$ar2) %>%
    select(-.id)
  
  # Determine % of simulations within the given range of observed values:
  sim_metrics <- sim_metrics %>%
    summarise(pt1 = length(pt1[pt1]),
              pt2 = length(pt2[pt2]),
              pi1 = length(pi1[pi1]),
              pi2 = length(pi2[pi2]),
              ar1 = length(ar1[ar1]),
              ar2 = length(ar2[ar2])) %>%
    mutate(season = seas)
  
  # Output results:
  return(sim_metrics)
  
}
