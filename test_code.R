# ---------------------------------------------------------------------------------------------------------------------
# Code to test model setup and results for expected behavior
# ---------------------------------------------------------------------------------------------------------------------

# Load libraries:
library(testthat)

# Functions:

check_transformations <- function(pomp_object) {
  # Function to check that parameters are correctly being transformed
  # params pomp_object: The pomp model object to be checked
  
  obj <- pomp_object
  x <- coef(obj, transform = TRUE)
  obj1 <- obj
  coef(obj1, transform = TRUE) <- x
  
  expect_true(all.equal(coef(obj), coef(obj1)),
              info = 'Parameters not correctly transformed')
}


check_params <- function(pomp_object) {
  # Function to check that setting parameter values works as expected
  # params pomp_object: The pomp model object to be checked
  
  select_parms <- sample(1:length(coef(pomp_object)), 3)
  coef(pomp_object, names(coef(pomp_object))[select_parms]) <- c(0.1, 0.2, 0.3)
  expect_true(all(coef(pomp_object)[select_parms] == c(0.1, 0.2, 0.3)))
}


check_correct_N_CONST <- function(pomp_object, true_n, n_sim = 10) {
  # Function to check that deterministic and stochastic simulations maintain correct N (when pop size constant)
  # params pomp_object: The pomp model object to be checked
  # params true_n: The expected population size
  # params n_sim: The number of stochastic simulations to run; if 0, don't run stochastic model
  
  sim_determ <- trajectory(object = pomp_object, format = 'data.frame') %>%
    rowwise() %>%
    mutate(Ncheck = sum(c_across(contains('X') | contains('V'))),
           Ntrue = true_n)
  expect_true(all.equal(sim_determ$Ncheck, sim_determ$Ntrue))
  
  if (n_sim > 0) {
    sim_stoch <- simulate(object = pomp_object, nsim = n_sim, format = 'data.frame') %>%
      rowwise() %>%
      mutate(Ncheck = sum(c_across(X_SS:X_RR)),
             Ntrue = true_n)
    expect_true(all.equal(sim_stoch$Ncheck, sim_stoch$Ntrue))
  }
  
}


check_correct_N_CONST_VACC <- function(sim_res, true_n) {
  # Function to check that deterministic simulations maintain correct N in the presence of vaccination (when pop size constant)
  # params sim_res: The deterministic simulation results from a single run of the model with vaccination
  # params true_n: The expected population size
  
  check_n <- sim_res %>%
    rowwise() %>%
    mutate(Ncheck = round(sum(c_across(contains('X') | contains('V'))), 5)) %>%
    pull(Ncheck) %>%
    unique()
  expect_true(length(check_n) == 1)
  expect_true(check_n == true_n)
  
}


check_obs_lessthan_samples <- function(pomp_object, n_sim = 10, test_diff = FALSE) {
  # Function to check that simulated case numbers never exceed the total number of tests performed
  # params pomp_object: The pomp model object to be checked
  # params n_sim: The number of stochastic simulations to run
  # param test_diff: are there a different number of tests performed for the two viruses being modeled?
  # returns: Plot of # of samples, observations, and simulated positive tests
  
  sim_stoch <- simulate(object = pomp_object, nsim = n_sim, format = 'data.frame')
  sim_stoch <- sim_stoch %>%
    as_tibble() %>%
    left_join(dat_pomp, by = 'time') %>%
    rename(sim_nP1 = n_P1.x,
           sim_nP2 = n_P2.x,
           obs_nP1 = n_P1.y,
           obs_nP2 = n_P2.y) %>%
    filter(!is.na(sim_nP1))
  
  if (test_diff) {
    
    expect_true(all(sim_stoch$sim_nP1 <= sim_stoch$n_T1))
    expect_true(all((sim_stoch$sim_nP2 <= sim_stoch$n_T2)[!is.na(sim_stoch$sim_nP2 <= sim_stoch$n_T2)]))
    
    sim_stoch <- sim_stoch %>%
      select(time:.id, sim_nP1:sim_nP2, obs_nP1:obs_nP2, n_T1:n_T2) %>%
      pivot_longer(sim_nP1:sim_nP2, names_to = 'sim_names', values_to = 'sim_vals') %>%
      pivot_longer(obs_nP1:obs_nP2, names_to = 'obs_names', values_to = 'obs_vals') %>%
      mutate(sim_names = str_sub(sim_names, 5, 7),
             obs_names = str_sub(obs_names, 5, 7)) %>%
      filter(sim_names == obs_names)
    
    p_temp <- ggplot(data = sim_stoch) + geom_line(aes(x = time, y = n_T1), lwd = 1.5) +
      geom_line(aes(x = time, y = n_T2), lwd = 1.5, linetype = 2) +
      geom_line(aes(x = time, y = sim_vals, group = .id, col = .id)) +
      geom_point(aes(x = time, y = obs_vals)) +
      facet_wrap(~sim_names, ncol = 1) +
      labs(x = 'Time (Weeks)', y = 'Simulated Observations') + theme_classic() +
      scale_color_viridis(discrete = TRUE)
    
  } else {
    
    expect_true(all(sim_stoch$sim_nP1 <= sim_stoch$n_T))
    expect_true(all(sim_stoch$sim_nP2 <= sim_stoch$n_T))
    
    sim_stoch <- sim_stoch %>%
      select(time:.id, sim_nP1:sim_nP2, obs_nP1:obs_nP2, n_T) %>%
      pivot_longer(sim_nP1:sim_nP2, names_to = 'sim_names', values_to = 'sim_vals') %>%
      pivot_longer(obs_nP1:obs_nP2, names_to = 'obs_names', values_to = 'obs_vals') %>%
      mutate(sim_names = str_sub(sim_names, 5, 7),
             obs_names = str_sub(obs_names, 5, 7)) %>%
      filter(sim_names == obs_names)
    
    p_temp <- ggplot(data = sim_stoch) + geom_line(aes(x = time, y = n_T), lwd = 1.5) +
      geom_line(aes(x = time, y = sim_vals, group = .id, col = .id)) +
      geom_point(aes(x = time, y = obs_vals)) +
      facet_wrap(~sim_names, ncol = 1) +
      labs(x = 'Time (Weeks)', y = 'Simulated Observations') + theme_classic() +
      scale_color_viridis(discrete = TRUE)
    
  }
  
  return(p_temp)
}


check_independent_dynamics <- function(pomp_object) {
  # Function to check that, when no interaction is specified, virus dynamics are independent
  # params pomp_object: The pomp model object to be checked
  # returns: Plot of epidemic dynamics with and without interaction
  
  p_mat <- parmat(params = coef(pomp_object), nrep = 3)
  p_mat['I10', ] <- c(1e-5, 1e-5, 0)
  p_mat['I20', ] <- c(1e-5, 0, 1e-5)
  
  sim_determ <- trajectory(object = pomp_object,
                           params = p_mat,
                           format = 'data.frame') %>%
    pivot_longer(cols = -c('time', '.id'), names_to = 'var_nm', values_to = 'val') %>%
    mutate(var_type = if_else(str_detect(var_nm, 'X_'), 'state', 'accum')) %>%
    filter(var_nm %in% c('H1', 'H2'))
  
  expect_true(all.equal(sim_determ$val[sim_determ$.id == 1 & sim_determ$var_nm == 'H1'],
                        sim_determ$val[sim_determ$.id == 2 & sim_determ$var_nm == 'H1']))
  expect_true(all.equal(sim_determ$val[sim_determ$.id == 1 & sim_determ$var_nm == 'H2'],
                        sim_determ$val[sim_determ$.id == 3 & sim_determ$var_nm == 'H2']))
  
  p_temp <- ggplot(data = sim_determ %>% filter(.id == 1), aes(x = time, y = 100 * val)) +
    geom_point(aes(color = var_nm)) +
    geom_line(data = sim_determ %>% filter(.id == 2, var_nm == 'H1'), color = 'pink') +
    geom_line(data = sim_determ %>% filter(.id == 3, var_nm == 'H2'), color = 'purple') +
    labs(x = 'Time (Weeks)', y = 'Incidence (%)') +
    theme_classic()
  
  return(p_temp)
}


quick_explore_interaction <- function(pomp_object, int_vals, n_sim = 10) {
  # Function to get a quick idea of how the interaction strength impacts dynamics
  # params pomp_object: The pomp model object to use for simulations
  # params int_vals: A vector containing values of interaction strength to test
  # params n_sim: The number of stochastic simulations to run
  # returns: List of plots showing the effect of interactions of different strengths
  
  int_parms <- c('theta_lambda1', 'theta_lambda2', 'theta_rho1', 'theta_rho2')
  p_list <- vector('list', length(int_parms))
  
  for (i in 1:length(int_parms)) {
    
    p_mat <- parmat(params = coef(pomp_object), nrep = length(int_vals))
    p_mat[int_parms[i], ] <- int_vals
    
    sim_stoch <- simulate(object = pomp_object,
                          params = p_mat,
                          nsim = n_sim,
                          format = 'data.frame') %>%
      select(time:.id, n_P1:n_P2) %>%
      # select(time:.id, H1_tot:H2_tot) %>%
      pivot_longer(cols = -c('time', '.id'), names_to = 'vir', values_to = 'val') %>%
      mutate(id.parm = as.numeric(str_split(.id, '_') %>% map_chr(., 1)),
             id.sim = as.numeric(str_split(.id, '_') %>% map_chr(., 2))) %>%
      mutate(int_parm_val = int_vals[id.parm]) %>%
      mutate(int_parm_val = as.character(int_parm_val))
    
    p_temp <- ggplot(data = sim_stoch, aes(x = time, y = val, col = int_parm_val, linetype = vir, group = paste(.id, vir))) +
      geom_line() + theme_classic() +
      labs(x = 'Time', y = 'Cases', linetype = 'Virus', col = 'Param. Value', title = int_parms[i]) +
      scale_color_viridis(discrete = TRUE)
    
    p_list[[i]] <- p_temp
    
  }
  
  return(p_list)
}


check_independent_dynamics_VACC <- function(dat, t_vacc, mod_parms, Ri_max1, Ri_max2, d2_max, debug_bool) {
  # Function to check that, when the vaccine has no interaction effect, virus dynamics are independent
  # param dat: Virological, ILI, and covariate data (tibble)
  # param t_vacc: The week at which to vaccinate some proportion of the susceptible population (numeric)
  # param mod_parms: Either a named vector or a matrix (created using parmat) of parameter values for model runs
  # param Ri1_max: Upper bound of initial reproduction no of virus 1 (double, passed as global argument in the C script)
  # param Ri2_max: Upper bound of initial reproduction no of virus  2 (double, passed as global argument in the C script)
  # param d2_max: Upper bound of multiplicative difference between delta1 and delta2 (double, passed as global argument in the C script) 
  # param debug_bool: Should debugging info be printed? (boolean)
  # returns: Plot of flu and RSV when modeled together and in isolation
  
  sim_determ <- run_simulation_with_vaccination(dat, t_vacc, mod_parms, Ri_max1, Ri_max2, d2_max, debug_bool) %>%
    pivot_longer(cols = -c('time', '.id'), names_to = 'var_nm', values_to = 'val') %>%
    mutate(var_type = if_else(str_detect(var_nm, 'X_'), 'state', 'accum')) %>%
    filter(var_nm %in% c('H1', 'H2'))
  
  expect_true(all.equal(sim_determ$val[sim_determ$.id == 1 & sim_determ$var_nm == 'H1'],
                        sim_determ$val[sim_determ$.id == 2 & sim_determ$var_nm == 'H1']))
  expect_true(all.equal(sim_determ$val[sim_determ$.id == 1 & sim_determ$var_nm == 'H2'],
                        sim_determ$val[sim_determ$.id == 3 & sim_determ$var_nm == 'H2']))
  
  # print(all.equal(sim_determ$val[sim_determ$.id == 1 & sim_determ$var_nm == 'H1'],
  #                 sim_determ$val[sim_determ$.id == 2 & sim_determ$var_nm == 'H1']))
  # print(all.equal(sim_determ$val[sim_determ$.id == 1 & sim_determ$var_nm == 'H2'],
  #                 sim_determ$val[sim_determ$.id == 3 & sim_determ$var_nm == 'H2']))
  
  p_temp <- ggplot(data = sim_determ %>% filter(.id == 1), aes(x = time, y = 100 * val)) +
    geom_point(aes(color = var_nm)) +
    geom_line(data = sim_determ %>% filter(.id == 2, var_nm == 'H1'), color = 'pink') +
    geom_line(data = sim_determ %>% filter(.id == 3, var_nm == 'H2'), color = 'purple') +
    labs(x = 'Time (Weeks)', y = 'Incidence (%)') +
    theme_classic()
  
  return(p_temp)
}


check_single_virus_impact <- function(dat, t_vacc, mod_parms, Ri_max1, Ri_max2, d2_max, debug_bool, sens = 'main', test_diff = FALSE) {
  # Function to check whether a vaccine impacting only flu/RSV actually only impacts the outbreak of flu/RSV
  # param dat: Virological, ILI, and covariate data (tibble)
  # param t_vacc: The week at which to vaccinate some proportion of the susceptible population (numeric)
  # param mod_parms: Either a named vector or a matrix (created using parmat) of parameter values for model runs
  # param Ri1_max: Upper bound of initial reproduction no of virus 1 (double, passed as global argument in the C script)
  # param Ri2_max: Upper bound of initial reproduction no of virus  2 (double, passed as global argument in the C script)
  # param d2_max: Upper bound of multiplicative difference between delta1 and delta2 (double, passed as global argument in the C script) 
  # param debug_bool: Should debugging info be printed? (boolean)
  # param sens: main analysis or the name of a specific sensitivity analysis
  # param test_diff: are there a different number of tests performed for the two viruses being modeled?
  # returns: Plot of flu and RSV in populations with and without vaccination
  ##check_single_virus_impact(dat_pomp, t_vacc, model_params, Ri_max1, Ri_max2, d2_max, debug_bool, sens, test_diff = test_diff_use)
  # dat=dat_pomp
  # mod_parms=model_params
  # t_vacc=10 
  # sens = 'main'
  # test_diff = FALSE
  
  sim_determ <- run_simulation_with_vaccination(dat, t_vacc, mod_parms, Ri_max1, Ri_max2, d2_max, debug_bool, sens, test_diff) %>%
    pivot_longer(cols = -c('time', '.id'), names_to = 'var_nm', values_to = 'val') %>%
    mutate(var_type = if_else(str_detect(var_nm, 'X_'), 'state', 'accum')) %>%
    filter(var_nm %in% c('H1', 'H2'))
  
  if (unique(mod_parms['vacc_eff', ]) > 0) {
    
    expect_true(!all(sim_determ$val[sim_determ$var_nm == 'H1' & sim_determ$.id == 1] ==
                       sim_determ$val[sim_determ$var_nm == 'H1' & sim_determ$.id == 2]))
    expect_true(all.equal(sim_determ$val[sim_determ$var_nm == 'H2' & sim_determ$.id == 1],
                          sim_determ$val[sim_determ$var_nm == 'H2' & sim_determ$.id == 2]))
    
  } else {
    
    expect_true(all.equal(sim_determ$val[sim_determ$var_nm == 'H1' & sim_determ$.id == 1],
                          sim_determ$val[sim_determ$var_nm == 'H1' & sim_determ$.id == 2]))
    expect_true(!all(sim_determ$val[sim_determ$var_nm == 'H2' & sim_determ$.id == 1] ==
                       sim_determ$val[sim_determ$var_nm == 'H2' & sim_determ$.id == 2]))
    
  }
  
  p_temp <- ggplot(data = sim_determ %>% filter(.id == 1), aes(x = time, y = 100 * val, color = var_nm)) +
    geom_point() + geom_line(data = sim_determ %>% filter(.id == 2)) +
    labs(x = 'Time (Weeks)', y = 'Incidence (%)') +
    theme_classic()
  
  return(p_temp)
}
