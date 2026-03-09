# ---------------------------------------------------------------------------------------------------------------------
# Functions to assist with flu/RSV model
# ---------------------------------------------------------------------------------------------------------------------

create_SITRxSITR_mod <- function(dat, Ri1_max = 3.0, Ri2_max = 3.0, d2_max = 10.0, debug_bool = F, sens = 'main', loc = 'hk') {
  # Function to create pomp object 
  # param dat: ILI and virological data (data frame, first time point must be 1) 
  # param Ri1_max: upper bound of initial reproduction no of virus 1 (double, passed as global argument in the C script)
  # param Ri2_max: upper bound of initial reproduction no of virus  2 (double, passed as global argument in the C script)
  # param d2_max: upper bound of multiplicative difference between delta1 and delta2 (double, passed as global argument in the C script) 
  # param debug_bool: should debugging info be printed? (boolean)
  # param sens: main analysis or the name of a specific sensitivity analysis
  # param loc: where are the data from? ('hk' or 'canada')
  
  # Same number of tests for flu and RSV?
  if (loc == 'hk') {
    test_diff <- FALSE
  } else {
    test_diff <- TRUE
  }
  
  # Read model C code:
  mod_code <- readLines('/home/yfxu/virus_interference/code/resp_interaction_model.c')
  components_nm <- c('globs', 'toest', 'fromest', 'dmeas_orig', 'dmeas_testdiff', 'rmeas_orig', 'rmeas_testdiff', 'rinit', 'skel', 'rsim')
  components_l <- vector(mode = 'list', length = length(components_nm))
  names(components_l) <- components_nm
  
  for (nm in components_nm) {
    components_l[[nm]] <- mod_code[str_which(mod_code, paste0('start_', nm)):str_which(mod_code, paste0('end_', nm))] %>%
      str_flatten(collapse = '\n')
    
    if(nm == 'globs') {
      components_l[[nm]] <- paste(components_l[[nm]], 
                                  sprintf('static int debug = %d; \nstatic double Ri1_max = %f; \nstatic double Ri2_max = %f; \nstatic double d2_max = %f; \nstatic char sens[] = "%s"; \nstatic char loc[] = "%s";', 
                                          as.integer(debug_bool),
                                          as.numeric(Ri1_max),
                                          as.numeric(Ri2_max),
                                          as.numeric(d2_max),
                                          as.character(sens),
                                          as.character(loc)),
                                  sep = '\n')
    }
    
    components_l[[nm]] <- Csnippet(text = components_l[[nm]])
  }
  
  # Check that population size is the same at all timepoints:
  expect_true(length(unique(dat$pop)) == 1)
  
  # Create pomp object:
  if (test_diff) {
    po <- pomp(data = dat[, c('time', 'n_P1', 'n_P2')],
               times = 'time',
               t0 = 0,
               covar = covariate_table(dat[, c('time', 'i_ILI', 'n_T1', 'n_T2', 'temp', 'ah', 'h3_inc', 'rhino_inc')], times = 'time'),
               accumvars = c('H1_tot', 'H2_tot', 'H1', 'H2'),
               obsnames = c('n_P1', 'n_P2'),
               statenames = c('X_SS', 'X_IS', 'X_TS', 'X_RS', 
                              'X_SI', 'X_II', 'X_TI', 'X_RI', 
                              'X_ST', 'X_IT', 'X_TT', 'X_RT',
                              'X_SR', 'X_IR', 'X_TR', 'X_RR', 
                              'H1_tot', 'H2_tot', 
                              'H1', 'H2'),
               paramnames = c('Ri1', 'Ri2', # initial effective reproductive numbers
                              'gamma1', 'gamma2', # 1 / average infectious periods
                              # 'delta', # 1 / average refractory period (assume same duration for flu and RSV)
                              'delta1', 'd2', #'delta2', # 1 / average refractory periods; relative length of refractory period for RSV->flu
                              'theta_lambda1', 'theta_lambda2', # interaction effects on susceptibility to infection
                              'rho1', 'rho2', # probs. infection leads to ILI consultation
                              'alpha', 'phi', # amplitude and phase of seasonality of all-cause consultations
                              'theta_rho1', 'theta_rho2', # interaction effects on severity of infections
                              'eta_temp1', 'eta_temp2', # temperature forcing on virus 1 and 2
                              'eta_ah1', 'eta_ah2', # absolute humidity on virus 1 and 2
                              'b1', 'b2', 'phi1', 'phi2', # amplitudes and phase shifts for sinusoidal seasonality
                              'beta_h3', # effect of H3 on RSV force of infection
                              'beta_rhino', # effect of rhinovirus on influenza force of infection
                              'beta_sd1', 'beta_sd2', # extrademographic stochasticity (k-value) for virus 1 and 2
                              'N', # population size
                              'I10', 'I20', # props. infectious at outbreak start
                              'R10', 'R20', 'R120'), # props. recovered at outbreak start
               params = c(Ri1 = 1.5, Ri2 = 1.5,
                          gamma1 = 7 / gamma_vir1, gamma2 = 7 / gamma_vir2, # or 4 for flu?
                          # delta = 7 / 5,
                          delta1 = 7 / 5, d2 = 1.0, #delta2 = 7 / 5,
                          theta_lambda1 = 1.0, theta_lambda2 = 1.0,
                          rho1 = 0.5, rho2 = 0.15,
                          alpha = 0, phi = 0,
                          theta_rho1 = 1.0, theta_rho2 = 1.0,
                          eta_temp1 = 0, eta_temp2 = 0,
                          eta_ah1 = 0, eta_ah2 = 0,
                          b1 = 0, b2 = 0, phi1 = 0, phi2 = 0,
                          beta_h3 = 0, beta_rhino = 0,
                          beta_sd1 = 0, beta_sd2 = 0,
                          N = unique(dat$pop),
                          I10 = 1e-5, I20 = 1e-5,
                          R10 = 0, R20 = 0, R120 = 0),
               globals = components_l[['globs']],
               dmeasure = components_l[['dmeas_testdiff']],
               rmeasure = components_l[['rmeas_testdiff']],
               skeleton = vectorfield(components_l[['skel']]),
               rprocess = euler(step.fun = components_l[['rsim']], delta.t = 0.01),
               partrans = parameter_trans(toEst = components_l[['toest']], fromEst = components_l[['fromest']]),
               rinit = components_l[['rinit']]
    )
    
  } else {
    po <- pomp(data = dat[, c('time', 'n_P1', 'n_P2')],
               times = 'time',
               t0 = 0,
               covar = covariate_table(dat[, c('time', 'i_ILI', 'n_T', 'temp', 'ah', 'h3_inc', 'rhino_inc')], times = 'time'),
               accumvars = c('H1_tot', 'H2_tot', 'H1', 'H2'),
               obsnames = c('n_P1', 'n_P2'),
               statenames = c('X_SS', 'X_IS', 'X_TS', 'X_RS', 
                              'X_SI', 'X_II', 'X_TI', 'X_RI', 
                              'X_ST', 'X_IT', 'X_TT', 'X_RT',
                              'X_SR', 'X_IR', 'X_TR', 'X_RR', 
                              'H1_tot', 'H2_tot', 
                              'H1', 'H2'),
               paramnames = c('Ri1', 'Ri2', # initial effective reproductive numbers
                              'gamma1', 'gamma2', # 1 / average infectious periods
                              # 'delta', # 1 / average refractory period (assume same duration for flu and RSV)
                              'delta1', 'd2', #'delta2', # 1 / average refractory periods; relative length of refractory period for RSV->flu
                              'theta_lambda1', 'theta_lambda2', # interaction effects on susceptibility to infection
                              'rho1', 'rho2', # probs. infection leads to ILI consultation
                              'alpha', 'phi', # amplitude and phase of seasonality of all-cause consultations
                              'theta_rho1', 'theta_rho2', # interaction effects on severity of infections
                              'eta_temp1', 'eta_temp2', # temperature forcing on virus 1 and 2
                              'eta_ah1', 'eta_ah2', # absolute humidity on virus 1 and 2
                              'b1', 'b2', 'phi1', 'phi2', # amplitudes and phase shifts for sinusoidal seasonality
                              'beta_h3', # effect of H3 on RSV force of infection
                              'beta_rhino', # effect of rhinovirus on influenza force of infection
                              'beta_sd1', 'beta_sd2', # extrademographic stochasticity (k-value) for virus 1 and 2
                              'N', # population size
                              'I10', 'I20', # props. infectious at outbreak start
                              'R10', 'R20', 'R120'), # props. recovered at outbreak start
               params = c(Ri1 = 1.5, Ri2 = 1.5,
                          gamma1 = 7 / gamma_vir1, gamma2 = 7 / gamma_vir2, # or 4 for flu?
                          # delta = 7 / 5,
                          delta1 = 7 / 5, d2 = 1.0, #delta2 = 7 / 5,
                          theta_lambda1 = 1.0, theta_lambda2 = 1.0,
                          rho1 = 0.5, rho2 = 0.15,
                          alpha = 0, phi = 0,
                          theta_rho1 = 1.0, theta_rho2 = 1.0,
                          eta_temp1 = 0, eta_temp2 = 0,
                          eta_ah1 = 0, eta_ah2 = 0,
                          b1 = 0, b2 = 0, phi1 = 0, phi2 = 0,
                          beta_h3 = 0, beta_rhino = 0,
                          beta_sd1 = 0, beta_sd2 = 0,
                          N = unique(dat$pop),
                          I10 = 1e-5, I20 = 1e-5,
                          R10 = 0, R20 = 0, R120 = 0),
               globals = components_l[['globs']],
               dmeasure = components_l[['dmeas_orig']],
               rmeasure = components_l[['rmeas_orig']],
               skeleton = vectorfield(components_l[['skel']]),
               rprocess = euler(step.fun = components_l[['rsim']], delta.t = 0.01),
               partrans = parameter_trans(toEst = components_l[['toest']], fromEst = components_l[['fromest']]),
               rinit = components_l[['rinit']]
    )
  }
  
  return(po)
}


create_SITRxSITR_mod_VACC <- function(dat, Ri1_max = 3.0, Ri2_max = 3.0, d2_max = 10.0, t0_eff = 0, debug_bool = F, sens = 'main', test_diff = FALSE) {
  # Function to create pomp object with influenza vaccination
  # param dat: ILI and virological data (data frame) 
  # param Ri1_max: upper bound of initial reproduction no of virus 1 (double, passed as global argument in the C script)
  # param Ri2_max: upper bound of initial reproduction no of virus  2 (double, passed as global argument in the C script)
  # param d2_max: upper bound of multiplicative difference between delta1 and delta2 (double, passed as global argument in the C script) 
  # param t0_eff: time at which the simulation should be started
  # param debug_bool: should debugging info be printed? (boolean)
  # param sens: main analysis or the name of a specific sensitivity analysis
  # param test_diff: are there a different number of tests performed for the two viruses being modeled?
  
  # Are data from Canada?:
  if (test_diff) {
    loc <- 'canada'
  } else {
    loc <- 'hk'
  }
  
  # Read model C code:
  mod_code <- readLines('/home/yfxu/virus_interference/code/vaccination_simulation_study/resp_interaction_model_VACC.c')
  components_nm <- c('globs', 'toest', 'fromest', 'dmeas_orig', 'dmeas_testdiff', 'rmeas_orig', 'rmeas_testdiff', 'rinit', 'skel')
  components_l <- vector(mode = 'list', length = length(components_nm))
  names(components_l) <- components_nm
  
  for (nm in components_nm) {
    components_l[[nm]] <- mod_code[str_which(mod_code, paste0('start_', nm)):str_which(mod_code, paste0('end_', nm))] %>%
      str_flatten(collapse = '\n')
    
    if(nm == 'globs') {
      components_l[[nm]] <- paste(components_l[[nm]], 
                                  sprintf('static int debug = %d; \nstatic double Ri1_max = %f; \nstatic double Ri2_max = %f; \nstatic double d2_max = %f; \nstatic char sens[] = "%s"; \nstatic char loc[] = "%s";', 
                                          as.integer(debug_bool),
                                          as.numeric(Ri1_max),
                                          as.numeric(Ri2_max),
                                          as.numeric(d2_max),
                                          as.character(sens),
                                          as.character(loc)), 
                                  sep = '\n')
    }
    
    components_l[[nm]] <- Csnippet(text = components_l[[nm]])
  }
  
  # Check that population size is the same at all timepoints:
  expect_true(length(unique(dat$pop)) == 1)
  
  # Create pomp object:
  if (test_diff) {
    
    po <- pomp(data = dat[, c('time', 'n_P1', 'n_P2')],
               times = 'time',
               t0 = t0_eff,
               covar = covariate_table(dat[, c('time', 'i_ILI', 'n_T1', 'n_T2', 'temp', 'ah')], times = 'time'),
               accumvars = c('H1', 'H2'),
               obsnames = c('n_P1', 'n_P2'),
               statenames = c('X_SS', 'X_IS', 'X_TS', 'X_RS', 
                              'X_SI', 'X_II', 'X_TI', 'X_RI', 
                              'X_ST', 'X_IT', 'X_TT', 'X_RT',
                              'X_SR', 'X_IR', 'X_TR', 'X_RR',
                              'V_SS1', 'V_SS2', 'V_SI', 'V_ST', 'V_SR', 'V_RS',
                              'H1', 'H2'),
               paramnames = c('Ri1', 'Ri2', # initial effective reproductive numbers
                              'gamma1', 'gamma2', # 1 / average infectious periods
                              # 'delta', # 1 / average refractory period (assume same duration for flu and RSV)
                              'delta1', 'd2', #'delta2', # 1 / average refractory periods; relative length of refractory period for RSV->flu
                              'theta_lambda1', 'theta_lambda2', # interaction effects on susceptibility to infection
                              'rho1', 'rho2', # probs. infection leads to ILI consultation
                              'alpha', 'phi', # amplitude and phase of seasonality of all-cause consultations
                              'theta_rho1', 'theta_rho2', # interaction effects on severity of infections
                              'eta_temp1', 'eta_temp2', # temperature forcing on virus 1 and 2
                              'eta_ah1', 'eta_ah2', # absolute humidity on virus 1 and 2
                              'b1', 'b2', 'phi1', 'phi2', # amplitudes and phase shifts for sinusoidal seasonality
                              'beta_sd1', 'beta_sd2', # extrademographic stochasticity (k-value) for virus 1 and 2
                              'N', # population size
                              'theta_lambda_vacc', # strength of indirect vaccine protection against virus 2
                              'delta_vacc', # 1 / average duration of indirect vaccine protection against virus 2
                              'vacc_eff', # vaccine efficacy against virus 1
                              'p_vacc', # prop. vaccination against virus 1
                              'I10', 'I20', # props. infectious at outbreak start
                              'R10', 'R20', 'R120', # props. recovered at outbreak start
                              'xss0', 'xis0', 'xts0', 'xrs0', 'xsi0', 'xii0', 'xti0', 'xri0',# NUMBER in each compartment at the time of vaccination
                              'xst0', 'xit0', 'xtt0', 'xrt0', 'xsr0', 'xir0', 'xtr0', 'xrr0'), # Note: Set to <0 to base initial conditions off proportions instead
               params = c(Ri1 = 2, Ri2 =1.5,
                          gamma1 = 7 / gamma_vir1, gamma2 = 7 / gamma_vir2, # or 4 for flu?
                          # delta = 7 / 5,
                          delta1 = 7 / 5, d2 = 1.0, #delta2 = 7 / 5,
                          theta_lambda1 = 1.0, theta_lambda2 = 1.0,
                          rho1 = 0.5, rho2 = 0.15,
                          alpha = 0, phi = 0,
                          theta_rho1 = 1.0, theta_rho2 = 1.0,
                          eta_temp1 = 0, eta_temp2 = 0,
                          eta_ah1 = 0, eta_ah2 = 0,
                          b1 = 0, b2 = 0, phi1 = 0, phi2 = 0,
                          beta_sd1 = 0, beta_sd2 = 0,
                          N = unique(dat$pop),
                          theta_lambda_vacc = 1.0,
                          delta_vacc = 7 / 5,
                          vacc_eff = 0,
                          p_vacc = 0,
                          I10 = 1e-5, I20 = 1e-5,
                          R10 = 0, R20 = 0, R120 = 0,
                          xss0 = -1, xis0 = -1, xts0 = -1, xrs0 = -1, xsi0 = -1, xii0 = -1, xti0 = -1, xri0 = -1,
                          xst0 = -1, xit0 = -1, xtt0 = -1, xrt0 = -1, xsr0 = -1, xir0 = -1, xtr0 = -1, xrr0 = -1),
               globals = components_l[['globs']],
               dmeasure = components_l[['dmeas_testdiff']],
               rmeasure = components_l[['rmeas_testdiff']],
               skeleton = vectorfield(components_l[['skel']]),
               partrans = parameter_trans(toEst = components_l[['toest']], fromEst = components_l[['fromest']]),
               rinit = components_l[['rinit']]
    )
    
  } else {
    
    po <- pomp(data = dat[, c('time', 'n_P1', 'n_P2')],
               times = 'time',
               t0 = t0_eff,
               covar = covariate_table(dat[, c('time', 'i_ILI', 'n_T', 'temp', 'ah')], times = 'time'),
               accumvars = c('H1', 'H2'),
               obsnames = c('n_P1', 'n_P2'),
               statenames = c('X_SS', 'X_IS', 'X_TS', 'X_RS', 
                              'X_SI', 'X_II', 'X_TI', 'X_RI', 
                              'X_ST', 'X_IT', 'X_TT', 'X_RT',
                              'X_SR', 'X_IR', 'X_TR', 'X_RR',
                              'V_SS1', 'V_SS2', 'V_SI', 'V_ST', 'V_SR', 'V_RS',
                              'H1', 'H2'),
               paramnames = c('Ri1', 'Ri2', # initial effective reproductive numbers
                              'gamma1', 'gamma2', # 1 / average infectious periods
                              # 'delta', # 1 / average refractory period (assume same duration for flu and RSV)
                              'delta1', 'd2', #'delta2', # 1 / average refractory periods; relative length of refractory period for RSV->flu
                              'theta_lambda1', 'theta_lambda2', # interaction effects on susceptibility to infection
                              'rho1', 'rho2', # probs. infection leads to ILI consultation
                              'alpha', 'phi', # amplitude and phase of seasonality of all-cause consultations
                              'theta_rho1', 'theta_rho2', # interaction effects on severity of infections
                              'eta_temp1', 'eta_temp2', # temperature forcing on virus 1 and 2
                              'eta_ah1', 'eta_ah2', # absolute humidity on virus 1 and 2
                              'b1', 'b2', 'phi1', 'phi2', # amplitudes and phase shifts for sinusoidal seasonality
                              'beta_sd1', 'beta_sd2', # extrademographic stochasticity (k-value) for virus 1 and 2
                              'N', # population size
                              'theta_lambda_vacc', # strength of indirect vaccine protection against virus 2
                              'delta_vacc', # 1 / average duration of indirect vaccine protection against virus 2
                              'vacc_eff', # vaccine efficacy against virus 1
                              'p_vacc', # prop. vaccination against virus 1
                              'I10', 'I20', # props. infectious at outbreak start
                              'R10', 'R20', 'R120', # props. recovered at outbreak start
                              'xss0', 'xis0', 'xts0', 'xrs0', 'xsi0', 'xii0', 'xti0', 'xri0',# NUMBER in each compartment at the time of vaccination
                              'xst0', 'xit0', 'xtt0', 'xrt0', 'xsr0', 'xir0', 'xtr0', 'xrr0'), # Note: Set to <0 to base initial conditions off proportions instead
               params = c(Ri1 = 2, Ri2 = 1.5,
                          gamma1 = 7 /gamma_vir1, gamma2 = 7 / gamma_vir2, # or 4 for flu?
                          # delta = 7 / 5,
                          delta1 = 7 / 5, d2 = 1.0, #delta2 = 7 / 5,
                          theta_lambda1 = 1.0, theta_lambda2 = 1.0,
                          rho1 = 0.5, rho2 = 0.15,
                          alpha = 0, phi = 0,
                          theta_rho1 = 1.0, theta_rho2 = 1.0,
                          eta_temp1 = 0, eta_temp2 = 0,
                          eta_ah1 = 0, eta_ah2 = 0,
                          b1 = 0, b2 = 0, phi1 = 0, phi2 = 0,
                          beta_sd1 = 0, beta_sd2 = 0,
                          N = unique(dat$pop),
                          theta_lambda_vacc = 1.0,
                          delta_vacc = 7 / 5,
                          vacc_eff = 0,
                          p_vacc = 0,
                          I10 = 1e-5, I20 = 1e-5,
                          R10 = 0, R20 = 0, R120 = 0,
                          xss0 = -1, xis0 = -1, xts0 = -1, xrs0 = -1, xsi0 = -1, xii0 = -1, xti0 = -1, xri0 = -1,
                          xst0 = -1, xit0 = -1, xtt0 = -1, xrt0 = -1, xsr0 = -1, xir0 = -1, xtr0 = -1, xrr0 = -1),
               globals = components_l[['globs']],
               dmeasure = components_l[['dmeas_orig']],
               rmeasure = components_l[['rmeas_orig']],
               skeleton = vectorfield(components_l[['skel']]),
               partrans = parameter_trans(toEst = components_l[['toest']], fromEst = components_l[['fromest']]),
               rinit = components_l[['rinit']]
    )
    
  }
  
  return(po)
}
