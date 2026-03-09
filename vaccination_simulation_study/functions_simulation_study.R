# ---------------------------------------------------------------------------------------------------------------------
# Functions to run flu/RSV model with instantaneous vaccination
# ---------------------------------------------------------------------------------------------------------------------

run_simulation_with_vaccination <- function(dat, t_vacc, mod_parms, Ri_max1, Ri_max2, d2_max, debug_bool, sens = 'main', test_diff = FALSE) {
  # Function to run deterministic simulations of the SITRxSITR model with instantaneous vaccination at any time point
  # param dat: Virological, ILI, and covariate data (tibble)
  # param t_vacc: The week at which to vaccinate some proportion of the susceptible population (numeric)
  # param mod_parms: Either a named vector or a matrix (created using parmat) of parameter values for model runs
  # param Ri1_max: Upper bound of initial reproduction no of virus 1 (double, passed as global argument in the C script)
  # param Ri2_max: Upper bound of initial reproduction no of virus  2 (double, passed as global argument in the C script)
  # param d2_max: Upper bound of multiplicative difference between delta1 and delta2 (double, passed as global argument in the C script) 
  # param debug_bool: Should debugging info be printed? (boolean)
  # param sens: main analysis or the name of a specific sensitivity analysis
  # param test_diff: are there a different number of tests performed for the two viruses being modeled?
  # returns: Data frame containing values of all states at all timepoints
  
  # # run_simulation_with_vaccination(dat_pomp, t_vacc, model_params, Ri_max1, Ri_max2, d2_max, debug_bool, sens, test_diff = test_diff_use)
  # dat=dat_pomp
  # t_vacc=5
  # mod_parms=model_params_ve0.5
  # sens = 'main'
  # test_diff = FALSE
  
  if (t_vacc == 0) {
    
    # Create pomp object to run simulation from beginning of season:
    resp_mod <- create_SITRxSITR_mod_VACC(dat = dat,
                                          Ri1_max = Ri_max1,
                                          Ri2_max = Ri_max2,
                                          d2_max = d2_max,
                                          t0_eff = 0,
                                          debug_bool = debug_bool,
                                          sens = sens,
                                          test_diff = test_diff)
    
    # Run model with input parameter values:
    if (is.null(ncol(mod_parms))) {
      
      # Set epidemiologic and vaccination parameters:
      resp_mod <- set_model_parameters(resp_mod, mod_parms, vaccinate = TRUE)
      if (debug_bool) print(coef(resp_mod))
      
      # Run model for full season:
      sim_determ <- trajectory(object = resp_mod, format = 'data.frame')
      
    } else {
      
      # Run model for full season:
      sim_determ <- trajectory(object = resp_mod, params = mod_parms, format = 'data.frame')
      
    }
    
  } else {
    
    # Create pomp object to run simulation from beginning of season:
    resp_mod <- create_SITRxSITR_mod_VACC(dat = dat,
                                          Ri1_max = Ri_max1,
                                          Ri2_max = Ri_max2,
                                          d2_max = d2_max,
                                          t0_eff = 0,
                                          debug_bool = debug_bool,
                                          sens = sens,
                                          test_diff = test_diff)
    
    # Run model with input parameter values:
    if (is.null(ncol(mod_parms))) {
      
      # Set epidemiologic parameters:
      resp_mod <- set_model_parameters(resp_mod, mod_parms, vaccinate = FALSE)
      if (debug_bool) print(coef(resp_mod))
      
      # Run model for full season:
      sim_determ_PRE <- trajectory(object = resp_mod, format = 'data.frame')
      
    } else {
      
      # Do not vaccinate yet:
      mod_parms_novacc <- mod_parms
      mod_parms_novacc['p_vacc', ] <- 0
      mod_parms_novacc['vacc_eff', ] <- 0
      
      # Run model for full season:
      sim_determ_PRE <- trajectory(object = resp_mod, params = mod_parms_novacc, format = 'data.frame')
      # sim_determ_PRE1 <- trajectory(object = resp_mod, params = mod_parms_novacc[,1], format = 'data.frame')
      # sim_determ_PRE2 <- trajectory(object = resp_mod, params = mod_parms_novacc[,2], format = 'data.frame')
      
      # sim_determ_PRE<-rbind(sim_determ_PRE1,sim_determ_PRE2)
      
     
    }
    
    # Get subset of data post-vaccination and create new pomp object:
    dat_vacc <- dat %>%
      filter(time >= t_vacc)
    resp_mod_vacc <- create_SITRxSITR_mod_VACC(dat = dat_vacc,
                                               Ri1_max = Ri_max1,
                                               Ri2_max = Ri_max2,
                                               d2_max = d2_max,
                                               t0_eff = t_vacc,
                                               debug_bool = debug_bool,
                                               sens = sens,
                                               test_diff = test_diff)
    
    # Set epidemiologic and vaccination parameters:
    if (is.null(ncol(mod_parms))) {
      
      # Set epidemiologic and vaccination parameters:
      resp_mod_vacc <- set_model_parameters(resp_mod_vacc, mod_parms, vaccinate = TRUE)
      if (debug_bool) print(coef(resp_mod_vacc))
      
      # Set initial conditions:
      coef(resp_mod_vacc)[c('xss0', 'xis0', 'xts0', 'xrs0', 'xsi0', 'xii0', 'xti0', 'xri0',
                            'xst0', 'xit0', 'xtt0', 'xrt0', 'xsr0', 'xir0', 'xtr0', 'xrr0')] <-
        c(sim_determ_PRE[t_vacc, 'X_SS'], sim_determ_PRE[t_vacc, 'X_IS'], sim_determ_PRE[t_vacc, 'X_TS'], sim_determ_PRE[t_vacc, 'X_RS'],
          sim_determ_PRE[t_vacc, 'X_SI'], sim_determ_PRE[t_vacc, 'X_II'], sim_determ_PRE[t_vacc, 'X_TI'], sim_determ_PRE[t_vacc, 'X_RI'],
          sim_determ_PRE[t_vacc, 'X_ST'], sim_determ_PRE[t_vacc, 'X_IT'], sim_determ_PRE[t_vacc, 'X_TT'], sim_determ_PRE[t_vacc, 'X_RT'],
          sim_determ_PRE[t_vacc, 'X_SR'], sim_determ_PRE[t_vacc, 'X_IR'], sim_determ_PRE[t_vacc, 'X_TR'], sim_determ_PRE[t_vacc, 'X_RR'])
      if (debug_bool) print(coef(resp_mod_vacc))
      
      # Run model from t_vacc onward:
      sim_determ_POST <- trajectory(object = resp_mod_vacc, format = 'data.frame')
      
    } else {
      
      # Set initial conditions:
      for (run in 1:ncol(mod_parms)) {
        
        mod_parms['xss0', run] <- sim_determ_PRE[sim_determ_PRE$.id == run & sim_determ_PRE$time == t_vacc, 'X_SS']
        mod_parms['xis0', run] <- sim_determ_PRE[sim_determ_PRE$.id == run & sim_determ_PRE$time == t_vacc, 'X_IS']
        mod_parms['xts0', run] <- sim_determ_PRE[sim_determ_PRE$.id == run & sim_determ_PRE$time == t_vacc, 'X_TS']
        mod_parms['xrs0', run] <- sim_determ_PRE[sim_determ_PRE$.id == run & sim_determ_PRE$time == t_vacc, 'X_RS']
        
        mod_parms['xsi0', run] <- sim_determ_PRE[sim_determ_PRE$.id == run & sim_determ_PRE$time == t_vacc, 'X_SI']
        mod_parms['xii0', run] <- sim_determ_PRE[sim_determ_PRE$.id == run & sim_determ_PRE$time == t_vacc, 'X_II']
        mod_parms['xti0', run] <- sim_determ_PRE[sim_determ_PRE$.id == run & sim_determ_PRE$time == t_vacc, 'X_TI']
        mod_parms['xri0', run] <- sim_determ_PRE[sim_determ_PRE$.id == run & sim_determ_PRE$time == t_vacc, 'X_RI']
        
        mod_parms['xst0', run] <- sim_determ_PRE[sim_determ_PRE$.id == run & sim_determ_PRE$time == t_vacc, 'X_ST']
        mod_parms['xit0', run] <- sim_determ_PRE[sim_determ_PRE$.id == run & sim_determ_PRE$time == t_vacc, 'X_IT']
        mod_parms['xtt0', run] <- sim_determ_PRE[sim_determ_PRE$.id == run & sim_determ_PRE$time == t_vacc, 'X_TT']
        mod_parms['xrt0', run] <- sim_determ_PRE[sim_determ_PRE$.id == run & sim_determ_PRE$time == t_vacc, 'X_RT']
        
        mod_parms['xsr0', run] <- sim_determ_PRE[sim_determ_PRE$.id == run & sim_determ_PRE$time == t_vacc, 'X_SR']
        mod_parms['xir0', run] <- sim_determ_PRE[sim_determ_PRE$.id == run & sim_determ_PRE$time == t_vacc, 'X_IR']
        mod_parms['xtr0', run] <- sim_determ_PRE[sim_determ_PRE$.id == run & sim_determ_PRE$time == t_vacc, 'X_TR']
        mod_parms['xrr0', run] <- sim_determ_PRE[sim_determ_PRE$.id == run & sim_determ_PRE$time == t_vacc, 'X_RR']
        
      }
      
      # Run model from t_vacc onward:
      sim_determ_POST <- trajectory(object = resp_mod_vacc, params = mod_parms, format = 'data.frame')
      
    }
    
    # Combine two tibbles into one outcome tibble:
    if (is.null(ncol(mod_parms))) {
      
      if (t_vacc == 1) {
        
        sim_determ_POST[1, c('H1', 'H2')] <- sim_determ_PRE[t_vacc, c('H1', 'H2')]
        sim_determ <- sim_determ_POST
        
      } else {
        
        sim_determ_POST[1, c('H1', 'H2')] <- sim_determ_PRE[t_vacc, c('H1', 'H2')]
        sim_determ <- bind_rows(sim_determ_PRE[1:(t_vacc - 1), ], sim_determ_POST)
        
      }
      
    } else {
      
      if (t_vacc == 1) {
        
        sim_determ_POST[sim_determ_POST$time == t_vacc, c('H1', 'H2')] <- sim_determ_PRE[sim_determ_PRE$time == t_vacc, c('H1', 'H2')]
        sim_determ <- sim_determ_POST
        
      } else {
        
        sim_determ_POST[sim_determ_POST$time == t_vacc, c('H1', 'H2')] <- sim_determ_PRE[sim_determ_PRE$time == t_vacc, c('H1', 'H2')]
        sim_determ <- bind_rows(sim_determ_PRE[sim_determ_PRE$time < t_vacc, ], sim_determ_POST)
        
        
      }
      
    }
    
    # sim_determ[which(sim_determ$.id == 1),]<-sim_determ_PRE[which(sim_determ_PRE$.id == 1),]
    
    # Ensure that all compartments, except those involved in vaccination, are the same at the end of simulation 1 and the beginning of simulation 2
    if (debug_bool) {
      all.equal(sim_determ_PRE[sim_determ_PRE$time == t_vacc, ],
                sim_determ_POST[sim_determ_POST$time == t_vacc, ])
      
      if (is.null(ncol(mod_parms))) {
        print(sim_determ_POST[1, ] - sim_determ_PRE[t_vacc, ])
      }
      
    }
    
  }
  
  return(sim_determ)
  
}


set_model_parameters <- function(po, mod_parms, vaccinate) {
  # Function to update the parameter values of a pomp object
  # param po: Pomp object whose parameters should be updated
  # param mod_parms: Named vector of desired parameter values
  # param vaccinate: Boolean determining whether vaccine-related parameters should be updated or not
  # returns: Pomp object with updated parameter values
  
  for (paramname in names(coef(po))) {
    
    if (paramname %in% names(mod_parms)) {
      
      if (paramname %in% c('theta_lambda_vacc', 'delta_vacc', 'vacc_eff', 'p_vacc')) {
        
        if (vaccinate) {
          coef(po, paramname) <- mod_parms[paramname]
        }else {coef(po, paramname) <- 0}
        
      } else {
        
        coef(po, paramname) <- mod_parms[paramname]
        
      }
      
    }
    
  }
  
  return(po)
  
}
