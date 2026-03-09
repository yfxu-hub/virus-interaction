# Quantifying Viral Interactions and Vaccination Impacts Among COVID-19, Influenza, and RSV in Post-Pandemic Hong Kong

Code to infer the strength and duration of interactions between COVID-19, flu and RSV

Directory Structure
-------------------
    * process_results
    * bootstrap
    * vaccination_simulation_study

Formatting data
---------------------------

Data from Hong Kong should be downloaded from the Centre for health Protection Department of health of the Government of the Hong Kong SAR (https://www.chp.gov.hk/en/statistics/submenu/641/index.html).
----------

* "resp_interaction_model.c": C code used to run the model, including parameter transformations, observation model, and deterministic skeleton
* "resp_interaction_model.R": Code to read in location- and season-specific data, create a pomp model, and run various model checks
* "functions_flu_RSV.R": Functions to create pomp objects used for running and fitting the model
* "test_code.R": Contains various functions used to check that model code is behaving as expected
* "setup_global_likelihood.R": Code to read in data from all seasons and load functions for evaluating the global log-likelihood

Fitting the interaction model to data
-------------------------------------

1.Initial parameter fitting (Round 1)
Run fit_traj_matching_round1.R to obtain initial estimates for all season‑specific parameters. Set sobol_size = 500, search_type = "broad", and for Hong Kong analyses set sens = "main".
Run process_results/01_check_missing_files_traj_matching_r1.R to identify any failed runs (by virus, season, or starting parameter set). The script will issue warnings for missing result files.
Run process_results/02_compile_results_traj_matching_r1.R to aggregate individual results into consolidated output files.

2.Full parameter fitting (Round 2)
Run fit_traj_matching_round2.R to fit both shared and season‑specific parameters and obtain maximum likelihood estimates. Use search_type = "round1_CIs", which_round = 1, sobol_size = 500, int_eff = "susc", and prof_lik = FALSE. Set sens according to the location being analyzed. The option run_parallel may be set to TRUE or FALSE depending on whether to run multiple starting parameter sets in parallel.
Run get_start_ranges_from_round2.R to generate the starting parameter ranges for the next fitting round.
If only one parameter set is statistically supported—i.e., its log‑likelihood lies within qchisq(0.95, df = number_of_parameters) / 2 of the MLE—perform another fitting round, as the MLE has likely not yet been reached.

3.Iterative refinement (subsequent rounds)
Re‑run fit_traj_matching_round2.R with search_type = "round2_CIs" and increment which_round to 2.
Run get_start_ranges_from_round2.R again to update parameter start ranges.
Repeat this procedure—incrementing which_round each time—until multiple parameter sets are statistically supported.

4.Final fitting round
Once multiple parameter sets are supported, run one final round of fit_traj_matching_round2.R with search_type = "round2_CIs".
Run get_start_ranges_from_round2.R, which will compute the final starting ranges for parametric bootstrapping and save the maximum likelihood estimates for all parameters.
The confidence intervals for the parameters are derived from all parameter estimates that pass the chi‑square likelihood ratio test.

5.Run parametric bootstrapping to get distributions for seasonal parameters.
In the bootstrap folder, sequentially run the R scripts beginning with bootstrap_01, bootstrap_02, and bootstrap_03 to generate the synthetic datasets and the corresponding distributions of the seasonal parameters.

Code to explore data/model fit
------------------------------

* "process_results/calculate_obs_and_sim_metrics.R:"
  * Code to calculate several outbreak metrics for influenza and RSV, including attack rates, week of peak activity, and duration of outbreaks, and to compare metrics calculated from the data to those calculated baed on simulations performed at the MLEs.
* "functions_evaluate_res.R":
  * Contains functions to obtain deterministic or stochastic simulations from the model at the MLE, and to calculate relevant outbreak metrics.

Simulation Study of Vaccine Impact
----------------------------------

Run "vaccination_simulation_study/run_vaccination_simulation_study.R" to get simulations of COVID-19 vaccine for all seasons, vaccine coverage levels, and vaccine timings.
Run "vaccination_simulation_study/run_vaccination_simulation_study_for_fluvaccine.R" to get simulations of flu vaccine for all seasons, vaccine coverage levels, and vaccine timings.

Data Sources
-------------------------------

* Hong Kong dataset:
  * Virologic data: https://www.chp.gov.hk/en/statistics/data/10/641/642/2274.html
  * ILI data: https://www.chp.gov.hk/en/statistics/data/10/26/44/292/7010.html
  * COVID-19 & Flu Express: https://www.chp.gov.hk/en/resources/29/100148.html

