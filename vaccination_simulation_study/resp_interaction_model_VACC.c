// Implementation of the SIR model for circulation of respiratory viruses

//start_globs
// Expit transform for parameters constrained in interval [a,b]
static double expitCons(double x, double a, double b) {
  double out = (a + b * exp(x)) / (1.0 + exp(x));
  if(ISNAN(out) | isinf(out)) out = (b + a * exp(-x)) / (1.0 + exp(-x)); // If x=+Inf, must return b
  return out;
}

// Logit transform for parameters constrained in interval [a,b]
static double logitCons(double x, double a, double b) {
  x = (x <= a) ? a : (x >= b ? b : x);
  double out = log((x - a) / (b - x));
  return out;
}

// Probability of transition for event of rate r during time step delta_t
// p = 1 - exp(-r * delta_t)
static double pTrans(double r, double delta_t) {
  
  // r: event (instantaneous rate)
  // delta_t: time step
  // Returns: probability of transition
  double out = 1.0 - exp(-r * delta_t);
  return out;
}
//end_globs

//start_toest
double sum_init = 0.0;

T_Ri1 = logitCons(Ri1, 1.0, Ri1_max); 
T_Ri2 = logitCons(Ri2, 1.0, Ri2_max); 
T_gamma1 = log(gamma1); 
T_gamma2 = log(gamma2);
//T_delta = log(delta);
T_delta1 = log(delta1);
//T_delta2 = log(delta2);
T_d2 = log(d2);
//T_d2 = logitCons(d2, 0.0, d2_max);
T_theta_lambda1 = logit(theta_lambda1);
T_theta_lambda2 = logit(theta_lambda2);
T_theta_rho1 = logit(theta_rho1);
T_theta_rho2 = logit(theta_rho2);

if (strcmp("canada", loc) == 0) {
  T_rho1 = log(rho1);
} else {
  T_rho1 = logit(rho1);
}
//T_rho1 = logit(rho1);

T_rho2 = logit(rho2);
//T_rho1 = log(rho1); 
//T_rho2 = log(rho2);
T_alpha = logit(alpha);
T_phi = logitCons(phi, 0.0, 52.25);
T_eta_temp1 = eta_temp1;
T_eta_temp2 = eta_temp2;
T_eta_ah1 = eta_ah1;
T_eta_ah2 = eta_ah2;
T_N = N;
T_beta_sd1 = log(beta_sd1);
T_beta_sd2 = log(beta_sd2);

T_theta_lambda_vacc = logit(theta_lambda_vacc);
T_delta_vacc = log(delta_vacc);
T_vacc_eff = logit(vacc_eff);
T_p_vacc = logit(p_vacc);

sum_init = I10 + I20 + R10 + R20 + R120;
//sum_init = I10 + R10 + R20 + R120;
T_I10 = log(I10 / (1.0 - sum_init));
T_I20 = log(I20 / (1.0 - sum_init));
T_R10 = log(R10 / (1.0 - sum_init));
T_R20 = log(R20 / (1.0 - sum_init));
T_R120 = log(R120 / (1.0 - sum_init));

//T_I20 = logit(I20);

//end_toest

//start_fromest
double sum_init = 0.0; 

Ri1 = expitCons(T_Ri1, 1.0, Ri1_max);
Ri2 = expitCons(T_Ri2, 1.0, Ri2_max);
gamma1 = exp(T_gamma1);
gamma2 = exp(T_gamma2);
//delta = exp(T_delta);
delta1 = exp(T_delta1);
//delta2 = exp(T_delta2);
d2 = exp(T_d2);
//d2 = expitCons(T_d2, 0.0, d2_max);
theta_lambda1 = expit(T_theta_lambda1);
theta_lambda2 = expit(T_theta_lambda2);
theta_rho1 = expit(T_theta_rho1);
theta_rho2 = expit(T_theta_rho2);

if (strcmp("canada", loc) == 0) {
  rho1 = exp(T_rho1);
} else {
  rho1 = expit(T_rho1);
}
//rho1 = expit(T_rho1);

rho2 = expit(T_rho2);
//rho1 = exp(T_rho1);
//rho2 = exp(T_rho2);
alpha = expit(T_alpha);
phi = expitCons(T_phi, 0.0, 52.25);
eta_temp1 = T_eta_temp1;
eta_temp2 = T_eta_temp2;
eta_ah1 = T_eta_ah1;
eta_ah2 = T_eta_ah2;
N = T_N;
beta_sd1 = exp(T_beta_sd1);
beta_sd2 = exp(T_beta_sd2);

theta_lambda_vacc = expit(T_theta_lambda_vacc);
delta_vacc = exp(T_delta_vacc);
vacc_eff = expit(T_vacc_eff);
p_vacc = expit(T_p_vacc);

sum_init = exp(T_I10) + exp(T_I20) + exp(T_R10) + exp(T_R20) + exp(T_R120);
//sum_init = exp(T_I10) + exp(T_R10) + exp(T_R20) + exp(T_R120);
I10 = exp(T_I10) / (1.0 + sum_init);
I20 = exp(T_I20) / (1.0 + sum_init);
R10 = exp(T_R10) / (1.0 + sum_init);
R20 = exp(T_R20) / (1.0 + sum_init);
R120 = exp(T_R120) / (1.0 + sum_init);

//I20 = expit(T_I20);

//end_fromest

//start_dmeas_orig
double fP1, fP2, ll;
double omega = (2 * M_PI) / 52.25;

double rho1_w = fmin2(1.0, rho1 * (1.0 + alpha * cos(omega * (t - phi))) * H1 / i_ILI); // Probability of detecting virus 1
double rho2_w = fmin2(1.0, rho2 * (1.0 + alpha * cos(omega * (t - phi))) * H2 / i_ILI); // Probability of detecting virus 2

if (rho1_w < 0) {
  rho1_w = 0.0;
}
if (rho2_w < 0) {
  rho2_w = 0.0;
}

fP1 = dbinom(n_P1, n_T, rho1_w, 1); // First likelihood component, natural scale
fP2 = dbinom(n_P2, n_T, rho2_w, 1); // Second likelihood component, natural scale

//fP1 = ISNA(n_P1) ? 0.0 : dbinom(n_P1, n_T, rho1_w, 1);
//fP2 = ISNA(n_P2) ? 0.0 : dbinom(n_P2, n_T, rho2_w, 1);
//if (isnan(fP1)) {
//  fP1 = dbinom(n_P1, n_T, 0.0, 1);
//}
//if (isnan(fP2)) {
//  fP2 = dbinom(n_P2, n_T, 0.0, 1);
//}

// If rho_w == 1, the resulting observation probability might be 0 (-Inf on log-scale)
// Replace by a big, but finite penalty if that's the case 
ll = fmax2(fP1 + fP2, -1e3);

// If data are NA, ll will be NA; in this case, set to zero
ll = ISNA(ll) ? 0.0 : ll;

if(debug) {
  Rprintf("t=%.1f, rho1_w=%.1f, rho2_w=%.1f, n_T=%.1f, fP1=%.1f, fP2=%.1f, sum=%.1f, ll=%.f\n", t, rho1_w, rho2_w, n_T, fP1, fP2, fP1 + fP2, ll);
} 
lik = (give_log) ? ll : exp(ll);

//end_dmeas_orig

//start_dmeas_testdiff
double fP1, fP2, ll;
double omega = (2 * M_PI) / 52.25;

double rho1_w = fmin2(1.0, rho1 * (1.0 + alpha * cos(omega * (t - phi))) * H1 / i_ILI); // Probability of detecting virus 1
double rho2_w = fmin2(1.0, rho2 * (1.0 + alpha * cos(omega * (t - phi))) * H2 / i_ILI); // Probability of detecting virus 2

if (rho1_w < 0) {
  rho1_w = 0.0;
}
if (rho2_w < 0) {
  rho2_w = 0.0;
}

fP1 = dbinom(n_P1, n_T1, rho1_w, 1); // First likelihood component, natural scale
fP2 = dbinom(n_P2, n_T2, rho2_w, 1); // Second likelihood component, natural scale

// If rho_w == 1, the resulting observation probability might be 0 (-Inf on log-scale)
// Replace by a big, but finite penalty if that's the case 
ll = fmax2(fP1 + fP2, -1e3);

// If data are NA, ll will be NA; in this case, set to zero
ll = ISNA(ll) ? 0.0 : ll;

if(debug) {
  Rprintf("t=%.1f, rho1_w=%.1f, rho2_w=%.1f, n_T1=%.1f, n_T2=%.1f, fP1=%.1f, fP2=%.1f, sum=%.1f, ll=%.f\n", t, rho1_w, rho2_w, n_T1, n_T2, fP1, fP2, fP1 + fP2, ll);
}
lik = (give_log) ? ll : exp(ll);

//end_dmeas_testdiff

//start_rmeas_orig
double omega = (2 * M_PI) / 52.25;

double rho1_w = fmin2(1.0, rho1 * (1.0 + alpha * cos(omega * (t - phi))) * H1 / i_ILI); // Probability of detecting virus 1
double rho2_w = fmin2(1.0, rho2 * (1.0 + alpha * cos(omega * (t - phi))) * H2 / i_ILI); // Probability of detecting virus 2

n_P1 = rbinom(n_T, rho1_w); // Generate of tests positive to virus 1
n_P2 = rbinom(n_T, rho2_w); // Generate of tests positive to virus 2
//end_rmeas_orig

//start_rmeas_testdiff
double omega = (2 * M_PI) / 52.25;

double rho1_w = fmin2(1.0, rho1 * (1.0 + alpha * cos(omega * (t - phi))) * H1 / i_ILI); // Probability of detecting virus 1
double rho2_w = fmin2(1.0, rho2 * (1.0 + alpha * cos(omega * (t - phi))) * H2 / i_ILI); // Probability of detecting virus 2

n_P1 = rbinom(n_T1, rho1_w); // Generate of tests positive to virus 1
n_P2 = rbinom(n_T2, rho2_w); // Generate of tests positive to virus 2
//end_rmeas_testdiff

//start_rinit
double X_SS_offset;

if (xss0 < 0) {
  
  X_SS = nearbyint((1.0 - I10 - I20 - R10 - R20 - R120) * N);
  X_IS = nearbyint(I10 * N);
  X_TS = 0;
  X_RS = nearbyint(R10 * N);
  X_SI = nearbyint(I20 * N);
  X_II = 0;
  X_TI = 0;
  X_RI = 0;
  X_ST = nearbyint((1.0 - p_vacc) * 0 * N);
  X_IT = 0;
  X_TT = 0;
  X_RT = 0;
  X_SR = nearbyint(R20 * N);
  X_IR = 0;
  X_TR = 0;
  X_RR = nearbyint(R120 * N);
  
  V_SS1 = nearbyint(p_vacc * X_SS);
  V_SR = nearbyint(p_vacc * X_SR);
  V_RS = nearbyint(p_vacc * X_RS);
  V_SI = nearbyint(p_vacc * X_SI);
  V_ST = 0;
  
  X_SS = X_SS - V_SS1;
  X_SR = X_SR - V_SR;
  X_RS = X_RS - V_RS;
  X_SI = X_SI - V_SI;
  
  X_SS_offset = X_IS + X_RS + X_SI + X_SR + X_RR + V_SS1 + V_SR + V_RS + V_SI;
  if ((X_SS + X_SS_offset) != N) {
    //Rprintf("X_SS=%.11f, sum=%.11f, N=%.1f\n", X_SS, X_SS + X_SS_offset, N);
    X_SS = nearbyint(N - X_SS_offset);
  }
  
} else {
  
  X_SS = (1.0 - p_vacc) * xss0;
  X_IS = xis0;
  X_TS = xts0;
  X_RS = (1.0 - p_vacc) * xrs0;
  X_SI = (1.0 - p_vacc) * xsi0;
  X_II = xii0;
  X_TI = xti0;
  X_RI = xri0;
  X_ST = (1.0 - p_vacc) * xst0;
  X_IT = xit0;
  X_TT = xtt0;
  X_RT = xrt0;
  X_SR = (1.0 - p_vacc) * xsr0;
  X_IR = xir0;
  X_TR = xtr0;
  X_RR = xrr0;
  
  V_SS1 = p_vacc * xss0;
  V_SR = p_vacc * xsr0;
  V_RS = p_vacc * xrs0;
  V_SI = p_vacc * xsi0;
  V_ST = p_vacc * xst0;
  
  X_SS_offset = X_IS + X_TS + X_RS + X_SI + X_II + X_TI + X_RI + X_ST + X_IT + X_TT + X_RT + X_SR + X_IR + X_TR + X_RR + V_SS1 + V_SR + V_RS + V_SI + V_ST;
  if ((X_SS + X_SS_offset) != N) {
    //Rprintf("X_SS=%.11f, sum=%.11f, N=%.1f\n", X_SS, X_SS + X_SS_offset, N);
    //X_SS = N - round((X_IS + X_TS + X_RS + X_SI + X_II + X_TI + X_RI + X_ST + X_IT + X_TT + X_RT + X_SR + X_IR + X_TR + X_RR + V_SS1 + V_SR + V_RS) * 100000000) / 100000000;
    X_SS = N - X_SS_offset;
    //Rprintf("X_SS=%.11f, sum=%.11f, N=%.1f\n", X_SS, X_SS + X_SS_offset, N);
  }
  
}

H1 = 0; 
H2 = 0;

V_SS2 = 0;
//V_SI = 0;
//V_ST = 0;

if ((X_SS + X_SS_offset) != N) {
  Rprintf("sum=%.10f, N=%.1f\n", X_SS + X_SS_offset, N);
}

if(debug) {
  Rprintf("%f, %f, %f, %f, %f, %f, %f\n", Ri1, Ri2, I10, I20, R10, R20, R120);
  Rprintf("%f, %f, %f, %f, %f, %f\n", X_SS, X_SR, X_RS, X_SI, X_ST, N);
  Rprintf("%f, %f, %f, %f, %f, %f\n", V_SS1, V_SR, V_RS, V_SI, V_ST, N);
  Rprintf("%f, %f, %f, %f, %f, %f\n", X_SS + V_SS1, X_SR + V_SR, X_RS + V_RS, X_SI + V_SI, X_ST + V_ST, N);
}

//end_rinit

//start_skel
double p1 = (X_IS + X_II + X_IT + X_IR) / N; // Prevalence of infection with virus 1
double p2 = (X_SI + X_II + X_TI + X_RI + V_SI) / N; // Prevalence of infection with virus 2

double omega = (2 * M_PI) / 52.25;

double beta1;
double beta2;

if (strcmp("sinusoidal_forcing", sens) == 0) {
  
  // Modelling seasonality with a sinusoidal wave:
  beta1 = Ri1 / (1.0 - (R10 + R120)) * (1 + b1 * cos(omega * (t - phi1))) * gamma1;
  beta2 = Ri2 / (1.0 - (R20 + R120)) * (1 + b2 * cos(omega * (t - phi2))) * gamma2;
  
} else {
  
  beta1 = Ri1 / (1.0 - (R10 + R120)) * exp(eta_ah1 * ah + eta_temp1 * temp) * gamma1; // Initial reproduction no of virus 1 (R10+R120: initial prop immune to v1)
  beta2 = Ri2 / (1.0 - (R20 + R120)) * exp(eta_ah2 * ah + eta_temp2 * temp) * gamma2; // Initial reproduction no of virus 2 (R20+R120: initial prop immune to v2)
  
}

double lambda1 = beta1 * p1; // Force of infection with virus 1
double lambda2 = beta2 * p2; // Force of infection with virus 2

double delta2 = d2 * delta1; // 1 / duration of refractory period (virus2 -> virus1)

// ODEs
DX_SS = -(lambda1 + lambda2) * X_SS; 
DX_IS = lambda1 * X_SS - (gamma1 + theta_lambda1 * lambda2) * X_IS + (1 - vacc_eff) * lambda1 * (V_SS1 + V_SS2); 
DX_TS = gamma1 * X_IS - (delta1 + theta_lambda1 * lambda2) * X_TS;
DX_RS = delta1 * X_TS - lambda2 * X_RS + delta_vacc * V_RS;
DX_SI = lambda2 * X_SS - (theta_lambda2 * lambda1 + gamma2) * X_SI;
DX_II = theta_lambda1 * lambda2 * X_IS + theta_lambda2 * lambda1 * X_SI - (gamma1 + gamma2) * X_II + (1 - vacc_eff) * theta_lambda2 * lambda1 * V_SI;
DX_TI = theta_lambda1 * lambda2 * X_TS + gamma1 * X_II - (delta1 + gamma2) * X_TI;
DX_RI = lambda2 * X_RS + delta1 * X_TI - gamma2 * X_RI + theta_lambda_vacc * lambda2 * V_RS; 
DX_ST = gamma2 * X_SI - (theta_lambda2 * lambda1 + delta2) * X_ST; 
DX_IT = gamma2 * X_II + theta_lambda2 * lambda1 * X_ST - (gamma1 + delta2) * X_IT + (1 - vacc_eff) * theta_lambda2 * lambda1 * V_ST;
DX_TT = gamma2 * X_TI + gamma1 * X_IT - (delta1 + delta2) * X_TT;
DX_RT = gamma2 * X_RI + delta1 * X_TT - delta2 * X_RT;
DX_SR = delta2 * X_ST - lambda1 * X_SR; 
DX_IR = delta2 * X_IT + lambda1 * X_SR - gamma1 * X_IR + (1 - vacc_eff) * lambda1 * V_SR; 
DX_TR = delta2 * X_TT + gamma1 * X_IR - delta1 * X_TR; 
DX_RR = delta2 * X_RT + delta1 * X_TR;

DV_SS1 = -(delta_vacc + theta_lambda_vacc * lambda2 + (1 - vacc_eff) * lambda1) * V_SS1;
DV_SS2 = delta_vacc * V_SS1 - ((1 - vacc_eff) * lambda1 + lambda2) * V_SS2;
DV_SI = theta_lambda_vacc * lambda2 * V_SS1 + lambda2 * V_SS2 - ((1 - vacc_eff) * theta_lambda2 * lambda1 + gamma2) * V_SI;
DV_ST = gamma2 * V_SI - ((1 - vacc_eff) * theta_lambda2 * lambda1 + delta2) * V_ST;
DV_SR = delta2 * V_ST - (1 - vacc_eff) * lambda1 * V_SR;
DV_RS = -(delta_vacc + theta_lambda_vacc * lambda2) * V_RS;

DH1 = gamma1 * p1; // Incidence rate of virus 1 infections (total)
DH2 = gamma2 * p2; // Incidence rate of virus 2 infections (total)
//DH1 = gamma1 * (X_IS + theta_rho2 * (X_II + X_IT) + X_IR) / N; // Incidence rate of reported virus 1 infections
//DH2 = gamma2 * (X_SI + theta_rho1 * (X_II + X_TI + V_SI) + X_RI) / N; // Incidence rate of reported virus 2 infections 
//end_skel
