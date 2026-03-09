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
T_theta_rho1 = log(theta_rho1);
T_theta_rho2 = log(theta_rho2);

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
T_b1 = logit(b1);
T_b2 = logit(b2);
T_phi1 = logitCons(phi1, 0.0, 52.25);
T_phi2 = logitCons(phi2, 0.0, 52.25);
T_beta_h3 = log(beta_h3);
T_beta_rhino = log(beta_rhino);
T_N = N;
T_beta_sd1 = log(beta_sd1);
T_beta_sd2 = log(beta_sd2);

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
theta_rho1 = exp(T_theta_rho1);
theta_rho2 = exp(T_theta_rho2);

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
b1 = expit(T_b1);
b2 = expit(T_b2);
phi1 = expitCons(T_phi1, 0.0, 52.25);
phi2 = expitCons(T_phi2, 0.0, 52.25);
beta_h3 = exp(T_beta_h3);
beta_rhino = exp(T_beta_rhino);
N = T_N;
beta_sd1 = exp(T_beta_sd1);
beta_sd2 = exp(T_beta_sd2);

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

double rho1_w;
double rho2_w;

if (strcmp("us", loc) == 0) {
  rho1_w = fmin2(1.0, rho1 * H1 / i_ILI); // Probability of detecting virus 1
  rho2_w = fmin2(1.0, rho2 * H2 / i_ILI); // Probability of detecting virus 2
} else {
  rho1_w = fmin2(1.0, rho1 * (1.0 + alpha * cos(omega * (t - phi))) * H1 / i_ILI); // Probability of detecting virus 1
  rho2_w = fmin2(1.0, rho2 * (1.0 + alpha * cos(omega * (t - phi))) * H2 / i_ILI); // Probability of detecting virus 2
}

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

double rho1_w;
double rho2_w;

if (strcmp("us", loc) == 0) {
  rho1_w = fmin2(1.0, rho1 * H1 / i_ILI); // Probability of detecting virus 1
  rho2_w = fmin2(1.0, rho2 * H2 / i_ILI); // Probability of detecting virus 2
} else {
  rho1_w = fmin2(1.0, rho1 * (1.0 + alpha * cos(omega * (t - phi))) * H1 / i_ILI); // Probability of detecting virus 1
  rho2_w = fmin2(1.0, rho2 * (1.0 + alpha * cos(omega * (t - phi))) * H2 / i_ILI); // Probability of detecting virus 2
}

n_P1 = rbinom(n_T1, rho1_w); // Generate of tests positive to virus 1
n_P2 = rbinom(n_T2, rho2_w); // Generate of tests positive to virus 2
//end_rmeas_testdiff

//start_rinit
X_SS = nearbyint((1.0 - I10 - I20 - R10 - R20 - R120) * N);
X_IS = nearbyint(I10 * N);
X_TS = 0;
X_RS = nearbyint(R10 * N);
X_SI = nearbyint(I20 * N);
X_II = 0;
X_TI = 0;
X_RI = 0;
X_ST = 0;
X_IT = 0;
X_TT = 0;
X_RT = 0;
X_SR = nearbyint(R20 * N);
X_IR = 0;
X_TR = 0;
X_RR = nearbyint(R120 * N);
H1_tot = 0;
H2_tot = 0;
H1 = 0; 
H2 = 0;

if ((X_SS + X_IS + X_RS + X_SI + X_SR + X_RR) != N) {
  X_SS = nearbyint(N - X_IS - X_RS - X_SI - X_SR - X_RR);
}

if ((X_SS + X_IS + X_RS + X_SI + X_SR + X_RR) != N) {
  Rprintf("SS=%f, IS=%f, RS=%f, SI=%f, SR=%f, RR=%f, sum=%f, N=%f\n", X_SS, X_IS, X_RS, X_SI, X_SR, X_RR, X_SS + X_IS + X_RS + X_SI + X_SR + X_RR, N);
}

if(debug) {
  Rprintf("%f, %f, %f, %f, %f, %f, %f\n", Ri1, Ri2, I10, I20, R10, R20, R120);
}

//end_rinit

//start_skel
double p1 = (X_IS + X_II + X_IT + X_IR) / N; // Prevalence of infection with virus 1
double p2 = (X_SI + X_II + X_TI + X_RI) / N; // Prevalence of infection with virus 2

double omega = (2 * M_PI) / 52.25;

double beta1;
double beta2;

if (strcmp("sinusoidal_forcing", sens) == 0) {
  
  // Modelling seasonality with a sinusoidal wave:
  beta1 = Ri1 / (1.0 - (R10 + R120)) * (1 + b1 * cos(omega * (t - phi1))) * gamma1;
  beta2 = Ri2 / (1.0 - (R20 + R120)) * (1 + b2 * cos(omega * (t - phi2))) * gamma2;
  
} else if (strcmp("h3_covar", sens) == 0) {
  
  // Include H3 as a covariate modulating RSV transmission
  // double beta_h3_t = fmax2(0.0, 1.0 - beta_h3 * h3_inc);
  double beta_h3_t = exp(-1 * beta_h3 * h3_inc);
  beta1 = Ri1 / (1.0 - (R10 + R120)) * exp(eta_ah1 * ah + eta_temp1 * temp) * gamma1;
  beta2 = Ri2 / (1.0 - (R20 + R120)) * beta_h3_t * exp(eta_ah2 * ah + eta_temp2 * temp) * gamma2;
  
} else if (strcmp("rhino_covar", sens) == 0) {
  
  // Include rhinovirus as a covariate modulating influenza transmission
  double beta_rhino_t = exp(-1 * beta_rhino * rhino_inc);
  beta1 = Ri1 / (1.0 - (R10 + R120)) * beta_rhino_t * exp(eta_ah1 * ah + eta_temp1 * temp) * gamma1;
  beta2 = Ri2 / (1.0 - (R20 + R120)) * exp(eta_ah2 * ah + eta_temp2 * temp) * gamma2;
  
} else {
  
  beta1 = Ri1 / (1.0 - (R10 + R120)) * exp(eta_ah1 * ah + eta_temp1 * temp) * gamma1; // Initial reproduction no of virus 1 (R10+R120: initial prop immune to v1)
  beta2 = Ri2 / (1.0 - (R20 + R120)) * exp(eta_ah2 * ah + eta_temp2 * temp) * gamma2; // Initial reproduction no of virus 2 (R20+R120: initial prop immune to v2)
  
  // Linear (vs. exponential) impact of climate:
  // beta1 = Ri1 / (1.0 - (R10 + R120)) * (1 + eta_ah1 * ah + eta_temp1 * temp) * gamma1; // Linear (vs. exponential) impact of climate
  // beta2 = Ri2 / (1.0 - (R20 + R120)) * (1 + eta_ah2 * ah + eta_temp2 * temp) * gamma2; // same
  
}

double lambda1 = beta1 * p1; // Force of infection with virus 1
double lambda2 = beta2 * p2; // Force of infection with virus 2

double delta2 = d2 * delta1; // 1 / duration of refractory period (virus2 -> virus1)

// ODEs
DX_SS = -(lambda1 + lambda2) * X_SS; 
DX_IS = lambda1 * X_SS - (gamma1 + theta_lambda1 * lambda2) * X_IS; 
DX_TS = gamma1 * X_IS - (delta1 + theta_lambda1 * lambda2) * X_TS;
DX_RS = delta1 * X_TS - lambda2 * X_RS; 
DX_SI = lambda2 * X_SS - (theta_lambda2 * lambda1 + gamma2) * X_SI;
DX_II = theta_lambda1 * lambda2 * X_IS + theta_lambda2 * lambda1 * X_SI - (gamma1 + gamma2) * X_II; 
DX_TI = theta_lambda1 * lambda2 * X_TS + gamma1 * X_II - (delta1 + gamma2) * X_TI;
DX_RI = lambda2 * X_RS + delta1 * X_TI - gamma2 * X_RI; 
DX_ST = gamma2 * X_SI - (theta_lambda2 * lambda1 + delta2) * X_ST; 
DX_IT = gamma2 * X_II + theta_lambda2 * lambda1 * X_ST - (gamma1 + delta2) * X_IT; 
DX_TT = gamma2 * X_TI + gamma1 * X_IT - (delta1 + delta2) * X_TT;
DX_RT = gamma2 * X_RI + delta1 * X_TT - delta2 * X_RT;
DX_SR = delta2 * X_ST - lambda1 * X_SR; 
DX_IR = delta2 * X_IT + lambda1 * X_SR - gamma1 * X_IR; 
DX_TR = delta2 * X_TT + gamma1 * X_IR - delta1 * X_TR;
DX_RR = delta2 * X_RT + delta1 * X_TR;

DH1_tot = gamma1 * p1; // Incidence rate of virus 1 infections (total)
DH2_tot = gamma2 * p2; // Incidence rate of virus 2 infections (total)
DH1 = gamma1 * (X_IS + theta_rho2 * (X_II + X_IT) + X_IR) / N; // Incidence rate of reported virus 1 infections
DH2 = gamma2 * (X_SI + theta_rho1 * (X_II + X_TI) + X_RI) / N; // Incidence rate of reported virus 2 infections 
//end_skel

//start_rsim

double omega = (2 * M_PI) / 52.25;

// Calculate prevalences:
double p1 = (X_IS + X_II + X_IT + X_IR) / N; // Prevalence of infection with virus 1
double p2 = (X_SI + X_II + X_TI + X_RI) / N; // Prevalence of infection with virus 2

// Calculate R0:
double R0_1;
double R0_2;

if (strcmp("sinusoidal_forcing", sens) == 0) {
  
  // Modelling seasonality with a sinusoidal wave:
  R0_1 = Ri1 / (1.0 - (R10 + R120)) * (1 + b1 * cos(omega * (t - phi1)));
  R0_2 = Ri2 / (1.0 - (R20 + R120)) * (1 + b2 * cos(omega * (t - phi2)));
  
} else if (strcmp("h3_covar", sens) == 0) {
  
  // Include H3 as a covariate modulating RSV transmission
  double beta_h3_t = exp(-1 * beta_h3 * h3_inc);
  R0_1 = Ri1 / (1.0 - (R10 + R120)) * exp(eta_ah1 * ah + eta_temp1 * temp);
  R0_2 = Ri2 / (1.0 - (R20 + R120)) * beta_h3_t * exp(eta_ah2 * ah + eta_temp2 * temp);
  
} else if (strcmp("rhino_covar", sens) == 0) {
  
  // Include rhinovirus as a covariate modulating influenza transmission
  double beta_rhino_t = exp(-1 * beta_rhino * rhino_inc);
  R0_1 = Ri1 / (1.0 - (R10 + R120)) * beta_rhino_t * exp(eta_ah1 * ah + eta_temp1 * temp);
  R0_2 = Ri2 / (1.0 - (R20 + R120)) * exp(eta_ah2 * ah + eta_temp2 * temp);
  
} else {
  
  R0_1 = Ri1 / (1.0 - (R10 + R120)) * exp(eta_ah1 * ah + eta_temp1 * temp); // Basic reproductive number (virus 1)
  R0_2 = Ri2 / (1.0 - (R20 + R120)) * exp(eta_ah2 * ah + eta_temp2 * temp); // Basic reproductive number (virus 2)
  
}

// Incorporate extrademographic stochasticity:
double beta1, beta2;
if (p1 > 0.0 && beta_sd1 > 0.0) {
  beta1 = rgammawn(sqrt(R0_1 / (p1 * N * beta_sd1 * dt)), R0_1 * gamma1);
} else {
  beta1 = R0_1 * gamma1;
}
if (p2 > 0.0 && beta_sd2 > 0.0) {
  beta2 = rgammawn(sqrt(R0_2 / (p2 * N * beta_sd2 * dt)), R0_2 * gamma2);
} else {
  beta2 = R0_2 * gamma2;
}

double lambda1 = beta1 * p1; // Force of infection with virus 1
double lambda2 = beta2 * p2; // Force of infection with virus 2

// Calculate duration of refractory period of virus 2:
double delta2 = d2 * delta1; // 1 / duration of refractory period (virus2 -> virus1)

// Calculate transitions:
double rates[18];
double fromSS[2], fromIS[2], fromTS[2], fromSI[2], fromII[2], fromTI[2], fromST[2], fromIT[2], fromTT[2];
double fromRS, fromRI, fromRT, fromSR, fromIR, fromTR;
//double reportII1, reportII2, reportIT, reportTI;

rates[0] = lambda1;
rates[1] = lambda2;
rates[2] = gamma1;
rates[3] = theta_lambda1 * lambda2;
rates[4] = delta1;
rates[5] = theta_lambda1 * lambda2;
rates[6] = theta_lambda2 * lambda1;
rates[7] = gamma2;
rates[8] = gamma1;
rates[9] = gamma2;
rates[10] = delta1;
rates[11] = gamma2;
rates[12] = theta_lambda2 * lambda1;
rates[13] = delta2;
rates[14] = gamma1;
rates[15] = delta2;
rates[16] = delta1;
rates[17] = delta2;

reulermultinom(2, X_SS, &rates[0], dt, &fromSS[0]);
reulermultinom(2, X_IS, &rates[2], dt, &fromIS[0]);
reulermultinom(2, X_TS, &rates[4], dt, &fromTS[0]);
reulermultinom(2, X_SI, &rates[6], dt, &fromSI[0]);
reulermultinom(2, X_II, &rates[8], dt, &fromII[0]);
reulermultinom(2, X_TI, &rates[10], dt, &fromTI[0]);
reulermultinom(2, X_ST, &rates[12], dt, &fromST[0]);
reulermultinom(2, X_IT, &rates[14], dt, &fromIT[0]);
reulermultinom(2, X_TT, &rates[16], dt, &fromTT[0]);

// Rprintf("first=%.1f, second=%.f, first=%.1f, second=%.f\n", fromSS[0], fromSS[1], fromIS[0], fromIS[1]);

fromRS = rbinom(X_RS, pTrans(lambda2, dt));
fromRI = rbinom(X_RI, pTrans(gamma2, dt));
fromRT = rbinom(X_RT, pTrans(delta2, dt));
fromSR = rbinom(X_SR, pTrans(lambda1, dt));
fromIR = rbinom(X_IR, pTrans(gamma1, dt));
fromTR = rbinom(X_TR, pTrans(delta1, dt));

//reportII1 = rbinom(fromII[0], theta_rho2);
//reportIT = rbinom(fromIT[0], theta_rho2);
//reportII2 = rbinom(fromII[1], theta_rho1);
//reportTI = rbinom(fromTI[1], theta_rho1);

// Balance equations:
X_SS += -fromSS[0] - fromSS[1];
X_IS += fromSS[0] - fromIS[0] - fromIS[1];
X_TS += fromIS[0] - fromTS[0] - fromTS[1];
X_RS += fromTS[0] - fromRS;

X_SI += fromSS[1] - fromSI[0] - fromSI[1];
X_II += fromIS[1] + fromSI[0] - fromII[0] - fromII[1];
X_TI += fromTS[1] + fromII[0] - fromTI[0] - fromTI[1];
X_RI += fromRS + fromTI[0] - fromRI;

X_ST += fromSI[1] - fromST[0] - fromST[1];
X_IT += fromII[1] + fromST[0] - fromIT[0] - fromIT[1];
X_TT += fromTI[1] + fromIT[0] - fromTT[0] - fromTT[1];
X_RT += fromRI + fromTT[0] - fromRT;

X_SR += fromST[1] - fromSR;
X_IR += fromIT[1] + fromSR - fromIR;
X_TR += fromTT[1] + fromIR - fromTR;
X_RR += fromRT + fromTR;

H1_tot += (fromIS[0] + fromII[0] + fromIT[0] + fromIR) / N;
H2_tot += (fromSI[1] + fromII[1] + fromTI[1] + fromRI) / N;
H1 += (fromIS[0] + theta_rho2 * (fromII[0] + fromIT[0]) + fromIR) / N;
H2 += (fromSI[1] + theta_rho1 * (fromII[1] + fromTI[1]) + fromRI) / N;
//H1 += (fromIS[0] + reportII1 + reportIT + fromIR) / N;
//H2 += (fromSI[1] + reportII2 + reportTI + fromRI) / N;

//end_rsim
