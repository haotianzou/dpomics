data{
  int n; // number of subjects
  int J;
  int nobs; // number of observations
  int nmiss[J]; 
  int id_long[nobs]; // ID
  int L0; // number of eigenvalues for U_i(t)
  int L1; // number of eigenvalues for W_ij(t)
  int P; // number of cubic B-spline functions
  int P_surv; // number of covariates in survival model
  
  real Y1[nobs]; // Y1
  real Y2[nobs]; // Y2
  real Y3[nobs];
  int miss_index1[nmiss[1]];  // missing index for Y1
  int miss_index2[nmiss[2]];
  int miss_index3[nmiss[3]]; 
  real time[nobs]; // observed time
  
  matrix[n, P_surv] x;
  
  real surv_time[n]; // survival time
  real status[n]; // censoring status
  matrix[nobs, P] b; // cubic B-spline matrix
  matrix[nobs, L0] phi; // estimated eigenfunctions for U_i(t)
  matrix[nobs, L1] psi; // estimated eigenfunctions for W_ij(t)
}
parameters{
  vector[P] A1; // coefficients for mu1
  vector[P] A2; // coefficients for mu2
  vector[P] A3;
  
  vector<lower=0>[L0] sqrt_d0; // standard deviation for FPC scores xi_il
  vector<lower=0>[L1] sqrt_d1; // standard deviation for zeta_ijl
  
  matrix[L0, n] xi; // FPC scores for U_i(t)
  matrix[L1, n] zeta_1; // FPC scores for W_i1(t)
  matrix[L1, n] zeta_2; // FPC scores for W_i2(t)
  matrix[L1, n] zeta_3;
  
  vector[J-1] beta; // beta
  real<lower=0> omega[J]; // scale parameter for error
  
  real logh0; 
  vector[P_surv] gamma_x; // coefficients for x[1]-x[4]
  row_vector[L0] gamma0; // coefficients for xi
  row_vector[L1] gamma11; // coefficients for zeta_1
  row_vector[L1] gamma12;
  row_vector[L1] gamma13;
  
  real Y1_imp[nmiss[1]]; // imputed Y1 
  real Y2_imp[nmiss[2]];
  real Y3_imp[nmiss[3]];
}
transformed parameters{
  real Y1_full[nobs] = Y1; // Full Y1 (observed + imputed)
  real Y2_full[nobs] = Y2;
  real Y3_full[nobs] = Y3;
  
  real mu1[nobs]; // m_{i1}(t)
  real mu2[nobs];
  real mu3[nobs];
  
  vector[L0] d0; // variance for xi
  vector[L1] d1; // variance for zeta
  vector[L0] tmp_xi; // xi[, i]
  
  vector[L1] tmp_zeta1;
  vector[L1] tmp_zeta2;
  vector[L1] tmp_zeta3;
  
  real rs_w[n]; // risk score for subject i
  real h[n]; // hazard function
  real H[n]; // cumulative hazard
  real LL[n]; // log survival likelihood
  
  Y1_full[miss_index1] = Y1_imp;
  Y2_full[miss_index2] = Y2_imp;
  Y3_full[miss_index3] = Y3_imp;
  
  for (i in 1:L0){
    d0[i] = pow(sqrt_d0[i], 2);
  }
  
  for (i in 1:L1){
    d1[i] = pow(sqrt_d1[i], 2);
  }
  
  for (i in 1:nobs){
    tmp_xi = col(xi, id_long[i]); // FPC scores xi_{i} for subject i
    
    tmp_zeta1 = col(zeta_1, id_long[i]); // FPC scores zeta_{i1} for subject i
    tmp_zeta2 = col(zeta_2, id_long[i]); // FPC scores zeta_{i2} for subject i
    tmp_zeta3 = col(zeta_3, id_long[i]);
    
    mu1[i] = b[i]*A1 + phi[i]*tmp_xi + psi[i]*tmp_zeta1;
    mu2[i] = b[i]*A2 + beta[1]*(phi[i]*tmp_xi + psi[i]*tmp_zeta2);
    mu3[i] = b[i]*A3 + beta[2]*(phi[i]*tmp_xi + psi[i]*tmp_zeta3);
  }
  
  for (i in 1:n){
    tmp_xi = col(xi, i);
    
    tmp_zeta1 = col(zeta_1, i);
    tmp_zeta2 = col(zeta_2, i);
    tmp_zeta3 = col(zeta_3, i);
    
    rs_w[i] = x[i]*gamma_x + 
              gamma0*tmp_xi + gamma11*tmp_zeta1 + gamma12*tmp_zeta2 + gamma13*tmp_zeta3;
    
    h[i] = exp(logh0 +  rs_w[i]);
    H[i] = exp(logh0 +  rs_w[i])*surv_time[i];
    
    LL[i] = status[i]*log(h[i]) + (-H[i]);
  }
}
model{
  A1 ~ normal(0, 10); // prior 
  A2 ~ normal(0, 10);
  A3 ~ normal(0, 10);
  
  for (i in 1:L0){
    sqrt_d0[i] ~ inv_gamma(0.1, 0.1);
    xi[i] ~ normal(0, sqrt_d0[i]);
  }
  
  for (i in 1:L1){
    sqrt_d1[i] ~ inv_gamma(0.1, 0.1);
    zeta_1[i] ~ normal(0, sqrt_d1[i]);
    zeta_2[i] ~ normal(0, sqrt_d1[i]);
    zeta_3[i] ~ normal(0, sqrt_d1[i]);
  }
  
  beta ~ normal(0, 10);
  omega ~ inv_gamma(0.1, 0.1);
  
  logh0 ~ normal(0, 10);
  gamma_x ~ normal(0, 10);
  gamma0 ~ normal(0, 10);
  gamma11 ~ normal(0, 10);
  gamma12 ~ normal(0, 10);
  gamma13 ~ normal(0, 10);
  
  target+=normal_lpdf(Y1_full | mu1, omega[1]); // likelihood
  target+=normal_lpdf(Y2_full | mu2, omega[2]);
  target+=normal_lpdf(Y3_full | mu3, omega[3]);
  target+=LL; // survival log likelihood
}
