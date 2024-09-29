data{
  int n; // number of subjects in testing dataset
  int J;
  int nobs; // number of observations prior or equal to time T in testing dataset
  int id_long[nobs]; // ID
  int L0; // number of eigenvalues for U_i(t)
  int L1; // number of eigenvalues for W_ij(t)
  int P; // number of cubic B-spline functions
  int P_surv;
  
  real Y1[nobs]; // Y1
  real Y2[nobs]; // Y2
  real Y3[nobs];
  real time[nobs]; // observed time
  matrix[n, P_surv] x;
  real surv_time[n]; // observed survival time
  matrix[nobs, P] b; // orthogonal cubic B-spline matrix
  matrix[nobs, L0] phi; // estimated eigenfunctions for U_i(t)
  matrix[nobs, L1] psi; // estimated eigenfunctions for W_ij(t)
  
  vector[P] A1; // estimated coefficients for mu1
  vector[P] A2; // estimated coefficients for mu2
  vector[P] A3;
  vector<lower=0>[L0] sqrt_d0; // estimated standard deviation for FPC scores xi_il
  vector<lower=0>[L1] sqrt_d1; // estimated standard deviation for FPC scores zeta_ijl

  vector[J-1] beta;  // estimated beta 
  vector[J] omega;
  
  real logh0; 
  vector[P_surv] gamma_x;
  row_vector[L0] gamma0; 
  row_vector[L1] gamma11;
  row_vector[L1] gamma12;
  row_vector[L1] gamma13;
}
parameters{
  matrix[L0, n] xi; // FPC scores for U_i(t)
  matrix[L1, n] zeta_1; // FPC scores for W_i1(t)
  matrix[L1, n] zeta_2; // FPC scores for W_i2(t)
  matrix[L1, n] zeta_3; // FPC scores for W_i3(t)
}
transformed parameters{
  real mu1[nobs]; // mean for Y1
  real mu2[nobs]; // mean for Y2
  real mu3[nobs];
  
  real ll1[nobs];
  real ll2[nobs];
  real ll3[nobs];
  
  vector[L0] tmp_xi; // xi[, i]: FPC score for U_i(t)
  vector[L1] tmp_zeta1; // zeta1[, i]: FPC scores for W_i1(t)
  vector[L1] tmp_zeta2; // zeta2[, i]: FPC scores for W_i2(t)
  vector[L1] tmp_zeta3; // zeta3[, i]: FPC scores for W_i3(t)
  
  real rs_w[n]; // risk score for subject i
  real H[n]; // cumulative hazard 
  real LL[n]; // log survival likelihood

  for (i in 1:nobs){
    tmp_xi = col(xi, id_long[i]); // FPC scores xi_{i} for subject i
    
    tmp_zeta1 = col(zeta_1, id_long[i]); // FPC scores zeta_{i1} for subject i
    tmp_zeta2 = col(zeta_2, id_long[i]); // FPC scores zeta_{i2} for subject i
    tmp_zeta3 = col(zeta_3, id_long[i]);
    
    mu1[i] = b[i]*A1 + phi[i]*tmp_xi + psi[i]*tmp_zeta1;
    mu2[i] = b[i]*A2 + beta[1]*(phi[i]*tmp_xi + psi[i]*tmp_zeta2);
    mu3[i] = b[i]*A3 + beta[2]*(phi[i]*tmp_xi + psi[i]*tmp_zeta3);
    
    if (Y1[i]==-100) ll1[i] = 0; else ll1[i] = normal_lpdf(Y1[i] | mu1[i], omega[1]);
    if (Y2[i]==-100) ll2[i] = 0; else ll2[i] = normal_lpdf(Y2[i] | mu2[i], omega[2]);
    if (Y3[i]==-100) ll3[i] = 0; else ll3[i] = normal_lpdf(Y3[i] | mu3[i], omega[3]);
  }
  
  for (i in 1:n){
    tmp_xi = col(xi, i);
    tmp_zeta1 = col(zeta_1, i);
    tmp_zeta2 = col(zeta_2, i);
    tmp_zeta3 = col(zeta_3, i);
    
    rs_w[i] = x[i]*gamma_x + gamma0*tmp_xi + gamma11*tmp_zeta1 + 
              gamma12*tmp_zeta2 + gamma13*tmp_zeta3;
    H[i] = exp(logh0 +  rs_w[i])*surv_time[i];
    
    LL[i] = (-H[i]);
  }
}
model{
  for (i in 1:L0){
    xi[i] ~ normal(0, sqrt_d0[i]);
  }
  
  for (i in 1:L1){
    zeta_1[i] ~ normal(0, sqrt_d1[i]);
    zeta_2[i] ~ normal(0, sqrt_d1[i]);
    zeta_3[i] ~ normal(0, sqrt_d1[i]);
  }
  
  target+=sum(ll1);
  target+=sum(ll2);
  target+=sum(ll3);
  target+=LL;
}
