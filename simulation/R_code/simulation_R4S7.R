options(mc.cores = parallel::detectCores())
#use_python("~/miniconda3/bin/python", required=TRUE)
#setwd("C:/Users/yg198/OneDrive - Duke University/2_cross_sectional_omics/0. Two step")
#mofa <- import("mofapy2")

library(tidyverse)
library(data.table)
library(splines)
library(survival)
library(rstan)
library(Matrix)
library(refund)
library(orthogonalsplinebasis)
library(mfaces)
library(statmod)
library(tdROC)
library(ipred)
library(Amelia)
rstan_options(auto_write = TRUE)

# seed = 1

source('source_code/functions.R')
source('source_code/data_gen.R')
source('source_code/data_gen_negbin.R')
source('source_code/sFPCA.R')
load('data_gen_x1.RData')
#tp$gamma_z = seq(0.3,0.7, length.out=15)
#tp$n_omics_factors = 15
#save(tp, file = "data_gen.RData")

#also changed the data_gen_ls function from source_code
#change factor numbers in survival model
N.train <- tp$N.train; N.test <- tp$N.test; N <- N.train + N.test
#N.train <- 300; N.test <- 100; N <- N.train + N.test;
ID <- 1:N; ID.train <- 1:N.train; ID.test <- (N.train+1):N


#tp$beta = c()
########First step: omics data
n_views = tp$n_omics_views
n_samples = N
n_features = tp$n_omics_features
n_factors = tp$n_omics_factors
n_total_factors = length(tp$gamma_z) + length(tp$gamma_x)

fname = paste0('MOFA_fit/R4S7/Z_', seed, '.RData')
load(fname)

##########Second step
dat <- data_gen_ls(seed, tp, N,Z_hat) ##generate longitudinal and survival data
long <- dat$long; long.train <- long[long$ID %in% ID.train, ]; long.test <- long[long$ID %in% ID.test, ]
surv <- dat$surv; surv.train <- surv[surv$ID %in% ID.train, ]; surv.test <- surv[surv$ID %in% ID.test, ]



J <- tp$J ## number of longitudinal outcomes
P <- 9 ## number of basis functions for B-spline
L0 <- tp$L0 ## npc for U_i(t)
L1 <- tp$L1 ## npc for W_ij(t)
tnew <- tp$tnew
tg <- length(tnew)


## fit Cox model for survival data ##
fit.cox <- coxph(Surv(surv_time, status) ~ x + 
                   z_hat.1 + z_hat.2 + z_hat.3 + z_hat.4 + z_hat.5 
                 #z_hat.6 + z_hat.7 + z_hat.8 + z_hat.9 + z_hat.10
                 , data = surv.train, method = 'breslow')
gamma_x_hat <- unname(fit.cox$coefficients)[1:n_total_factors]
#gamma_z_hat <- unname(fit.cox$coefficients)[-1]
baseline.cumulative.hazard <- basehaz(fit.cox)
lm.obj <- lm(hazard ~ time, data = baseline.cumulative.hazard)
logh0_hat <- log(unname(lm.obj$coefficients[2]))
x <- surv.train %>% dplyr::select(x, 
                                  z_hat.1:z_hat.5)
#z <- as.matrix(surv.train[,3:(3+n_factors-1)])
#############longitudinal data
## Impute the missing longitudinal outcomes, and average imputed datasets ##
M <- 5
long.impute <- amelia(x = long.train, m = M, idvars = 'ID', ts = 'time', splinetime = 6)
long.impute.dataset <- long.impute$imputations
long.i <- long.train
for (j in 1:J){
  tmp.y <- rep(0, nrow(long.train))
  for (m in 1:M){
    tmp.y <- tmp.y + long.impute.dataset[[m]][, j+2]
  }
  long.i[, j+2] <- tmp.y/M
}

Y.list.o <- list(Y1 = long.train$y1, Y2 = long.train$y2, Y3 = long.train$y3) ##original
Y.list.i <- list(Y1 = long.i$y1, Y2 = long.i$y2, Y3 = long.i$y3) #imputed

dat.mface <- list('y1' = data.frame('subj' = long.train$ID, 'argvals' = long.train$time, 'y' = Y.list.i[[1]]),
                  'y2' = data.frame('subj' = long.train$ID, 'argvals' = long.train$time, 'y' = Y.list.i[[2]]),
                  'y3' = data.frame('subj' = long.train$ID, 'argvals' = long.train$time, 'y' = Y.list.i[[3]]))
fit.mface <- mface.sparse(dat.mface, argvals.new = tnew, knots = 6, newdata = dat.mface)
#save(list = 'fit.mface', file = paste0('RData/fit_mface',"_seed",seed,'.RData'))

mu_est <- list()
for (i in 1:J){
  tmp.mu <- fit.mface$fit[[i]]$mu.new
  l <- bs.smooth(tmp.mu, tnew, tnew, nbasis = P)
  mu_est[[i]] <- list(value = tmp.mu, argvals = tnew, coefficient = l$coef.est)
}

## to calculate beta
C <- as.matrix(fit.mface$Chat.new)
C11 <- C[1:tg, 1:tg]
C12 <- C[1:tg, 1:tg+tg]
C13 <- C[1:tg, 1:tg+tg*2]
C22 <- C[1:tg+tg, 1:tg+tg]
C23 <- C[1:tg+tg, 1:tg+tg*2]
C33 <- C[1:tg+2*tg, 1:tg+2*tg]

beta_hat <- c(1, NA, NA)
beta_hat[2] <- sum(c(C13)*c(C23))/sum(c(C13)^2)
beta_hat[3] <- sum(c(C12)*c(C23))/sum(c(C12)^2)

coeff.C0 <- c(beta_hat[2], beta_hat[3],
              beta_hat[2]*beta_hat[3])
coeff.C1 <- beta_hat^2
C0_raw <- C1_raw <- matrix(NA, tg, tg)
for (t0 in 1:tg){
  for (t1 in 1:tg){
    C0_raw[t0, t1] <- sum(coeff.C0*c(C12[t0, t1], C13[t0, t1], 
                                     C23[t0, t1]))/sum(coeff.C0^2)
  }
}

C0_raw <- forceSymmetric(C0_raw)
C0_fit <- sFPCA_fit(C0_raw)
C0 <- C0_fit$S

for (t0 in 1:tg){
  for (t1 in 1:tg){
    C1_raw[t0, t1] <- sum(coeff.C1*c(C11[t0, t1] - beta_hat[1]^2*C0[t0, t1], 
                                     C22[t0, t1] - beta_hat[2]^2*C0[t0, t1],
                                     C33[t0, t1] - beta_hat[3]^2*C0[t0, t1]))/sum(coeff.C1^2)
  }
}

C1_fit <- sFPCA_fit(C1_raw)
C1 <- C1_fit$S
face.fit <- sFPCA_post(dat.mface, C0_fit, C1_fit, fit.mface, beta_hat, pve = 0.99)

d0_est <- face.fit$eigenvalues0[1:L0]
phi_est <- face.fit$eigenfunctions0[, 1:L0]
d1_est <- face.fit$eigenvalues1[1:L1]
psi_est <- matrix(face.fit$eigenfunctions1[, 1:L1], nrow = tg)


## Revert the sign of eigenfunctions for correct interpretation of eigenfunctions
phi_sign <- c(1, 1)
phi_index <- c(41, 1)   ## sqrt(2)*sin(pi*0.5)
phi <- matrix(NA, nrow(long.train), ncol = L0)
phi_surv <- matrix(NA, N.train, ncol = L0)
for (l in 1:L0){
  phi_est[, l] <- sign_eigen(phi_est[, l], phi_index[l], phi_sign[l])
  phi.smooth <- bs.smooth(phi_est[, l], tnew, argvals.new = long.train$time, nbasis = P)
  phi[, l] <- phi.smooth$est.value
  phi.smooth <- bs.smooth(phi_est[, l], tnew, argvals.new = surv.train$surv_time, nbasis = P)
  phi_surv[, l] <- phi.smooth$est.value
}

psi_sign <- c(1)
psi_index <- c(1)
psi <- matrix(NA, nrow(long.train), ncol = L1)
psi_surv <- matrix(NA, N.train, ncol = L1)
for (l in 1:L1){
  psi_est[, l] <- sign_eigen(psi_est[, l], psi_index[l], psi_sign[l])
  psi.smooth <- bs.smooth(psi_est[, l], tnew, argvals.new = long.train$time, nbasis = P)
  psi[, l] <- psi.smooth$est.value
  psi.smooth <- bs.smooth(psi_est[, l], tnew, argvals.new = surv.train$surv_time, nbasis = P)
  psi_surv[, l] <- psi.smooth$est.value
}

sigma_hat <- colMeans(sqrt(face.fit$var.error.new))
bs.mat <- create_bs(time.grid = tnew, pred.time = long.train$time, nbasis = P)
bs.mat.grid <- create_bs(time.grid = tnew, pred.time = tnew, nbasis = P)

## Determine missingness
Y_missing <- missing_determine(value = Y.list.o)
mean_Y <- list(mean1 = rep(mean(long.train$y1, na.rm = T), Y_missing$len.missing[1]), 
               mean2 = rep(mean(long.train$y2, na.rm = T), Y_missing$len.missing[2]),
               mean3 = rep(mean(long.train$y3, na.rm = T), Y_missing$len.missing[3]))
Y.list.new <- miss.index <- list()
for (j in 1:J){
  Y.list.new[[j]] <- Y_missing$new.value[[j]]
  miss.index[[j]] <- Y_missing$missing.index[[j]]
}

## weights; nodes; time vector of length N denoting interval (0, time[i]) ##
create_gq <- function(w, x, time){
  const <- rep(NA, length(time)) ## const_{i} = (b_{i} - a_{i})/2
  a <- b <- rep(NA, length(time))
  
  phi1_interval <- phi2_interval <- array(NA, dim = c(length(time), G))
  psi1_interval <- array(NA, dim = c(length(time), G))
  ##evaluated at G quadrature points for interval (0, surv_time[i])
  for (i in 1:length(time)){
    a[i] <- 0
    b[i] <- time[i]
    const[i] <- (b[i] - a[i])/2
    x_star <- (b[i] - a[i])/2*x + (b[i] + a[i])/2
    
    ## Evaluate interval of phi and psi on x_star
    phi.smooth <- bs.smooth(phi_est[, 1], tnew, argvals.new = x_star, nbasis = P)
    phi1_interval[i, ] <- phi.smooth$est.value
    phi.smooth <- bs.smooth(phi_est[, 2], tnew, argvals.new = x_star, nbasis = P)
    phi2_interval[i, ] <- phi.smooth$est.value
    
    psi.smooth <- bs.smooth(psi_est[, 1], tnew, argvals.new = x_star, nbasis = P)
    psi1_interval[i, ] <- psi.smooth$est.value
  }
  l <- list(const = const, 
            phi1_interval = phi1_interval, 
            phi2_interval = phi2_interval,
            psi1_interval = psi1_interval)
  return(l)
}

## Create Gaussian quadrature points on grid (a_{i}, b_{i}) ##
G <- 15 ## number of Gaussian quadrature points
weights <- gauss.quad(G, kind = 'legendre')$weights ## w_g^* in document
nodes <- gauss.quad(G, kind = 'legendre')$nodes  ## x_g^* in document

l <- create_gq(weights, nodes, surv.train$surv_time)
const <- l$const
phi1_interval <- l$phi1_interval
phi2_interval <- l$phi2_interval
psi1_interval <- l$psi1_interval


### stan sampling
md = stan_model('source_code/M1.stan')
stan_dat <- list(n = N.train, J = J, nobs = nrow(long.train), 
                 nmiss = Y_missing$len.missing, 
                 id_long = long.train$ID, 
                 L0 = L0, L1 = L1, P = P, 
                 P_surv = ncol(x), #P_z = ncol(Z_hat), 
                 Y1 = Y.list.new[[1]], Y2 = Y.list.new[[2]], Y3 = Y.list.new[[3]],
                 miss_index1 = as.array(miss.index[[1]]), 
                 miss_index2 = as.array(miss.index[[2]]), 
                 miss_index3 = as.array(miss.index[[3]]), 
                 time = long.train$time, 
                 x = x, #Z_hat = Z_hat,
                 surv_time = surv.train$surv_time, status = surv.train$status,
                 b = bs.mat, phi = phi, psi = psi)

rnd <- 1/3
inits1 <- list(A1 = g.inits(mu_est[[1]]$coefficient), A2 = g.inits(mu_est[[2]]$coefficient),
               A3 = g.inits(mu_est[[3]]$coefficient), 
               sqrt_d0 = as.array(g.inits.sd(sqrt(d0_est))), 
               sqrt_d1 = as.array(g.inits.sd(sqrt(d1_est))),
               beta = g.inits(beta_hat[2:J]), 
               omega = g.inits.sd(sigma_hat),
               logh0 = g.inits(logh0_hat), 
               gamma_x = as.array(g.inits(gamma_x_hat)), 
               #gamma_z = as.array(g.inits(gamma_z_hat)), 
               gamma0 = as.array(rep(0, L0)),
               gamma11 = as.array(rep(0, L1)), 
               gamma12 = as.array(rep(0, L1)), 
               gamma13 = as.array(rep(0, L1)), 
               Y1_imp = as.array(mean_Y[[1]]+rnd), 
               Y2_imp = as.array(mean_Y[[2]]+rnd), 
               Y3_imp = as.array(mean_Y[[3]]+rnd))
inits2 <- list(A1 = g.inits(mu_est[[1]]$coefficient), A2 = g.inits(mu_est[[2]]$coefficient),
               A3 = g.inits(mu_est[[3]]$coefficient), 
               sqrt_d0 = as.array(g.inits.sd(sqrt(d0_est))), 
               sqrt_d1 = as.array(g.inits.sd(sqrt(d1_est))),
               beta = g.inits(beta_hat[2:J]), 
               omega = g.inits.sd(sigma_hat),
               logh0 = g.inits(logh0_hat), 
               gamma_x = as.array(g.inits(gamma_x_hat)), 
               #gamma_z = as.array(g.inits(gamma_z_hat)),
               gamma0 = as.array(rep(0, L0)),
               gamma11 = as.array(rep(0, L1)), 
               gamma12 = as.array(rep(0, L1)), 
               gamma13 = as.array(rep(0, L1)), 
               Y1_imp = as.array(mean_Y[[1]]-rnd), 
               Y2_imp = as.array(mean_Y[[2]]-rnd), 
               Y3_imp = as.array(mean_Y[[3]]-rnd))
inits <- list(c1 = inits1, c2 = inits2)
#inits <- list(c1=inits1)
pars <- c('A1', 'A2', 'A3', 'd0', 'd1', 'beta', 
          'omega', 'logh0', 'gamma_x',#'gamma_z', 
          'gamma0', 
          'gamma11', 'gamma12', 'gamma13')
fitStan <- sampling(md, data = stan_dat, iter = 4000, warmup = 3000, 
                    chains = 2, thin=1, init = inits, pars = pars, seed = 2024*seed,
                    control = list(adapt_delta = 0.8, max_treedepth=10))
fname <- paste0("stan_fit/R4S7/stan_fit_", seed, '.RData',sep="")
save(list = 'fitStan', file = fname)
#load(fname)

Q <- 2000 
# A_sample <- extract(fitStan, pars = c('A1', 'A2', 'A3'), inc_warmup = F) ##coefficients for mu
# d0_sample <- extract(fitStan, pars = 'd0', inc_warmup = F)$d0
# d1_sample <- extract(fitStan, pars = 'd1', inc_warmup = F)$d1
# beta_sample <- extract(fitStan, pars = 'beta', inc_warmup = F)$beta
# omega_sample <- extract(fitStan, par = 'omega', inc_warmup = F)$omega
# logh0_sample <- extract(fitStan, pars = 'logh0', inc_warmup = F)$logh0
# gamma_x_sample <- extract(fitStan, pars = 'gamma_x', inc_warmup = F)$gamma_x
# gamma_z_sample <- extract(fitStan, pars = 'gamma_z', inc_warmup = F)$gamma_z
# gamma0_sample <- extract(fitStan, pars = 'gamma0', inc_warmup = F)$gamma0
# gamma1_sample <- extract(fitStan, pars = 'gamma1', inc_warmup = F)$gamma1

A_sample <- extract(fitStan, pars = c('A1', 'A2', 'A3')) ##coefficients for mu
d0_sample <- extract(fitStan, pars = 'd0')$d0
d1_sample <- extract(fitStan, pars = 'd1')$d1
beta_sample <- extract(fitStan, pars = 'beta')$beta
omega_sample <- extract(fitStan, par = 'omega')$omega
logh0_sample <- extract(fitStan, pars = 'logh0')$logh0
gamma_x_sample <- extract(fitStan, pars = 'gamma_x')$gamma_x
#gamma_z_sample <- extract(fitStan, pars = 'gamma_z')$gamma_z
gamma0_sample <- extract(fitStan, pars = 'gamma0')$gamma0
gamma11_sample <- extract(fitStan, pars = 'gamma11')$gamma11
gamma12_sample <- extract(fitStan, pars = 'gamma12')$gamma12
gamma13_sample <- extract(fitStan, pars = 'gamma13')$gamma13
mu1_hat <- mu2_hat <- mu3_hat <- matrix(0, Q, tg)
for (i in 1:Q){
  mu1_hat[i, ] <- (bs.mat.grid %*% A_sample[[1]][i, ])[, 1]
  mu2_hat[i, ] <- (bs.mat.grid %*% A_sample[[2]][i, ])[, 1]
  mu3_hat[i, ] <- (bs.mat.grid %*% A_sample[[3]][i, ])[, 1]
}

summary.parameters <- data.frame(mean = 0, sd = 0, lower_25_ci = 0, upper_975_ci = 0)
summary.parameters[1:L0, ] <- get_quantile_matrix(d0_sample)
summary.parameters[(L0+1):(L0+L1), ] <- get_quantile_matrix(d1_sample)
i1 <- L0 + L1
summary.parameters[(i1+1):(i1+J-1), ] <- get_quantile_matrix(beta_sample)
summary.parameters[(i1+J):(i1+J*2-1), ] <- get_quantile_matrix(omega_sample)

i2 <- i1 + J*2
summary.parameters[i2, ] <- get_quantile_vector(logh0_sample)
summary.parameters[(i2+1):(i2+n_total_factors),] <- get_quantile_matrix(gamma_x_sample)
#summary.parameters[(i2+2):(i2+2+n_factors-1), ] <- get_quantile_matrix(gamma_z_sample)
summary.parameters[(i2+n_total_factors+1):(i2+n_total_factors+L0), ] <- get_quantile_matrix(gamma0_sample)

i3 <- i2+n_total_factors+L0
summary.parameters[(i3+1):(i3+L1), ] <- get_quantile_matrix(gamma11_sample)
summary.parameters[(i3+L1+1):(i3+L1*2), ] <- get_quantile_matrix(gamma12_sample)
summary.parameters[(i3+L1*2+1):(i3+L1*3), ] <- get_quantile_matrix(gamma13_sample)

summary.mean <- data.frame(mean = 0)
summary.mean[1:tg, ] <- get_mean_function(mu1_hat)
summary.mean[(1+tg):(tg*2), ] <- get_mean_function(mu2_hat)
summary.mean[(1+tg*2):(tg*3), ] <- get_mean_function(mu3_hat)

## Dynamic Prediction ##
A <- lapply(1:J, FUN = function(x) colMeans(A_sample[[x]]))
d0 <- apply(d0_sample, 2, 'mean')
d1 <- apply(d1_sample, 2, 'mean')
beta <- apply(beta_sample, 2, 'mean')
omega <- sapply(1:J, FUN = function(x) mean(omega_sample[, x]))
logh0 <- mean(logh0_sample)
gamma_x <- apply(gamma_x_sample, 2, 'mean')
#gamma_z <- apply(gamma_z_sample, 2, 'mean')
gamma0 <- apply(gamma0_sample, 2, 'mean')
gamma11 <- apply(gamma11_sample, 2, 'mean')
gamma12 <- apply(gamma12_sample, 2, 'mean')
gamma13 <- apply(gamma13_sample, 2, 'mean')

## compute true survival probability at time t: exp(-cumulative hazard)
H_true <- function(tp, x, weights, nodes, 
                   xi, zeta1, zeta2, zeta3, t){
  f = exp(sum(tp$gamma0*xi) + sum(tp$gamma11*zeta1) + 
            sum(tp$gamma12*zeta2) + sum(tp$gamma13*zeta3))
  H <- exp(tp$logh0 + sum(x*c(tp$gamma_x, tp$gamma_z)))*t*f
  return(H)
}

## compute estimated survival probability at time t: exp(-cumulative hazard)
H_est <- function(x, weights, xi, zeta1, zeta2, zeta3, t){
  f = exp(sum(gamma0*xi) + sum(gamma11*zeta1) + sum(gamma12*zeta2)+ sum(gamma13*zeta3))
  H <- exp(logh0 + sum(x*gamma_x))*t*f
  return(H)
}

md_dp <- stan_model('source_code/M1_dp.stan')
pars_dp <- c('xi', 'zeta_1', 'zeta_2', 'zeta_3')

starting.time <- c(0.3, 0.4, 0.5, 0.55, 0.6) 
delta.time <- seq(0.1, 0.25, by = 0.01)
# par(mfrow = c(3, 3))
AUC <- BS <- matrix(NA, length(starting.time), length(delta.time))
AUC.true <- BS.true <- matrix(NA, length(starting.time), length(delta.time))

for (Tstart in starting.time){
  
  ## filter out subjects with observed survival time greater than starting time
  surv.test2 <- surv.test[surv.test$surv_time>Tstart, ]
  long.test2 <- long.test[long.test$ID %in% surv.test2$ID, ]
  surv.test2$newID <- 1:nrow(surv.test2)
  long.test2$newID <- rep(1:nrow(surv.test2), table(long.test2$ID))
  x.test <- surv.test2 %>% dplyr::select(x, 
                                         z_hat.1:z_hat.5)
  x.true = cbind(surv.test2$x, Z_true[surv.test2$ID, ])
  #z.test <- as.matrix(surv.test2[,3:(3+n_factors-1)])
  
  true_xi <- dat$xi[surv.test2$ID, ]
  true_zeta1 <- dat$zeta[[1]][surv.test2$ID, ]
  true_zeta2 <- dat$zeta[[2]][surv.test2$ID, ]
  true_zeta3 <- dat$zeta[[3]][surv.test2$ID, ]
  
  ## filter out longitudinal observations prior or equal to starting time
  long.test.prior <- long.test2[long.test2$time<=Tstart, ]
  long.test.posterior <- long.test2[long.test2$time>Tstart, ]
  tmp.ID.test <- unique(long.test2$newID)
  bs.mat.test <- create_bs(time.grid = tnew, pred.time = long.test.prior$time, nbasis = P)
  
  phi_test <- matrix(NA, nrow(long.test.prior), ncol = L0)
  for (l in 1:L0){
    phi.smooth <- bs.smooth(phi_est[, l], tnew, argvals.new = long.test.prior$time, nbasis = P)
    phi_test[, l] <- phi.smooth$est.value
  }
  
  psi_test <- matrix(NA, nrow(long.test.prior), ncol = L1)
  for (l in 1:L1){
    psi.smooth <- bs.smooth(psi_est[, l], tnew, argvals.new = long.test.prior$time, nbasis = P)
    psi_test[, l] <- psi.smooth$est.value
  }
  
  ## Determine missingness
  Y.list.test <- list(y1 = long.test.prior$y1, y2 = long.test.prior$y2, y3 = long.test.prior$y3)
  Y_missing.test <- missing_determine(value = Y.list.test)
  
  ## Create Gaussian quadrature points on grid (a_{i}, b_{i}) ##
  l_T0 <- create_gq(weights, nodes, time = Tstart)
  const_T0 <- l_T0$const
  phi1_interval_T0 <- l_T0$phi1_interval
  phi2_interval_T0 <- l_T0$phi2_interval
  psi1_interval_T0 <- l_T0$psi1_interval
  
  ## perform NUTS for posterior sampling of xi and zeta 
  stan_dat_dp <- list(n = length(tmp.ID.test), J = J, 
                      nobs = nrow(long.test.prior),
                      id_long = long.test.prior$newID,
                      L0 = L0, L1 = L1, P = P, 
                      P_surv = ncol(x.test), #P_z = ncol(z.test), 
                      Y1 = Y_missing.test$new.value[[1]], 
                      Y2 = Y_missing.test$new.value[[2]], 
                      Y3 = Y_missing.test$new.value[[3]],
                      time = long.test.prior$time, 
                      x = x.test, #z = z.test,
                      surv_time = rep(Tstart, length(tmp.ID.test)), 
                      b = bs.mat.test, phi = phi_test, psi = psi_test,
                      A1 = A[[1]], A2 = A[[2]], A3 = A[[3]],
                      sqrt_d0 = as.array(sqrt(d0)), sqrt_d1 = as.array(sqrt(d1)), 
                      beta = beta, omega = omega,
                      logh0 = logh0, 
                      gamma_x = as.array(gamma_x),
                      #gamma_z = as.array(gamma_z),
                      gamma0 = as.array(gamma0), 
                      gamma11 = as.array(gamma11), 
                      gamma12 = as.array(gamma12), 
                      gamma13 = as.array(gamma13))
  inits1_dp <- inits2_dp <- list(xi = matrix(0, L0, length(tmp.ID.test)), 
                                 zeta1 = matrix(0, L1, length(tmp.ID.test)), 
                                 zeta2 = matrix(0, L1, length(tmp.ID.test)), 
                                 zeta3 = matrix(0, L1, length(tmp.ID.test)))
  inits_dp <- list(c1 = inits1_dp, c2 = inits2_dp)
  fitStan_dp <- sampling(md_dp, data = stan_dat_dp, iter = 4000, warmup = 3000, 
                         chains = 2, thin=1, init = inits_dp, pars = pars_dp, seed = 123,
                         control = list(adapt_delta = 0.8, max_treedepth=10))
  
  xi_sample <- extract(fitStan_dp, pars = 'xi')$xi ## extract xi_{il}
  zeta1_sample <- extract(fitStan_dp, pars = 'zeta_1')$zeta_1 ## extract zeta1_{il}
  zeta2_sample <- extract(fitStan_dp, pars = 'zeta_2')$zeta_2 ## extract zeta2_{il}
  zeta3_sample <- extract(fitStan_dp, pars = 'zeta_3')$zeta_3
  
  ## calculate conditional survival probabilities at T+\delta_T
  for (Tdelta in delta.time){
    l_T1 <- create_gq(weights, nodes, time = Tstart+Tdelta)
    const_T1 <- l_T1$const
    phi1_interval_T1 <- l_T1$phi1_interval[1, ]
    phi2_interval_T1 <- l_T1$phi2_interval[1, ]
    psi1_interval_T1 <- l_T1$psi1_interval[1, ]
    
    tmp.surv.true <- rep(NA, length(tmp.ID.test))
    tmp.surv.predict <- rep(NA, length(tmp.ID.test))
    
    for (i in 1:length(tmp.ID.test)){
      xi_i_true <- true_xi[i, ]
      zeta1_i_true <- true_zeta1[i]
      zeta2_i_true <- true_zeta2[i]
      zeta3_i_true <- true_zeta3[i]
      
      Hi_true_t0 <- H_true(tp, x.test[i, ], #z.test[i, ], 
                           weights, nodes, 
                           xi_i_true, zeta1_i_true, zeta2_i_true, zeta3_i_true, 
                           t = Tstart)
      
      Hi_true_t1 <- H_true(tp, x.test[i, ], #z.test[i, ], 
                           weights, nodes, 
                           xi_i_true, zeta1_i_true, zeta2_i_true, zeta3_i_true, 
                           t = Tstart+Tdelta)
      
      tmp.surv.true[i] <- exp(-Hi_true_t1 + Hi_true_t0)
      
      ## estimated survival probability ##
      xi_i <- xi_sample[, , i] 
      zeta1_i <- zeta1_sample[, , i]
      zeta2_i <- zeta2_sample[, , i]
      zeta3_i <- zeta3_sample[, , i]
      Hi_est_t0 <- Hi_est_t1 <- rep(NA, Q)
      for (q in 1:Q){
        Hi_est_t0[q] <- H_est(x.test[i, ], #z.test[i, ], 
                              weights = weights, 
                              xi_i[q, ], zeta1_i[q], zeta2_i[q], zeta3_i[q], 
                              t = Tstart)
        Hi_est_t1[q] <- H_est(x.test[i, ], #z.test[i, ], 
                              weights = weights, 
                              xi_i[q, ], zeta1_i[q], zeta2_i[q], zeta3_i[q], 
                              t = Tstart + Tdelta)
      }
      
      cond_surv_prob <- exp(-Hi_est_t1 + Hi_est_t0)
      tmp.surv.predict[i] <- mean(cond_surv_prob)
    }
    
    ROC.est <- tdROC(X = 1 - tmp.surv.predict, Y = surv.test2$surv_time,
                     delta = surv.test2$status, tau = Tstart + Tdelta,
                     span = 0.1, alpha = 0.05,
                     n.grid = 1000, cut.off = 0.5)
    AUC[which(Tstart==starting.time), which(Tdelta==delta.time)] <- ROC.est$main_res$AUC.empirical
    
    ROC.true <- tdROC(X = 1 - tmp.surv.true, Y = surv.test2$surv_time, 
                      delta = surv.test2$status, tau = Tstart + Tdelta,
                      span = 0.1, alpha = 0.05,
                      n.grid = 1000, cut.off = 0.5)
    AUC.true[which(Tstart==starting.time), which(Tdelta==delta.time)] <- ROC.true$main_res$AUC.empirical
    
    surv.obj <- Surv(surv.test2$surv_time, surv.test2$status)
    BS[which(Tstart==starting.time), which(Tdelta==delta.time)] <- sbrier(surv.obj, tmp.surv.predict, btime = Tstart + Tdelta)
    BS.true[which(Tstart==starting.time), which(Tdelta==delta.time)] <- sbrier(surv.obj, tmp.surv.true, btime = Tstart + Tdelta)
  }
}

# fname <- paste0('result_M2_2/', "10050seed",seed, '.RData')
# save(list = c('summary.parameters', 'summary.mean', 'phi_est', 'psi_est',
#               'AUC', 'BS', 'AUC.true', 'BS.true'), file = fname)
# datatt<- load(fname)

result <- list(summary.parameters = summary.parameters, 
               summary.mean = summary.mean, 
               phi_est = phi_est, psi_est = psi_est,
               AUC = AUC, BS = BS, 
               AUC.true = AUC.true, BS.true = BS.true)

fname <- paste0("result/R4S7/", seed, '.RData',sep="")
save(result,file=fname)
#load(fname)
