library(tidyverse)
library(gridExtra)
library(latex2exp)
library(xtable)
library(data.table)
library(rstan)


J <- 3
P <- 9 

datadir <- 'result_M2_2/'
nrep = 100
resname0 <- "1241_4000_3000_200_cph11_nc"
resname <- paste0("1241_4000_3000_200_cph11_nc",1:nrep)
##actually 400 now
##nc, cw: with _; cp_nc, cp_cw: no _

load('data_gen.RData')
tp$omega <- c(1.2,1,1); tp$d0 <- c(5,2); tp$d1 <- 5
tp$gamma0 <- 0.4
tp$beta <- c(1,-1.5,0.4)
tp$N.train = 300; tp$N.test = 100
tp$n_omics_factors = 5
tp$gamma_z <- c(rep(0.5,2),rep(1,3))

mu <- list(
  f1 = function(t) 10 + 5*t + 5*t^2,
  f2 = function(t) 10 + 10*sqrt(t) + 5*sin(2*pi*t),
  f3 = function(t) sin(2*pi*t) + cos(2*pi*t) + log(t+1)
)

result.parameters <- result.mean <- result.phi <- result.psi <- list()
result.AUC <- result.AUC.true <- result.BS <-  result.BS.true <- list()
fitStan <- fitStan_dp <- list()
index.res <- 1
#fname <- paste0(resname, '.RData')
fname <- paste0(datadir, resname, '.RData')
for (index.res in 1:nrep){
 # if (file.exists(fname)){
    load(fname[index.res])
    fitStan[[index.res]] <- result[[1]]$fitStan
    fitStan_dp[[index.res]] <- result[[1]]$fitStan_dp
    result.parameters[[index.res]] <- result[[1]]$summary.parameters
    result.mean[[index.res]] <- result[[1]]$summary.mean
    result.phi[[index.res]] <- result[[1]]$phi_est
    result.psi[[index.res]] <- result[[1]]$psi_est
    result.AUC[[index.res]] <- result[[1]]$AUC
    result.BS[[index.res]] <- result[[1]]$BS
    result.AUC.true[[index.res]] <- result[[1]]$AUC.true
    result.BS.true[[index.res]] <- result[[1]]$BS.true
    #if (index.res==100) break;
    index.res <- index.res + 1
  #}
}

##check diagnostic: 
#https://cran.r-project.org/web/packages/rstan/vignettes/stanfit-objects.html
#https://jrnold.github.io/bayesian_notes/mcmc-diagnostics.html
# a<-fitStan[[1]]
# fit_summary <- rstan::summary(a)
# print(names(fit_summary))
# print(fit_summary$summary)
# #range(fit_summary$summary[,10])
# traceplot(a,pars=c(gamma_z))
# ggsave(file="sim_tracep_1241_4000_2001.eps",
#        width=6,height = 4)


#Sampler diagnostics
# sampler_params <- get_sampler_params(a, inc_warmup = FALSE)
# sampler_params_chain1 <- sampler_params[[1]]
# colnames(sampler_params_chain1)
# 
# mean_accept_stat_by_chain <- sapply(sampler_params, function(x) mean(x[, "accept_stat__"]))
# print(mean_accept_stat_by_chain)
# 
# max_treedepth_by_chain <- sapply(sampler_params, function(x) max(x[, "treedepth__"]))
# print(max_treedepth_by_chain)


get_summary <- function(result, row.index, true.value){
  tmp.mean <- rep(NA, length(result))
  tmp.se <- rep(NA, length(result))
  tmp.cp <- 0
  for (i in 1:length(result)){
    tmp.mean[i] <- result[[i]][row.index, 1]
    tmp.se[i] <- result[[i]][row.index, 2]
    if (result[[i]][row.index, 3]<true.value & result[[i]][row.index, 4]>true.value){
      tmp.cp <- tmp.cp + 1
    }
  }
  v <- c(mean = mean(tmp.mean), 
         bias = mean(tmp.mean) - true.value, 
         sd = sd(tmp.mean), 
         se = sqrt(sum(tmp.se^2)/length(result)),
         cp = tmp.cp/length(result))
  return(v)
}



true_value <- c(tp$d0, tp$d1, tp$beta[2:3], 
                tp$omega, tp$logh0, 
                tp$gamma_x, 
                tp$gamma_z, 
                tp$gamma0, tp$gamma11, tp$gamma12, tp$gamma13)
summ.dat <- data.frame(true.values = 0, mean = 0, bias = 0, sd = 0, se = 0, cp = 0)
for (i in 1:length(true_value)){
  summ.dat[i, 2:6] <- get_summary(result.parameters, i, true_value[i])
}
summ.dat$true.values <- true_value

rownames(summ.dat) <- c('d01', 'd02', 'd11', 'beta_2', 'beta_3',
                        'omega_1', 'omega_2', 'omega_3', 
                        'logh0', 
                        'gamma_x',paste0("gamma_z",1:tp$n_omics_factors),
                        'gamma0',  
                        'gamma11', 'gamma12', 'gamma13')

## Summary of functions ##
tnew <- (0:100)/100
tg <- length(tnew)

L0 <- tp$L0
L1 <- tp$L1

mu_tnew <- list(mu1 = mu[[1]](tnew), mu2 = mu[[2]](tnew), mu3 = mu[[3]](tnew))


## get estimated function and 95% CI, and plot #
func_est <- function(result, time.grid, function.index, index, B, lab.x, lab.y, range.y){
  mse <- rep(0, length(result))
  B_est <- matrix(NA, length(result), length(time.grid))
  B_summary <- data.frame(mean = rep(0, length(time.grid)), 
                          lower = rep(0, length(time.grid)), 
                          upper = rep(0, length(time.grid)))
  for (i in 1:length(result)){
    B_est[i, ] <- result[[i]][index, function.index]
    mse[i] <- sum((B - B_est[i, ])^2)
  }
  for (s in 1:length(time.grid)){
    B_summary$mean[s] <- mean(B_est[, s])
    B_summary$lower[s] <- quantile(B_est[, s], 0.025)
    B_summary$upper[s] <- quantile(B_est[, s], 0.975)
  }
  B_summary$true <- B
  p <- ggplot(data = B_summary, aes(x = time.grid, y = true)) + 
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = 'grey70') + 
    geom_line(aes(y = true), linetype = 'solid', size = 0.7) + 
    geom_line(aes(y = mean), linetype = 'dashed', col = 'blue', size = 0.7) + 
    theme_minimal() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "black", fill=NA, size=1),
          plot.margin = margin(0.5, 0.25, 0, 0.25, "cm")) +
    xlab(TeX(lab.x)) +
    ylab(TeX(lab.y)) + 
    ylim(range.y[1], range.y[2])
  l <- list(p = p, mse = mean(mse)/length(result))
  return(l)
}

l1 <- func_est(result.mean, tnew, 1, 1:tg, mu_tnew[[1]], lab.x = 'Time', lab.y = '$mu_1(t)$', c(9, 22))
l2 <- func_est(result.mean, tnew, 1, 1:tg+tg, mu_tnew[[2]], lab.x = 'Time', lab.y = '$mu_2(t)$', c(9, 25))
l3 <- func_est(result.mean, tnew, 1, 1:tg+tg*2, mu_tnew[[3]], lab.x = 'Time', lab.y = '$mu_3(t)$', c(-2, 3))

amse.mu1 <- l1[[2]]
amse.mu2 <- l2[[2]]
amse.mu3 <- l3[[2]]

setEPS()
postscript(paste0("plot/sim_mu_2_",resname0,".eps"), height = 4)
grid.arrange(l1[[1]], l2[[1]], l3[[1]], nrow = 1)
dev.off()

phi <- list(
  f1 = function(t) sqrt(2)*sin(pi*t),
  f2 = function(t) sqrt(2)*cos(pi*t)
)

phi_obs <- list(f1 = phi[[1]](tnew), f2 = phi[[2]](tnew))

psi <- list(
  f1 = function(t) sqrt(2)*cos(2*pi*t)
)

psi_obs <- list(f1 = psi[[1]](tnew))

phi1 <- func_est(result.phi, tnew, 1, 1:tg, phi_obs[[1]], lab.x = 'Time', lab.y = '$\\phi_1(t)$', c(-1, 2))
phi2 <- func_est(result.phi, tnew, 2, 1:tg, phi_obs[[2]], lab.x = 'Time', lab.y = '$\\phi_2(t)$', c(-2, 2))
amse.phi1 <- phi1$mse
amse.phi2 <- phi2$mse

psi1 <- func_est(result.psi, tnew, 1, 1:tg, psi_obs[[1]], lab.x = 'Time', lab.y = '$\\psi_1(t)$', c(-2, 2))
amse.psi1 <- psi1$mse


setEPS()
postscript(paste0("plot/sim_fpred_2_",resname0,".eps"), height = 4)
grid.arrange(phi1[[1]], phi2[[1]], psi1[[1]], nrow=1)
dev.off()


summ.dat[1:6, 7] <- c('AMSE.mu1', 'AMSE.mu2', 'AMSE.mu3', 'AMSE.phi1', 'AMSE.phi2', 'AMSE.psi1')
summ.dat[1:6, 8] <- c(amse.mu1, amse.mu2, amse.mu3, amse.phi1, amse.phi2, amse.psi1)
write.csv(summ.dat, file = paste0("RData/mfmm_jm_simulation_parameters_2_",resname0,".csv"))
#xtable(summ.dat, type = "latex", digits = c(0, rep(3, 6), 0, 3))

## iAUC is calculated via integration over \delta_t: (0.1, 0.11, ..., 0.25)
## Integration is calculated by Simpson's Rule:
## For 13 points: (0.1, ..., 0.22), we adopt Composite Simpson's Rule
## For 4 points: (0.22, 0.23 ,0.24, 0.25), we adopt 3/8 Simpson's Rule
composite.simpson <- function(vec, grid.points){
  n <- length(grid.points) - 1
  h <- grid.points[2] - grid.points[1]
  int.approx <- 0
  
  ## vector for R starts with index 1. So for even index, multiply by 4.
  for (j in 1:(n/2)){
    int.approx <- int.approx + 4*vec[2*j]
  }
  for (j in 1:(n/2-1)){
    int.approx <- int.approx + 2*vec[2*j+1]
  }
  int.approx <- int.approx + vec[1] + vec[n+1]
  int.approx <- h/3*int.approx
  return(int.approx)
}


# For the remaining 4 points #
composite.simpson.2 <- function(vec.remain, grid.points){
  h <- grid.points[2] - grid.points[1]
  int.approx <- 3/8*h*(vec.remain[1] + 3*vec.remain[2] + 3*vec.remain[3] + vec.remain[4])
  return(int.approx) 
}

starting.time <- c(0.3, 0.4, 0.5, 0.55, 0.6) 
delta.time <- seq(0.1, 0.25, by = 0.01)
l.time <- length(starting.time)
l.delta.time <- length(delta.time)
tg.1 <- delta.time[1:(l.delta.time-3)]
tg.2 <- delta.time[(l.delta.time-3):(l.delta.time)]

AUC.matrix <- BS.matrix <- matrix(0, l.time, l.delta.time)
AUC.true.matrix <- BS.true.matrix <- matrix(0, l.time, l.delta.time)
integrated_AUC_BS <- matrix(0, l.time, 4)
index1 <- 0
for (i in 1:(index.res-1)){
  tmp.AUC <- result.AUC[[i]]
  tmp.AUC.true <- result.AUC.true[[i]]
  tmp.BS <- result.BS[[i]]
  tmp.BS.true <- result.BS.true[[i]]
  
  if (sum(is.na(tmp.AUC))==0 & sum(is.na(tmp.AUC.true))==0){
    index1 <- index1 + 1
    AUC.matrix <- AUC.matrix + tmp.AUC
    BS.matrix <- BS.matrix + tmp.BS
    AUC.true.matrix <- AUC.true.matrix + tmp.AUC.true
    BS.true.matrix <- BS.true.matrix + tmp.BS.true
    
    for (j in 1:l.time){
      integrated_AUC_BS[j, 1] <- integrated_AUC_BS[j, 1] + 
        composite.simpson(tmp.AUC.true[j, 1:(l.delta.time-3)], tg.1) + 
        composite.simpson.2(tmp.AUC.true[j, (l.delta.time-3):(l.delta.time)], tg.2)
      integrated_AUC_BS[j, 2] <- integrated_AUC_BS[j, 2] + 
        composite.simpson(tmp.AUC[j, 1:(l.delta.time-3)], tg.1) + 
        composite.simpson.2(tmp.AUC[j, (l.delta.time-3):(l.delta.time)], tg.2)
      integrated_AUC_BS[j, 3] <- integrated_AUC_BS[j, 3] + 
        composite.simpson(tmp.BS.true[j, 1:(l.delta.time-3)], tg.1) + 
        composite.simpson.2(tmp.BS.true[j, (l.delta.time-3):(l.delta.time)], tg.2)
      integrated_AUC_BS[j, 4] <- integrated_AUC_BS[j, 4] + 
        composite.simpson(tmp.BS[j, 1:(l.delta.time-3)], tg.1) + 
        composite.simpson.2(tmp.BS[j, (l.delta.time-3):(l.delta.time)], tg.2)
    }
    
  }
  
}

AUC.matrix <- AUC.matrix/index1
AUC.true.matrix <- AUC.true.matrix/index1
BS.matrix <- BS.matrix/index1
BS.true.matrix <- BS.true.matrix/index1
integrated_AUC_BS <- integrated_AUC_BS/index1/(max(delta.time) - min(delta.time))

rownames(AUC.matrix) <- rownames(AUC.true.matrix) <- rownames(BS.matrix) <- rownames(BS.true.matrix) <- paste0('Landmark Time = ', starting.time)
colnames(AUC.matrix) <- colnames(AUC.true.matrix) <- colnames(BS.matrix) <- colnames(BS.true.matrix) <- paste0('Window = ', delta.time)
rownames(integrated_AUC_BS) <- paste0('T = ', starting.time)
colnames(integrated_AUC_BS) <- c('True iAUC', 'Estimated iAUC', 'True iBS', 'Estimated iBS')

xtable(integrated_AUC_BS, digits = c(0, rep(3, 4)))

write.table(AUC.true.matrix, file = paste0("RData/AUC_BS_2_",resname0,".csv"), append = FALSE, row.names = TRUE, col.names = TRUE, sep = ',')
write.table(AUC.matrix, file = paste0("RData/AUC_BS_2_",resname0,".csv"), append = TRUE, row.names = TRUE, col.names = FALSE, sep = ',')
write.table(BS.true.matrix, file = paste0("RData/AUC_BS_2_",resname0,".csv"), append = TRUE, row.names = TRUE, col.names = FALSE, sep = ',')
write.table(BS.matrix, file = paste0("RData/AUC_BS_2_",resname0,".csv"), append = TRUE, row.names = TRUE, col.names = FALSE, sep = ',')
resname0 <- "400_4000_3000_100_cph11_nc"
write.table(integrated_AUC_BS, file = paste0("RData/integrated_AUC_BS_2_",resname0,".csv"),
            append = F, row.names = TRUE, col.names = TRUE, sep = ',')
