load('data_gen.RData')

tp$N.train = 300
tp$N.test = 100

f1 = function(t) 10 + 5*t + t^2
f2 = function(t) 10 + sqrt(t)
f3 = function(t) sin(2*t) + cos(2*t) + log(t+1)

tp$mu[[1]] = f1(tp$tnew)
tp$mu[[2]] = f2(tp$tnew)
tp$mu[[3]] = f3(tp$tnew)

tp$d0 = c(5, 2)
tp$d1 = c(5)

tp$alpha = NULL

tp$gamma0 = c(0.4, 0.3)
tp$gamma_z = c(0.5, 0.5, 1, 1, 1)

tp$n_omics_views = 3
tp$n_omics_factors = 5
tp$n_omics_features = 1000

save(tp, file = 'data_gen_x1.RData')

##
load('data_gen_x1.RData')
tp$n_omics_features = 5000
save(tp, file = 'data_gen_x2.RData')


##
load('data_gen.RData')

tp$N.train = 100
tp$N.test = 100

f1 = function(t) 10 + 5*t + t^2
f2 = function(t) 10 + sqrt(t)
f3 = function(t) sin(2*t) + cos(2*t) + log(t+1)

tp$mu[[1]] = f1(tp$tnew)
tp$mu[[2]] = f2(tp$tnew)
tp$mu[[3]] = f3(tp$tnew)

tp$d0 = c(5, 2)
tp$d1 = c(5)

tp$alpha = NULL

tp$gamma0 = c(0.4, 0.3)
tp$gamma_z = c(0.5, 0.5, 1, 1, 1)

tp$n_omics_views = 3
tp$n_omics_factors = 5
tp$n_omics_features = 1000

save(tp, file = 'data_gen_R4S1.RData')

##
load('data_gen.RData')

tp$N.train = 200
tp$N.test = 100

f1 = function(t) 10 + 5*t + t^2
f2 = function(t) 10 + sqrt(t)
f3 = function(t) sin(2*t) + cos(2*t) + log(t+1)

tp$mu[[1]] = f1(tp$tnew)
tp$mu[[2]] = f2(tp$tnew)
tp$mu[[3]] = f3(tp$tnew)

tp$d0 = c(5, 2)
tp$d1 = c(5)

tp$alpha = NULL

tp$gamma0 = c(0.4, 0.3)
tp$gamma_z = c(0.5, 0.5, 1, 1, 1)

tp$n_omics_views = 3
tp$n_omics_factors = 5
tp$n_omics_features = 1000

save(tp, file = 'data_gen_R4S2.RData')



##
load('data_gen.RData')

tp$N.train = 300
tp$N.test = 100

f1 = function(t) 10 + 5*t + t^2
f2 = function(t) 10 + sqrt(t)
f3 = function(t) sin(2*t) + cos(2*t) + log(t+1)

tp$mu[[1]] = f1(tp$tnew)
tp$mu[[2]] = f2(tp$tnew)
tp$mu[[3]] = f3(tp$tnew)

tp$d0 = c(5, 2)
tp$d1 = c(5)

tp$alpha = NULL

tp$gamma0 = c(0.4, 0.3)
tp$gamma_z = c(0.5, 0.5, 1, 1, 1)

tp$n_omics_views = 3
tp$n_omics_factors = 5
tp$n_omics_features = 1000

tp$logh0 = -2

save(tp, file = 'data_gen_R4S3.RData')

##
load('data_gen.RData')

tp$N.train = 300
tp$N.test = 100

f1 = function(t) 10 + 5*t + t^2
f2 = function(t) 10 + sqrt(t)
f3 = function(t) sin(2*t) + cos(2*t) + log(t+1)

tp$mu[[1]] = f1(tp$tnew)
tp$mu[[2]] = f2(tp$tnew)
tp$mu[[3]] = f3(tp$tnew)

tp$d0 = c(5, 2)
tp$d1 = c(5)

tp$alpha = NULL

tp$gamma0 = c(0.4, 0.3)
tp$gamma_z = c(0.5, 0.5, 1, 1, 1)

tp$n_omics_views = 3
tp$n_omics_factors = 5
tp$n_omics_features = 1000

tp$logh0 = -1

save(tp, file = 'data_gen_R4S4.RData')

##
load('data_gen.RData')

tp$N.train = 300
tp$N.test = 100

f1 = function(t) 10 + 5*t + t^2
f2 = function(t) 10 + sqrt(t)
f3 = function(t) sin(2*t) + cos(2*t) + log(t+1)

tp$mu[[1]] = f1(tp$tnew)
tp$mu[[2]] = f2(tp$tnew)
tp$mu[[3]] = f3(tp$tnew)

tp$d0 = c(5, 2)
tp$d1 = c(5)

tp$alpha = NULL

tp$gamma0 = c(0.4, 0.3)
tp$gamma_z = c(0.5, 0.5, 1, 1, 1)

tp$n_omics_views = 3
tp$n_omics_factors = 5
tp$n_omics_features = 1000

tp$beta = c(1.0, -0.2, 0.2)

save(tp, file = 'data_gen_R4S5.RData')

##
load('data_gen.RData')

tp$N.train = 300
tp$N.test = 100

f1 = function(t) 10 + 5*t + t^2
f2 = function(t) 10 + sqrt(t)
f3 = function(t) sin(2*t) + cos(2*t) + log(t+1)

tp$mu[[1]] = f1(tp$tnew)
tp$mu[[2]] = f2(tp$tnew)
tp$mu[[3]] = f3(tp$tnew)

tp$d0 = c(5, 2)
tp$d1 = c(5)

tp$alpha = NULL

tp$gamma0 = c(0.4, 0.3)
tp$gamma_z = c(0.5, 0.5, 1, 1, 1)

tp$n_omics_views = 3
tp$n_omics_factors = 5
tp$n_omics_features = 1000

tp$beta = c(1.0, -2.0, 2.0)

save(tp, file = 'data_gen_R4S6.RData')
