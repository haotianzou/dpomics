
source('source_code/functions.R')
source('source_code/data_gen.R')
source('source_code/data_gen_negbin.R')
source('source_code/sFPCA.R')
load('data_gen_x1.RData')
#tp$gamma_z = seq(0.3,0.7, length.out=15)
#tp$n_omics_factors = 15
#save(tp, file = "data_gen.RData")

library(reticulate)
library(MOFA2)
library(MOFAdata)
library(readxl)
library(foreach)
#library(doParallel)
#cl <- makeCluster(5)
#registerDoParallel(cl)

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

start.time <- Sys.time()

result <- list()


for (seed in 1:110){
  
  
  set.seed(seed*2024)
  ############
  omics.all <- data_gen_omics_negbin(
    n_views = n_views, 
    n_samples = n_samples, 
    n_features = n_features, 
    n_factors = n_factors, 
    likelihood = c('gaussian', 'negbin', 'bernoulli')
  )
  
  omics = omics.all[[1]]
  Z_true = omics.all$Z[[1]]
  
  #lapply(data_omics,dim)
  
  #omics <- data_omics; omics.train <- omics[omics$ID %in% ID.train, ]; omics.test <- omics[omics$ID %in% ID.test, ]
  
  MOFAobject <- create_mofa(omics)
  #plot_data_overview(MOFAobject)
  
  ##############define options
  ###define data options
  data_opts <- get_default_data_options(MOFAobject)
  #head(data_opts)
  
  ## define model options
  model_opts <- get_default_model_options(MOFAobject)
  model_opts$num_factors = n_factors
  model_opts$likelihoods = c('gaussian', 'poisson', 'bernoulli')
  names(model_opts$likelihoods) = c('view_1', 'view_2', 'view_3')
  #head(model_opts)
  
  ## define train options
  train_opts <- get_default_training_options(MOFAobject)
  train_opts$seed = seed*2024
  #head(train_opts)
  
  #################Build and train the MOFA object
  MOFAobject <- prepare_mofa(
    object = MOFAobject,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
  )
  
  outfile = file.path(paste0("MOFA_fit/R4S7/", seed,".hdf5"))
  #outfile = file.path(paste0("./MOFA_fit_2/1241_4000_3000_200_cph11_nc_", seed,".hdf5"))
  MOFAobject.trained <- run_mofa(MOFAobject, outfile,use_basilisk = F) ## F
  MOFAobject.trained <- load_model(outfile,remove_inactive_factors=F)
  #saveRDS(MOFAobject.trained,paste0('MOFA_fit_2/',"seed",seed,"_MOFAobject.trained.rds"))
  
  #MOFAobject.trained <- readRDS("MOFAobject.trained.rds")
  ######
  # slotNames(MOFAobject.trained)
  # names(MOFAobject.trained@data)
  # dim(MOFAobject.trained@data$Methylation$group1)
  # names(MOFAobject.trained@expectations)
  # # Dimensionality of the factor matrix: 200 samples, 15 factors
  # dim(MOFAobject.trained@expectations$Z$group1)
  # # Dimensionality of the mRNA Weight matrix: 5000 features, 15 factors
  # dim(MOFAobject.trained@expectations$W$view_1)
  # dim(MOFAobject.trained@expectations$W$view_2)
  
  ###matrix
  Z_hat <- matrix(unlist(MOFAobject.trained@expectations$Z$group1),n_samples,n_factors,byrow=F)
  
  fname = paste0('MOFA_fit/R4S7/Z_', seed, '.RData')
  save(list = c('Z_hat', 'Z_true'), file = fname)
  
  cat(sprintf("Seed = %d\n", seed))
}
