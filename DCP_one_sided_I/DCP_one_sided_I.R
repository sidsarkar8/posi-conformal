###### CLOUD ######
######## plots for dcp
######## One-sided intervalss
### functions now return a correct value of alpha to take directly quantile with

set.seed(2024)

source("/home/siddhaas/Conf_redone/conf_I_delta/dumbgun_wellner_one_side/get_bds_dw_one_side.R")
source("/home/siddhaas/Conf_redone/conf_I_delta/dkw_one_side/get_bds_dkw_one_side.R")

library(dplyr)
library(quantregForest)
library(foreach)
library(doParallel)
#####################

data_gen = function(n)
{
  x = runif(n, min = 1, max = 5)
  y = numeric(n)
  
  eps1 = rnorm(n)
  eps2 = rnorm(n)
  unif_ind = as.numeric( runif(n) < 0.01 )
  y = rpois(n, lambda = (sin(x))^2 + 0.1) + 0.03*x*eps1 + 1*unif_ind*eps2
  
  return(data.frame(x = x, y = y))
}

#####################

data_gen_x = function(x)
{
  n = length(x)
  y = numeric(n)
  
  eps1 = rnorm(n)
  eps2 = rnorm(n)
  unif_ind = as.numeric( runif(n) < 0.01 )
  y = rpois(n, lambda = (sin(x))^2 + 0.1) + 0.03*x*eps1 + 1*unif_ind*eps2
  
  return(data.frame(x = x, y = y))
}

#####################

cond_quant_data_gen = function(x, n_runs = 100000, alpha )
{
  quants = data.frame( x = rep(x, each = length(alpha)), alpha = rep(alpha, length(x)),
                       q_low = rep(NA,length(alpha)*length(x)), 
                       q_up = rep(NA,length(alpha)*length(x)) )
  
  for( i in 1:length(x))
  {
    eps1 = rnorm(n_runs)
    eps2 = rnorm(n_runs)
    unif_ind = as.numeric( runif(n_runs) < 0.01 )
    y_i = rpois(n_runs, lambda = (sin(x[i]))^2 + 0.1) + 0.03*x[i]*eps1 + 1*unif_ind*eps2
    quants[ ((i-1)*length(alpha) + 1):(i*length(alpha)), 3] = quantile(y_i, alpha/2)
    quants[ ((i-1)*length(alpha) + 1):(i*length(alpha)), 4] = quantile(y_i, 1-(alpha/2))
  }
  return(quants)
}
#####################

n_train = 50000
n_calib = 10000
M = 500
n_test = 20000

delta_n = 0.1

#source("/home/siddhaarth/DCP_calib_fin/get_bds_dw.R")
#source("/home/siddhaarth/DCP_calib_fin/get_bds_dkw.R")

alpha_preset = c(0.5,1)

alpha_seq = start(from = 0.51, to = 0.9, length.out = 19)

# dw_quant_val1 = dw_quant_one_side(delta = delta_n,
#                                   I = alpha_preset,
#                                   n = n_calib,
#                                   n_runs = 500,
#                                   type = "Interval")
# 
# dkw_quant_val1 = dkw_quant_one_side(delta = delta_n,
#                                     I = alpha_preset,
#                                     n = n_calib,
#                                     n_runs = 500,
#                                     type = "Interval")

dw_quant_val1 = read.csv("/home/siddhaas/Conf_redone/conf_I_delta/dumbgun_wellner_one_side/dw_quant_10000_I.csv")[,2]
dkw_quant_val1 = read.csv("/home/siddhaas/Conf_redone/conf_I_delta/dkw_one_side/dkw_quant_10000_I.csv")[,2]

##### score computation
########################################################
########## coverage and width
########################################################

##################################################################
##################################################################
## baseline width and coverage

data_test_base = data_gen(n_test*5)
alpha_cov_base = numeric(19)
alpha_width_base = numeric(19)

print("conditional quant data gen for baseline dataset")
pred_test_base_full =  cond_quant_data_gen( data_test_base$x, 
                                            alpha = alpha_seq, n_runs = 10000 )

for( idx in 1:19 )
{
  alpha_n = alpha_seq[idx]
  
  pred_test_base = pred_test_base_full %>% filter(alpha == alpha_n)
  
  cov_base = data_test_base$y < pred_test_base$q_low | data_test_base$y > pred_test_base$q_up
  alpha_cov_base[idx] = 1- mean(cov_base)
  width_base = pred_test_base$q_up - pred_test_base$q_low
  alpha_width_base[idx] = mean(width_base)
}

##############################################################
##############################################################

data_test = data_gen(n_test)

parallel_fin_func = function(simul_idx)
{
  set.seed(simul_idx)
  
  alpha_cov_overfit = numeric(19)
  alpha_cov_calib = numeric(19)
  alpha_cov_calib_dw = numeric(19)
  alpha_cov_calib_dkw = numeric(19)
  
  alpha_width_overfit = numeric(19)
  alpha_width_calib = numeric(19)
  alpha_width_calib_dw = numeric(19)
  alpha_width_calib_dkw = numeric(19)
  
  data_train = data_gen(n_train)
  data_calib = data_gen(n_calib)
  
  print(paste("training for simul_idx:", simul_idx))  
  multi_rqfit <- quantregForest( x = data_train[,1] %>% as.matrix(ncol = 1 ), 
                                 y = data_train$y,
                                 nodesize = 300)
  
  ################# quant on data calib
  print(paste("predicting on data_calib for simul_idx:", simul_idx))  
  quant_ecdf = predict(multi_rqfit, newdata = data_calib[,1] %>%
                         as.matrix(ncol = 1 ), 
                       what = ecdf)
  score_n = numeric(n_calib)
  
  for(i in 1:n_calib)
  {
    score_n[i] = abs(quant_ecdf[[i]](data_calib$y[i]) - 0.5)
  }
  
  ################# quant on data train  
  print(paste("predicting on data_train for simul_idx:", simul_idx))  
  quant_ecdf_train = predict(multi_rqfit, newdata = data_train[,1] %>% as.matrix(ncol = 1 ), 
                             what = ecdf)
  score_n_overfit = numeric(n_train)
  
  for(i in 1:n_train)
  {
    score_n_overfit[i] = abs(quant_ecdf_train[[i]](data_train$y[i]) - 0.5)
  }
  
  ##### cutoff for CDFs
  
  print("cutoffs and test data")
  
  ##### computing cutoffs
  for(idx in 1:19)
  {  
    #print(idx)
    alpha_n = alpha_seq[idx]
    
    ### Duembgun Wellner
    ### which point will give us the lower bound of F(x) >= 1- alpha
    
    alpha_corr_dw = alpha_dw_one_side(alpha = alpha_n,
                                      I = alpha_preset,
                                      n = n_calib,
                                      dw_quant_val = dw_quant_val1,
                                      type = "Interval")
    
    s_cutoff_dw = quantile(score_n, probs = 1 - alpha_corr_dw)
    
    ### DKW simple
    ### which point will give us the lower bound of F(x) >= 1- alpha
    alpha_corr_dkw = alpha_dkw_one_side(alpha = alpha_n,
                                        I = alpha_preset,
                                        n = n_calib,
                                        dkw_quant_val = dkw_quant_val1,
                                        type = "Interval")
    
    s_cutoff_dkw = quantile(score_n, probs = 1 - alpha_corr_dkw)
    
    ### standard split conf
    s_cutoff = quantile(score_n, probs = 1- alpha_n)
    
    ### cutoff based on training data
    s_cutoff_overfit = quantile(score_n_overfit, probs = 1 - alpha_n)
    
    ######### 
    
    print(paste("predicting on data_calib for simul_idx:", simul_idx, ", alpha:", alpha_n))  
    pred_overfit = predict(multi_rqfit, newdata = data_test[,1] %>% as.matrix(nrow = 1),
                           what = c( 0.5 - s_cutoff_overfit  , 0.5 + s_cutoff_overfit))
    pred_calib = predict(multi_rqfit, newdata = data_test[,1] %>% as.matrix(nrow = 1),
                         what = c( 0.5 - s_cutoff, 0.5 + s_cutoff ))
    pred_calib_dw = predict(multi_rqfit, newdata = data_test[,1] %>% as.matrix(nrow = 1),
                            what = c( 0.5 - s_cutoff_dw, 0.5 + s_cutoff_dw ))
    pred_calib_dkw = predict(multi_rqfit, newdata = data_test[,1] %>% as.matrix(nrow = 1),
                             what = c( 0.5 - s_cutoff_dkw, 0.5 + s_cutoff_dkw ))
    
    cov_overfit = data_test$y < pred_overfit[,1] | data_test$y > pred_overfit[,2]
    cov_calib = data_test$y < pred_calib[,1]  | data_test$y > pred_calib[,2]
    cov_calib_dw = data_test$y < pred_calib_dw[,1]  | data_test$y > pred_calib_dw[,2]
    cov_calib_dkw = data_test$y < pred_calib_dkw[,1]  | data_test$y > pred_calib_dkw[,2]
    
    alpha_cov_overfit[idx] = 1- mean(cov_overfit)
    alpha_cov_calib[idx] = 1- mean(cov_calib)
    alpha_cov_calib_dw[idx] = 1- mean(cov_calib_dw)
    alpha_cov_calib_dkw[idx] = 1- mean(cov_calib_dkw)
    
    width_overfit = pred_overfit[,2] - pred_overfit[,1]
    width_calib = pred_calib[,2] - pred_calib[,1] 
    width_calib_dw = pred_calib_dw[,2] - pred_calib_dw[,1]
    width_calib_dkw = pred_calib_dkw[,2] - pred_calib_dkw[,1]
    
    alpha_width_overfit[idx] = mean(width_overfit)
    alpha_width_calib[idx] = mean(width_calib)
    alpha_width_calib_dw[idx] = mean(width_calib_dw)
    alpha_width_calib_dkw[idx] = mean(width_calib_dkw)
  }
  
  df_dcp = data.frame( alpha = rep(alpha_seq,  4), 
                       Method = rep(c("training data quantile", "split dcp","DKW univ", 
                                      "DW univ"), each = 19),
                       coverage = rep(NA, 4*19),
                       width = rep(NA, 4*19),
                       simul_id = rep(NA, 4*19),
                       width_base = rep(NA, 4*19) )
  
  
  df_dcp$coverage[ (1:19)] = alpha_cov_overfit
  df_dcp$coverage[ (1:19) + (19)] = alpha_cov_calib  
  df_dcp$coverage[ (1:19) + (38)] = alpha_cov_calib_dkw
  df_dcp$coverage[ (1:19) + (57)] = alpha_cov_calib_dw
  
  df_dcp$width[ (1:19)] = alpha_width_overfit
  df_dcp$width[ (1:19) + (19)] = alpha_width_calib  
  df_dcp$width[ (1:19) + (38)] = alpha_width_calib_dkw
  df_dcp$width[ (1:19) + (57)] = alpha_width_calib_dw
  
  df_dcp$simul_id[ (1:19)] = simul_idx
  df_dcp$simul_id[ (1:19) + (19)] = simul_idx 
  df_dcp$simul_id[ (1:19) + (38)] = simul_idx
  df_dcp$simul_id[ (1:19) + (57)] = simul_idx
  
  df_dcp$width_base[ (1:19)] = alpha_width_base 
  df_dcp$width_base[ (1:19) + (19)] = alpha_width_base
  df_dcp$width_base[ (1:19) + (38)] = alpha_width_base
  df_dcp$width_base[ (1:19) + (57)] = alpha_width_base
  
  df_dcp = df_dcp %>% mutate(diff_coverage = coverage + alpha - 1)
  df_dcp = df_dcp %>% mutate( ratio_width = width/width_base)
  print(paste("Done with", simul_idx))  
  return(df_dcp)
}

##################################################################
##################################################################

numCores = detectCores()
registerDoParallel(numCores-1) 

print("run sims")
df_dcp_fin = foreach( simul_idx = 1:M,  .combine = rbind ) %dopar%{
  parallel_fin_func(simul_idx)
}


write.csv(df_dcp_fin, file = "/home/siddhaas/Conf_redone/DCP_one_sided_I/dcp_files/df_dcp_I.csv")

stopImplicitCluster()