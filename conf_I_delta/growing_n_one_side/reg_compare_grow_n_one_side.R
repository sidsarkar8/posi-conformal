#### larger n

source("/home/siddhaas/conf_I_delta/dumbgun_wellner_one_side/get_bds_dw_one_side.R")
source("/home/siddhaas/conf_I_delta/dkw_one_side/get_bds_dkw_one_side.R")
source("/home/siddhaas/conf_I_delta/two_stage/dumbgun_wellner_new/get_bds_dw_int1.R")
source("/home/siddhaas/conf_I_delta/two_stage/dkw/get_bds_dkw_int.R")

set.seed(2024)

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

n_calib_seq = c(500,1000,2000,5000,10000)
delta_n = 0.1
alpha_preset = c(0.02,0.5)
alpha_fin = seq(from = alpha_preset[1], to = alpha_preset[2], by = 0.01) 


##################################################################
##################################################################
## baseline width and coverage

print("baseline")
data_test_base = data_gen(10000)
alpha_cov_base = numeric(length(alpha_fin))
alpha_width_base = numeric(length(alpha_fin))

pred_test_base_full =  cond_quant_data_gen( data_test_base$x, alpha = alpha_fin, n_runs = 10000 )
for( idx in 1:length(alpha_fin) )
{
  print(idx)
  alpha_n = alpha_fin[idx]
  
  pred_test_base = pred_test_base_full %>% filter(alpha == alpha_n)
  
  cov_base = data_test_base$y < pred_test_base$q_low | data_test_base$y > pred_test_base$q_up
  alpha_cov_base[idx] = 1- mean(cov_base)
  width_base = pred_test_base$q_up - pred_test_base$q_low
  alpha_width_base[idx] = mean(width_base)
}


######################################

for( n_idx in 1:length(n_calib_seq))
{
  n_train = 10000
  n_calib = n_calib_seq[n_idx]
  M = 100
  n_test = 20000
  
  
  #####################
  
  print("generating bounds")
  
  ########## 2 stage procedures
  ########## Dumben wellner
  
  bds_dw_full = get_bds_dw_int(alpha = delta_n, n = n_calib,
                               I = 1:n_calib) 
  
  I_dw_1 = which(bds_dw_full$low >= 1- alpha_preset[2] &
                   bds_dw_full$low <= 1- alpha_preset[1])  
  
  
  bds_dw_int = get_bds_dw_int(alpha = delta_n, n = n_calib, I = I_dw_1) 
  
  alpha_corr_dw_2stage = function(al)
  {
    if(length(which(bds_dw_int$low >= 1- al))!= 0 )
    {
      low_idx_dw = min( which(bds_dw_int$low >= 1- al))
    }else{
      low_idx_dw = n_calib
    }
    
    return(1 - low_idx_dw/n_calib)
  }
  
  ##### DKW ##############
  
  bds_dkw_full = get_bds_dkw_simp_int(alpha = delta_n, n = n_calib, I = 1:n_calib) 
  
  I_dkw_1 = which(bds_dkw_full$low >= 1- alpha_preset[2] &
                    bds_dkw_full$low <= 1- alpha_preset[1])  
  
  bds_dkw_int = get_bds_dkw_simp_int(alpha = delta_n, n = n_calib, I = I_dkw_1) 
  
  
  alpha_corr_dkw_2stage = function(al)
  {
    if(length(which(bds_dkw_int$low >= 1- al))!= 0 )
    {
      low_idx_dkw = min( which(bds_dkw_int$low >= 1- al))
    }else{
      low_idx_dkw = n_calib
    }
    return(1 - low_idx_dkw/n_calib)
  }
  
  ############ new I,delta methods #####
  
  dw_quant_val1 = dw_quant_one_side(delta = delta_n,
                                       I = alpha_preset,
                                       n = n_calib,
                                       n_runs = 50000,
                                       type = "Interval")
  
  dkw_quant_val1 = dkw_quant_one_side(delta = delta_n,
                                         I = alpha_preset,
                                         n = n_calib,
                                         n_runs = 50000,
                                         type = "Interval")
  
  
  ####################################################################
  ####################################################################
  
  ########## coverage and width
  ########################################################
  
  data_test = data_gen(n_test)
  
  
  parallel_fin_func = function(simul_idx)
  {
    set.seed(simul_idx)
    
    alpha_cov_overfit =rep(NA,length(alpha_fin))
    alpha_cov_calib =rep(NA,length(alpha_fin))
    alpha_cov_calib_dw =rep(NA,length(alpha_fin))
    alpha_cov_calib_dkw =rep(NA,length(alpha_fin))
    alpha_cov_calib_dkw_new =rep(NA,length(alpha_fin))
    alpha_cov_calib_dw_new =rep(NA,length(alpha_fin))
    
    alpha_width_overfit =rep(NA,length(alpha_fin))
    alpha_width_calib =rep(NA,length(alpha_fin))
    alpha_width_calib_dw =rep(NA,length(alpha_fin))
    alpha_width_calib_dkw =rep(NA,length(alpha_fin))
    alpha_width_calib_dkw_new =rep(NA,length(alpha_fin))
    alpha_width_calib_dw_new =rep(NA,length(alpha_fin))
    
    
    ##### gen data
    
    data_train = data_gen(n_train)
    data_calib = data_gen(n_calib)
    
    multi_rqfit <- quantregForest( x = data_train[,1] %>% as.matrix(ncol = 1 ), 
                                   y = data_train$y,
                                   nodesize = 300)
    
    ################# quant on data calib
    quant_ecdf = predict(multi_rqfit, newdata = data_calib[,1] %>% as.matrix(ncol = 1 ), 
                         what = ecdf)
    score_n = numeric(n_calib)
    
    for(i in 1:n_calib)
    {
      score_n[i] = abs(quant_ecdf[[i]](data_calib$y[i]) - 0.5)
    }
    
    ################# quant on data train  
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
    for(idx in 1:length(alpha_fin))
    {  
      #print(idx)
      alpha_n = alpha_fin[idx]
      
      ### Duembgun Wellner
      ### which point will give us the lower bound of F(x) >= 1- alpha
      
      s_cutoff_dw = quantile(score_n, probs = 1 - alpha_corr_dw_2stage(alpha_n))
      
      ### DKW simple
      ### which point will give us the lower bound of F(x) >= 1- alpha
      
      s_cutoff_dkw = quantile(score_n, probs = 1 - alpha_corr_dkw_2stage(alpha_n))
      
      ### DW new
      s_cutoff_dw_new = quantile(score_n, 
                                 probs = 1 - alpha_dw_one_side(alpha = alpha_n,
                                                          I = alpha_preset,
                                                          n = n_calib,
                                                          dw_quant_val = dw_quant_val1,
                                                          type = "Interval"))
      
      ### DKW new
      s_cutoff_dkw_new = quantile(score_n, 
                                  probs = 1 - alpha_dkw_one_side(alpha = alpha_n,
                                                            I = alpha_preset,
                                                            n = n_calib,
                                                            dkw_quant_val = dkw_quant_val1,
                                                            type = "Interval"))
      
      ### standard split conf PAC
      #s_cutoff = quantile(score_n, probs = 1- alpha_n)
      s_cutoff = quantile(score_n, 
                          probs = qbinom(1-delta_n, size = n_calib, 
                                         prob = 1- alpha_n)/n_calib)
      
      ### cutoff based on training data
      #s_cutoff_overfit = quantile(score_n_overfit, probs = 1 - alpha_n)
      s_cutoff_overfit = quantile(score_n_overfit, 
                                  probs = qbinom(1-delta_n, size = n_calib, 
                                                 prob = 1- alpha_n)/n_calib)
      
      ######### 
      
      pred_overfit = predict(multi_rqfit, newdata = data_test[,1] %>% as.matrix(nrow = 1),
                             what = c( 0.5 - s_cutoff_overfit  , 0.5 + s_cutoff_overfit))
      pred_calib = predict(multi_rqfit, newdata = data_test[,1] %>% as.matrix(nrow = 1),
                           what = c( 0.5 - s_cutoff, 0.5 + s_cutoff ))
      pred_calib_dw = predict(multi_rqfit, newdata = data_test[,1] %>% as.matrix(nrow = 1),
                              what = c( 0.5 - s_cutoff_dw, 0.5 + s_cutoff_dw ))
      pred_calib_dkw = predict(multi_rqfit, newdata = data_test[,1] %>% as.matrix(nrow = 1),
                               what = c( 0.5 - s_cutoff_dkw, 0.5 + s_cutoff_dkw ))
      pred_calib_dkw_new = predict(multi_rqfit, newdata = data_test[,1] %>% as.matrix(nrow = 1),
                                   what = c( 0.5 - s_cutoff_dkw_new, 0.5 + s_cutoff_dkw_new ))
      pred_calib_dw_new = predict(multi_rqfit, newdata = data_test[,1] %>% as.matrix(nrow = 1),
                                  what = c( 0.5 - s_cutoff_dw_new, 0.5 + s_cutoff_dw_new ))
      
      cov_overfit = data_test$y < pred_overfit[,1] | data_test$y > pred_overfit[,2]
      cov_calib = data_test$y < pred_calib[,1]  | data_test$y > pred_calib[,2]
      cov_calib_dw = data_test$y < pred_calib_dw[,1]  | data_test$y > pred_calib_dw[,2]
      cov_calib_dkw = data_test$y < pred_calib_dkw[,1]  | data_test$y > pred_calib_dkw[,2]
      cov_calib_dkw_new = data_test$y < pred_calib_dkw_new[,1]| data_test$y > pred_calib_dkw_new[,2]
      cov_calib_dw_new = data_test$y < pred_calib_dw_new[,1]| data_test$y > pred_calib_dw_new[,2]
      
      alpha_cov_overfit[idx] = 1- mean(cov_overfit)
      alpha_cov_calib[idx] = 1- mean(cov_calib)
      alpha_cov_calib_dw[idx] = 1- mean(cov_calib_dw)
      alpha_cov_calib_dkw[idx] = 1- mean(cov_calib_dkw)
      alpha_cov_calib_dkw_new[idx] = 1- mean(cov_calib_dkw_new)
      alpha_cov_calib_dw_new[idx] = 1- mean(cov_calib_dw_new)
      
      width_overfit = pred_overfit[,2] - pred_overfit[,1]
      width_calib = pred_calib[,2] - pred_calib[,1] 
      width_calib_dw = pred_calib_dw[,2] - pred_calib_dw[,1]
      width_calib_dkw = pred_calib_dkw[,2] - pred_calib_dkw[,1]
      width_calib_dkw_new = pred_calib_dkw_new[,2] - pred_calib_dkw_new[,1]
      width_calib_dw_new = pred_calib_dw_new[,2] - pred_calib_dw_new[,1]
      
      alpha_width_overfit[idx] = mean(width_overfit)
      alpha_width_calib[idx] = mean(width_calib)
      alpha_width_calib_dw[idx] = mean(width_calib_dw)
      alpha_width_calib_dkw[idx] = mean(width_calib_dkw)
      alpha_width_calib_dkw_new[idx] = mean(width_calib_dkw_new)
      alpha_width_calib_dw_new[idx] = mean(width_calib_dw_new)
      
    }
    
    ####
    df_dcp = data.frame( alpha = rep(alpha_fin,  6), 
                         Method = rep(c("training data quantile", "split dcp","DKW", "DW", "DKW new", "DW new"), each = length(alpha_fin)),
                         coverage = rep(NA, 6*length(alpha_fin)),
                         width = rep(NA, 6*length(alpha_fin)),
                         simul_id = rep(NA, 6*length(alpha_fin)),
                         width_base = rep(NA, 6*length(alpha_fin)) )
    
    df_dcp$coverage[ (1 : length(alpha_fin))] = alpha_cov_overfit
    df_dcp$coverage[ (1 : length(alpha_fin)) + length(alpha_fin)] = alpha_cov_calib 
    df_dcp$coverage[ (1 : length(alpha_fin)) + (2*length(alpha_fin))] = alpha_cov_calib_dkw
    df_dcp$coverage[ (1 : length(alpha_fin)) + (3*length(alpha_fin))] = alpha_cov_calib_dw
    df_dcp$coverage[ (1 : length(alpha_fin)) + (4*length(alpha_fin))] = alpha_cov_calib_dkw_new
    df_dcp$coverage[ (1 : length(alpha_fin)) + (5*length(alpha_fin))] = alpha_cov_calib_dw_new
    
    df_dcp$width[ (1 : length(alpha_fin))] = alpha_width_overfit
    df_dcp$width[ (1 : length(alpha_fin)) + length(alpha_fin)] = alpha_width_calib 
    df_dcp$width[ (1 : length(alpha_fin)) + (2*length(alpha_fin))] = alpha_width_calib_dkw
    df_dcp$width[ (1 : length(alpha_fin)) + (3*length(alpha_fin))] = alpha_width_calib_dw
    df_dcp$width[ (1 : length(alpha_fin)) + (4*length(alpha_fin))] = alpha_width_calib_dkw_new
    df_dcp$width[ (1 : length(alpha_fin)) + (5*length(alpha_fin))] = alpha_width_calib_dw_new
    
    df_dcp$simul_id[ (1 : length(alpha_fin))] = simul_idx
    df_dcp$simul_id[ (1 : length(alpha_fin)) + length(alpha_fin)] = simul_idx 
    df_dcp$simul_id[ (1 : length(alpha_fin)) + (2*length(alpha_fin))] = simul_idx
    df_dcp$simul_id[ (1 : length(alpha_fin)) + (3*length(alpha_fin))] = simul_idx
    
    df_dcp$width_base[ (1 : length(alpha_fin))] = alpha_width_base 
    df_dcp$width_base[ (1 : length(alpha_fin)) + length(alpha_fin)] = alpha_width_base
    df_dcp$width_base[ (1 : length(alpha_fin)) + (2*length(alpha_fin))] = alpha_width_base
    df_dcp$width_base[ (1 : length(alpha_fin)) + (3*length(alpha_fin))] = alpha_width_base
    df_dcp$width_base[ (1 : length(alpha_fin)) + (4*length(alpha_fin))] = alpha_width_base
    df_dcp$width_base[ (1 : length(alpha_fin)) + (5*length(alpha_fin))] = alpha_width_base
    
    
    
    df_dcp = df_dcp %>% mutate(diff_coverage = coverage + alpha - 1)
    df_dcp = df_dcp %>% mutate( ratio_width = width/width_base)
    
    print(paste("Done with", simul_idx))  
    return(df_dcp)
  }
  
  
  numCores = detectCores()
  registerDoParallel(numCores) 
  
  print("run sims")
  df_dcp_fin = foreach( simul_idx = 1:M,  .combine = rbind ) %dopar%{
    parallel_fin_func(simul_idx)
  }
  
  
  write.csv(df_dcp_fin, 
            file = paste("/home/siddhaas/conf_I_delta/growing_n_one_side/df_dcp_I_delta_",n_calib,".csv",sep = "" ))
}
stopImplicitCluster()
####################################################################


####################################################################
####################################################################

#write.csv(df_prop_cross, file = "/home/siddhaarth/DCP_calib_fin/prop_files2/df_prop_0.1_1.csv")
#write.csv(df_dcp, file = "/home/siddhaarth/DCP_calib_fin/dcp_files2/df_dcp_1.csv")

# df_dcp %>% 
#   ggplot(mapping = aes(x = alpha, y = diff_coverage, color = Method)) +
#   geom_line(aes(group = simul_id),alpha = 0.3, 
#             show.legend = FALSE) +
#   facet_grid(~ Method) + 
#   scale_x_continuous(breaks = seq(1/20, length(alpha_fin)/20, by = 0.2)) +
#   geom_hline(yintercept = 0,linetype = 2) + 
#   ylab(expression(paste("Coverage - (1 -", alpha, ")", sep=""))) + 
#   ggtitle(expression(paste("Coverage of methods with ",
#                            delta," = 0.1" ,sep = "")),
#           subtitle = expression(paste("(difference with (1-",alpha,"))", sep = ""))) +
#   theme(text = element_text(size = 13)) 
# 
# ggsave("/home/siddhaarth/DCP_calib_fin/cov2/dcp_cov_1.pdf",
# .      width = 12, height = 3.5, dpi = 300, units = "in" )

####################################################################

# df_dcp %>% 
#   ggplot(mapping = aes(x = alpha, y = ratio_width, color = Method)) +
#   geom_line(aes(group = simul_id), alpha = 0.3, 
#             show.legend = FALSE) +
#   facet_grid(~ Method) + 
#   geom_hline(yintercept = 1,linetype = 2) + 
#   scale_x_continuous(breaks = seq(1/20, length(alpha_fin)/20, by = 0.2)) +
#   ylab("Width / Width(true quantiles)") + 
#   ggtitle(expression(paste("Width of prediction set across methods ",
#                            delta,"= 0.1", sep = "" )),
#           subtitle = "(ratio to true quantile widths)") +
#   theme(text = element_text(size = 13)) 
# 
# ggsave("/home/siddhaarth/DCP_calib_fin/width2/dcp_width_1.pdf",
#       width = 12, height = 3.5, dpi = 300, units = "in" )
