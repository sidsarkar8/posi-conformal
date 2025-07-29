##### adding bates confidence intervals

set.seed(2022)
library(dplyr)
library(ranger)
library(quantregForest)
library(tidyverse)
##########################################
##########################################

##### bounds generator
source("/Users/sidsarkar/Desktop/Conf/dkw/get_bds_dkw_int.R")
source("/Users/sidsarkar/Desktop/Conf/dumbgun_wellner_new/get_bds_dw_int1.R")

###################### Asymptotic Eicker ######################
log2 = function(x)
{return(log(log(x)))}
log3 = function(x)
{return(log(log(log(x))))}


###################### EIC with all indices ######################

get_bds_eic1 = function(alpha = 0.05, n)
{
  bds_eic = data.frame(low = rep(NA, n),
                       up = rep(NA, n))
  
  t_eic = log( -2/(sqrt(pi)*log(1-alpha)) )
  
  t_eic_n = (2*log2(n))^(1/2)*( 1 + ((log3(n) + 2*t_eic)/(4*log2(n))) )
  
  bds_eic$low = (1:n)/(n+1) -  t_eic_n*sqrt( (1:n)*(n + 1 - (1:n)))*(n+1)^(-1)*(n+2)^(-0.5)
  bds_eic$up = (1:n)/(n+1) + t_eic_n*sqrt( (1:n)*(n + 1 - (1:n)))*(n+1)^(-1)*(n+2)^(-0.5)
  
  bds_eic$low[ bds_eic$low < 0] = 0
  bds_eic$up[ bds_eic$up > 1] = 1
  
  return(bds_eic)
}

###################### Simes ######################

get_bds_simes = function(alpha = 0.05, n , k = NA)
{
  if(is.na(k)){k = n/2}
  #if(is.na(k_up)){k_up = n/2}
  
  bds_simes = data.frame(low = rep(NA, n),
                         up = rep(NA, n))
  
  for (i in 1:n)
  {
    #print(i)
    if(i < k){
      bds_simes$up[n - i + 1] = 1
      bds_simes$low[i] = 0
    }else{
      simes_seq = (i:(i - k + 1)) / (n:(n - k + 1))
      temp_log_simes = ( sum(log(simes_seq))/k + log(alpha/2)/k )
      
      #bds_simes$low[i] = (alpha/2)^(1/k)*(prod(simes_seq))^(1/k)
      bds_simes$low[i] = as.double(exp(temp_log_simes))  
      
      #bds_simes$up[n - i + 1] = 1 - (alpha / 2) ^ (1 / k) * (prod(simes_seq))^(1/k)
      bds_simes$up[n - i + 1] = 1 - as.double( exp(temp_log_simes))
    }
  }
  
  bds_simes$low[bds_simes$low<0] = 0
  bds_simes$up[bds_simes$up>1] = 1
  
  return(bds_simes)
}

##################################

#####################

airf = read.table("/Users/sidsarkar/Desktop/Conf/bike/airfoil_self_noise.dat")
colnames(airf)[6] = "y"
n_train = 700
n_calib = 300
n_test = 503


shuffle_idx = sample(nrow(airf),rep = F)
airf_train = airf[shuffle_idx[1:n_train],]
airf_calib = airf[shuffle_idx[(n_train + 1):(n_train + n_calib)],]
airf_test = airf[shuffle_idx[- (1:(n_train + n_calib))],]

delta_n = 0.1
alpha_preset = c(0.1,0.5)

print("generating bounds")

##################################

bds_dw_full = get_bds_dw_int(alpha = delta_n, n = n_calib, I = 1:n_calib) 

I_dw_1 = which(bds_dw_full$low >= 1- alpha_preset[2] &
                 bds_dw_full$low <= 1- alpha_preset[1])  

bds_dw_int = get_bds_dw_int(alpha = delta_n, n = n_calib, I = I_dw_1) 

bds_dkw_simp_full = get_bds_dkw_simp_int(alpha = delta_n, n = n_calib, I = 1:n_calib) 
I_dkw_simp_1 = which(bds_dkw_simp_full$low >= 1- alpha_preset[2] &
                       bds_dkw_simp_full$low <= 1- alpha_preset[1])  

bds_dkw_simp_int = get_bds_dkw_simp_int(alpha = delta_n, n = n_calib, I = I_dkw_simp_1) 

### new bates bounds
bds_simes = get_bds_simes(alpha = delta_n, n = n_calib) 
bds_eic = get_bds_eic1(alpha = delta_n, n = n_calib) 

##################################

airf_rqfit <- quantregForest( y = airf_train$y, x = dplyr::select(airf_train, -y))


################# quant on data calib
quant_ecdf = predict(airf_rqfit, newdata = dplyr::select(airf_calib, -y), 
                     what = ecdf)

score_n = numeric(n_calib)

for(i in 1:n_calib)
{
  score_n[i] = abs(quant_ecdf[[i]](airf_calib$y[i]) - 0.5)
}

################# quant on data train  
quant_ecdf_train = predict(airf_rqfit, newdata = dplyr::select(airf_train, -y) , 
                           what = ecdf)
score_n_overfit = numeric(n_train)

for(i in 1:n_train)
{
  score_n_overfit[i] = abs(quant_ecdf_train[[i]](airf_train$y[i]) - 0.5)
}

##### cutoff for CDFs
print("cutoffs and test data")
##### computing cutoffs

alpha_fin = seq(from = alpha_preset[1], to = alpha_preset[2], by = 0.01) 

alpha_cov_overfit = numeric(length(alpha_fin))
alpha_cov_calib = numeric(length(alpha_fin))
alpha_cov_calib_dw = numeric(length(alpha_fin))
alpha_cov_calib_dkw_simp = numeric(length(alpha_fin))
alpha_cov_calib_dkw_complex = numeric(length(alpha_fin))
alpha_cov_calib_simes = numeric(length(alpha_fin))
alpha_cov_calib_eic = numeric(length(alpha_fin))

alpha_width_overfit = numeric(length(alpha_fin))
alpha_width_calib = numeric(length(alpha_fin))
alpha_width_calib_dw = numeric(length(alpha_fin))
alpha_width_calib_dkw_simp = numeric(length(alpha_fin))
alpha_width_calib_dkw_complex = numeric(length(alpha_fin))
alpha_width_calib_simes = numeric(length(alpha_fin))
alpha_width_calib_eic = numeric(length(alpha_fin))

for(idx in 1:length(alpha_fin))
{  
  print(idx)
  alpha_n = alpha_fin[idx]
  
  ### Duembgun Wellner
  ### which point will give us the lower bound of F(x) >= 1- alpha
  low_idx_dw = min( which(bds_dw_int$low >= 1- alpha_n))
  F_s_cutoff_dw =  low_idx_dw/n_calib
  F_s_upperbd_dw = bds_dw_int$up[low_idx_dw] 
  
  s_cutoff_dw = quantile(score_n, probs = F_s_cutoff_dw)
  
  ### DKW simple
  ### which point will give us the lower bound of F(x) >= 1- alpha
  low_idx_dkw_simp = min( which(bds_dkw_simp_int$low >= 1- alpha_n))
  F_s_cutoff_dkw_simp =  low_idx_dkw_simp/n_calib
  F_s_upperbd_dkw_simp = bds_dkw_simp_int$up[low_idx_dkw_simp] 
  
  s_cutoff_dkw_simp = quantile(score_n, probs = F_s_cutoff_dkw_simp)
  
  ### standard split conf
  s_cutoff = quantile(score_n, probs = 1- alpha_n)
  
  ### cutoff based on training data
  s_cutoff_overfit = quantile(score_n_overfit, probs = 1 - alpha_n)
  
  ### Generalized Simes Inequality with k = n/2
  ### which point will give us the lower bound of F(x) >= 1- alpha
  low_idx_simes = min( which(bds_simes$low >= 1- alpha_n))
  F_s_cutoff_simes =  low_idx_simes/n_calib
  F_s_upperbd_simes = bds_simes$up[low_idx_simes] 
  
  s_cutoff_simes = quantile(score_n, probs = F_s_cutoff_simes)
  
  ### Eicker statistic
  ### which point will give us the lower bound of F(x) >= 1- alpha
  low_idx_eic = min( which(bds_eic$low >= 1- alpha_n))
  F_s_cutoff_eic =  low_idx_eic/n_calib
  F_s_upperbd_eic = bds_eic$up[low_idx_eic] 
  
  s_cutoff_eic = quantile(score_n, probs = F_s_cutoff_eic)
  ######### 
  
  pred_overfit = predict(airf_rqfit, newdata = dplyr::select(airf_test, -y) ,
                         what = c( 0.5 - s_cutoff_overfit  , 0.5 + s_cutoff_overfit))
  pred_calib = predict(airf_rqfit, newdata = dplyr::select(airf_test, -y) ,
                       what = c( 0.5 - s_cutoff, 0.5 + s_cutoff ))
  pred_calib_dw = predict(airf_rqfit, newdata = dplyr::select(airf_test, -y) ,
                          what = c( 0.5 - s_cutoff_dw, 0.5 + s_cutoff_dw ))
  pred_calib_dkw_simp = predict(airf_rqfit, newdata = dplyr::select(airf_test, -y) ,
                                what = c( 0.5 - s_cutoff_dkw_simp, 0.5 + s_cutoff_dkw_simp ))
  pred_calib_simes = predict(airf_rqfit, newdata = dplyr::select(airf_test, -y) ,
                          what = c( 0.5 - s_cutoff_simes, 0.5 + s_cutoff_simes ))
  pred_calib_eic = predict(airf_rqfit, newdata = dplyr::select(airf_test, -y) ,
                          what = c( 0.5 - s_cutoff_eic, 0.5 + s_cutoff_eic ))
  
  
  
  cov_overfit = airf_test$y < pred_overfit[,1] | airf_test$y > pred_overfit[,2]
  cov_calib = airf_test$y < pred_calib[,1]  | airf_test$y > pred_calib[,2]
  cov_calib_dw = airf_test$y < pred_calib_dw[,1]  | airf_test$y > pred_calib_dw[,2]
  cov_calib_dkw_simp = airf_test$y < pred_calib_dkw_simp[,1]  | airf_test$y > pred_calib_dkw_simp[,2]
  cov_calib_simes = airf_test$y < pred_calib_simes[,1]  | airf_test$y > pred_calib_simes[,2]
  cov_calib_eic = airf_test$y < pred_calib_eic[,1]  | airf_test$y > pred_calib_eic[,2]
  
  alpha_cov_overfit[idx] = 1- mean(cov_overfit)
  alpha_cov_calib[idx] = 1- mean(cov_calib)
  alpha_cov_calib_dw[idx] = 1- mean(cov_calib_dw)
  alpha_cov_calib_dkw_simp[idx] = 1- mean(cov_calib_dkw_simp)
  alpha_cov_calib_simes[idx] = 1- mean(cov_calib_simes)
  alpha_cov_calib_eic[idx] = 1- mean(cov_calib_eic)
  
  width_overfit = pred_overfit[,2] - pred_overfit[,1]
  width_calib = pred_calib[,2] - pred_calib[,1] 
  width_calib_dw = pred_calib_dw[,2] - pred_calib_dw[,1]
  width_calib_dkw_simp = pred_calib_dkw_simp[,2] - pred_calib_dkw_simp[,1]
  width_calib_simes = pred_calib_simes[,2] - pred_calib_simes[,1]
  width_calib_eic = pred_calib_eic[,2] - pred_calib_eic[,1]
  
  alpha_width_overfit[idx] = mean(width_overfit)
  alpha_width_calib[idx] = mean(width_calib)
  alpha_width_calib_dw[idx] = mean(width_calib_dw)
  alpha_width_calib_dkw_simp[idx] = mean(width_calib_dkw_simp)
  alpha_width_calib_simes[idx] = mean(width_calib_simes)
  alpha_width_calib_eic[idx] = mean(width_calib_eic)
}

df_airf = data.frame( alpha = rep(alpha_fin,  6), 
                      Method = rep(c("training data quantile", "split dcp","DKW simul", "DW simul",
                                     "Simes simul", "Eicker simul"), each = length(alpha_fin)),
                      coverage = rep(NA, 6*length(alpha_fin)),
                      width = rep(NA, 6*length(alpha_fin)))

df_airf$coverage[ (1:length(alpha_fin))] = alpha_cov_overfit
df_airf$coverage[ (1:length(alpha_fin)) + (length(alpha_fin))] = alpha_cov_calib  
df_airf$coverage[ (1:length(alpha_fin)) + (2*length(alpha_fin))] = alpha_cov_calib_dkw_simp
df_airf$coverage[ (1:length(alpha_fin)) + (3*length(alpha_fin))] = alpha_cov_calib_dw
df_airf$coverage[ (1:length(alpha_fin)) + (4*length(alpha_fin))] = alpha_cov_calib_simes
df_airf$coverage[ (1:length(alpha_fin)) + (5*length(alpha_fin))] = alpha_cov_calib_eic

df_airf$width[ (1:length(alpha_fin))] = alpha_width_overfit
df_airf$width[ (1:length(alpha_fin)) + (length(alpha_fin))] = alpha_width_calib  
df_airf$width[ (1:length(alpha_fin)) + (2*length(alpha_fin))] = alpha_width_calib_dkw_simp
df_airf$width[ (1:length(alpha_fin)) + (3*length(alpha_fin))] = alpha_width_calib_dw
df_airf$width[ (1:length(alpha_fin)) + (4*length(alpha_fin))] = alpha_width_calib_simes
df_airf$width[ (1:length(alpha_fin)) + (5*length(alpha_fin))] = alpha_width_calib_eic


alpha_width_overfit = df_airf %>% filter( Method == "training data quantile") %>% 
  ungroup() %>% .$width  %>% as.numeric()

df_airf = df_airf %>% group_by(Method) %>% mutate(diff_width = width - alpha_width_overfit)
df_airf = df_airf %>% group_by(Method) %>% mutate(diff_coverage = coverage + alpha - 1)

p1 = df_airf%>%
  ggplot(mapping = aes(x = alpha, y = diff_coverage, color = Method)) +
  geom_line() + 
  geom_hline(yintercept = 0,linetype = 2, color = "grey") + 
  ylab(expression(paste("Coverage - (1 -", alpha, ")", sep=""))) + 
  ggtitle(expression(paste("Airfoil dataset: Coverage probability ",sep = "")),
          subtitle = expression(paste("with ", delta," = 0.1, I = [0,0.5] ", "(difference with (1-",alpha,"))", sep = ""))) +
  theme(legend.position="bottom")

#save_plot("/Users/sidsarkar/Desktop/Conf/bike/airfoil_plots/airf_cov_int_cow.pdf",p1)

p2 = df_airf%>%
  ggplot(mapping = aes(x = alpha, y = width, color = Method)) +
  geom_line() + 
  ylab("Width") + 
  ggtitle(expression(paste("Airfoil dataset: Average width of prediction set ",sep = "")),
          subtitle = expression(paste("with ", delta," = 0.1, I = [0,0.5] ", sep = ""))) +
  theme(legend.position="bottom")

#save_plot("/Users/sidsarkar/Desktop/Conf/bike/airfoil_plots/airf_width_int_cow.pdf",p2)


################### combined



p1_comb = df_airf%>%
  ggplot(mapping = aes(x = alpha, y = diff_coverage, color = Method)) +
  geom_line() + 
  geom_hline(yintercept = 0,linetype = 2, color = "grey") + 
  ylab(expression(paste("Coverage - (1 -", alpha, ")", sep=""))) + 
  ggtitle(expression(paste("Coverage probability",sep = "")),
          subtitle = expression(paste("with ", delta," = 0.1, I = [0,0.5] ", " (difference with (1-",alpha,"))",  sep = ""))) +
  theme(legend.position="bottom") +
  theme(text = element_text(size = 13)) +
  theme(legend.text=element_text(size=13))

p2_comb = df_airf%>%
  ggplot(mapping = aes(x = alpha, y = width, color = Method)) +
  geom_line() + 
  ylab("Width") + 
  ggtitle(expression(paste("Average width of prediction set",sep = "")),
          subtitle = expression(paste("with ", delta," = 0.1, I = [0,0.5] ", sep = ""))) +
  theme(legend.position="bottom") +
  theme(text = element_text(size = 13)) +
  theme(legend.text=element_text(size=13))

p21 = ggarrange(p1_comb , 
                p2_comb + theme(text = element_text(size = 13)) , 
                ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

p21 = annotate_figure(p21, top = text_grob("Airfoil dataset",size = 17))
p21
ggsave(p21 , file = "/Users/sidsarkar/Desktop/Conf/bike/airfoil_plots/airf_comb_bates.pdf",  width = 11, height = 3.7, dpi = 300, units = "in" )
#save_plot("/Users/sidsarkar/Desktop/Conf/bike/airfoil_plots/airf_comb_bates_cow.pdf",p21)

#####

p1_airf = p1_comb
p2_airf = p2_comb

p21_real_airf = ggarrange(p1_airf, p2_airf, p1_real, p1_real, 
                          ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
p21_real_airf = annotate_figure(p21_real_airf, top = text_grob("Airfoil(top) and F-MNIST(bottom) dataset",size = 17))
ggsave(p21_real_airf , file = "/Users/sidsarkar/Desktop/Conf/bike/airfoil_plots/p21_airf_real_comb_bates.pdf",  width = 11, height = 6, dpi = 300, units = "in" )
