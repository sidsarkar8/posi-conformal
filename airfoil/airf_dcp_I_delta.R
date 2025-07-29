set.seed(2024)
library(dplyr)
library(ranger)
library(quantregForest)
library(tidyverse)
library(cowplot)
library(ggpubr)
##########################################
##########################################


setwd("/Users/sidsarkar/Documents/Projects/Conf_redone/bike")

#####################

airf = read.table("/Users/sidsarkar/Documents/Projects/Conf_redone/bike/airfoil_self_noise.dat")
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

source("/Users/sidsarkar/Documents/Projects/Conf_redone/conf_I_delta/dumbgun_wellner_one_side/get_bds_dw_one_side.R")
source("/Users/sidsarkar/Documents/Projects/Conf_redone/conf_I_delta/dkw_one_side/get_bds_dkw_one_side.R")
 
dw_quant_val1 = read.csv("/Users/sidsarkar/Documents/Projects/Conf_redone/conf_I_delta/dumbgun_wellner_one_side/dw_quant_300_I.csv")[,2]
dkw_quant_val1 = read.csv("/Users/sidsarkar/Documents/Projects/Conf_redone/conf_I_delta/dkw_one_side/dkw_quant_300_I.csv")[,2]


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


#### cov

alpha_cov_overfit = numeric(length(alpha_fin))
alpha_cov_calib = numeric(length(alpha_fin))
alpha_cov_calib_dw = numeric(length(alpha_fin))
alpha_cov_calib_dkw_simp = numeric(length(alpha_fin))

low_cov_overfit = numeric(length(alpha_fin))
low_cov_calib = numeric(length(alpha_fin))
low_cov_calib_dw = numeric(length(alpha_fin))
low_cov_calib_dkw_simp = numeric(length(alpha_fin))

up_cov_overfit = numeric(length(alpha_fin))
up_cov_calib = numeric(length(alpha_fin))
up_cov_calib_dw = numeric(length(alpha_fin))
up_cov_calib_dkw_simp = numeric(length(alpha_fin))

#### width 

alpha_width_overfit = numeric(length(alpha_fin))
alpha_width_calib = numeric(length(alpha_fin))
alpha_width_calib_dw = numeric(length(alpha_fin))
alpha_width_calib_dkw_simp = numeric(length(alpha_fin))

low_width_overfit = numeric(length(alpha_fin))
low_width_calib = numeric(length(alpha_fin))
low_width_calib_dw = numeric(length(alpha_fin))
low_width_calib_dkw_simp = numeric(length(alpha_fin))

up_width_overfit = numeric(length(alpha_fin))
up_width_calib = numeric(length(alpha_fin))
up_width_calib_dw = numeric(length(alpha_fin))
up_width_calib_dkw_simp = numeric(length(alpha_fin))



for(idx in 1:length(alpha_fin))
{  
  print(idx)
  alpha_n = alpha_fin[idx]
  
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
  
  s_cutoff_dkw_simp = quantile(score_n, probs = 1 - alpha_corr_dkw)
  
  ### standard CQR
  
  #s_cutoff = quantile(score_n, probs = 1- alpha_n)
  alpha_corr_split = 1 - (qbinom(1-delta_n, 
                                 size = n_calib, 
                                 prob = 1- alpha_n)/n_calib)
  
  s_cutoff = quantile(score_n, probs = 1- alpha_corr_split)
  
  ### cutoff based on training data
  s_cutoff_overfit = quantile(score_n_overfit, probs = 1 - alpha_n)
  
  
  pred_overfit = predict(airf_rqfit, newdata = dplyr::select(airf_test, -y) ,
                         what = c( 0.5 - s_cutoff_overfit  , 0.5 + s_cutoff_overfit))
  pred_calib = predict(airf_rqfit, newdata = dplyr::select(airf_test, -y) ,
                       what = c( 0.5 - s_cutoff, 0.5 + s_cutoff ))
  pred_calib_dw = predict(airf_rqfit, newdata = dplyr::select(airf_test, -y) ,
                          what = c( 0.5 - s_cutoff_dw, 0.5 + s_cutoff_dw ))
  pred_calib_dkw_simp = predict(airf_rqfit, newdata = dplyr::select(airf_test, -y) ,
                                what = c( 0.5 - s_cutoff_dkw_simp, 0.5 + s_cutoff_dkw_simp ))
  
  cov_overfit = airf_test$y < pred_overfit[,1] | airf_test$y > pred_overfit[,2]
  cov_calib = airf_test$y < pred_calib[,1]  | airf_test$y > pred_calib[,2]
  cov_calib_dw = airf_test$y < pred_calib_dw[,1]  | airf_test$y > pred_calib_dw[,2]
  cov_calib_dkw_simp = airf_test$y < pred_calib_dkw_simp[,1]  | airf_test$y > pred_calib_dkw_simp[,2]
  
  
  ### coverage
  alpha_cov_overfit[idx] = 1- mean(cov_overfit)
  alpha_cov_calib[idx] = 1- mean(cov_calib)
  alpha_cov_calib_dw[idx] = 1- mean(cov_calib_dw)
  alpha_cov_calib_dkw_simp[idx] = 1- mean(cov_calib_dkw_simp)
  
  up_cov_overfit[idx] = quantile( 1-(cov_overfit), probs = 0.95)
  up_cov_calib[idx] =  quantile( 1-(cov_calib), probs = 0.95)
  up_cov_calib_dw[idx] =  quantile( 1-(cov_calib_dw), probs = 0.95)
  up_cov_calib_dkw_simp[idx] =  quantile( 1-(cov_calib_dkw_simp), probs = 0.95) 
  
  low_cov_overfit[idx] =  quantile( 1-(cov_overfit), probs = 0.05)
  low_cov_calib[idx] = quantile( 1-(cov_calib), probs = 0.05)
  low_cov_calib_dw[idx] = quantile( 1-(cov_calib_dw), probs = 0.05)
  low_cov_calib_dkw_simp[idx] =  quantile( 1-(cov_calib_dkw_simp), probs = 0.05) 
  
  #### width 
  
  width_overfit = pred_overfit[,2] - pred_overfit[,1]
  width_calib = pred_calib[,2] - pred_calib[,1] 
  width_calib_dw = pred_calib_dw[,2] - pred_calib_dw[,1]
  width_calib_dkw_simp = pred_calib_dkw_simp[,2] - pred_calib_dkw_simp[,1]
  
  alpha_width_overfit[idx] = mean(width_overfit)
  alpha_width_calib[idx] = mean(width_calib)
  alpha_width_calib_dw[idx] = mean(width_calib_dw)
  alpha_width_calib_dkw_simp[idx] = mean(width_calib_dkw_simp)
  
  up_width_overfit[idx] = quantile(width_overfit, probs = 0.95)
  up_width_calib[idx] = quantile(width_calib, probs = 0.95)
  up_width_calib_dw[idx] = quantile(width_calib_dw, probs = 0.95)
  up_width_calib_dkw_simp[idx] = quantile(width_calib_dkw_simp, probs = 0.95)
  
  low_width_overfit[idx] = quantile(width_overfit, probs = 0.05)
  low_width_calib[idx] = quantile(width_calib, probs = 0.05)
  low_width_calib_dw[idx] = quantile(width_calib_dw, probs = 0.05)
  low_width_calib_dkw_simp[idx] = quantile(width_calib_dkw_simp, probs = 0.05)
}

df_airf = data.frame( alpha = rep(alpha_fin,  4), 
                      Method = rep(c("training data quantile", "split dcp","DKW simul", "DW simul"), each = length(alpha_fin)),
                      coverage = rep(NA, 4*length(alpha_fin)),
                      low_coverage = rep(NA, 4*length(alpha_fin)),
                      up_coverage = rep(NA, 4*length(alpha_fin)),
                      width = rep(NA, 4*length(alpha_fin)),
                      low_width = rep(NA, 4*length(alpha_fin)),
                      up_width = rep(NA, 4*length(alpha_fin)))

##### cov

df_airf$coverage[ (1:length(alpha_fin))] = alpha_cov_overfit
df_airf$coverage[ (1:length(alpha_fin)) + (length(alpha_fin))] = alpha_cov_calib  
df_airf$coverage[ (1:length(alpha_fin)) + (2*length(alpha_fin))] = alpha_cov_calib_dkw_simp
df_airf$coverage[ (1:length(alpha_fin)) + (3*length(alpha_fin))] = alpha_cov_calib_dw

df_airf$low_coverage[ (1:length(alpha_fin))] = low_cov_overfit
df_airf$low_coverage[ (1:length(alpha_fin)) + (length(alpha_fin))] = low_cov_calib  
df_airf$low_coverage[ (1:length(alpha_fin)) + (2*length(alpha_fin))] = low_cov_calib_dkw_simp
df_airf$low_coverage[ (1:length(alpha_fin)) + (3*length(alpha_fin))] = low_cov_calib_dw

df_airf$up_coverage[ (1:length(alpha_fin))] = up_cov_overfit
df_airf$up_coverage[ (1:length(alpha_fin)) + (length(alpha_fin))] = up_cov_calib  
df_airf$up_coverage[ (1:length(alpha_fin)) + (2*length(alpha_fin))] = up_cov_calib_dkw_simp
df_airf$up_coverage[ (1:length(alpha_fin)) + (3*length(alpha_fin))] = up_cov_calib_dw

#### width 

df_airf$width[ (1:length(alpha_fin))] = alpha_width_overfit
df_airf$width[ (1:length(alpha_fin)) + (length(alpha_fin))] = alpha_width_calib  
df_airf$width[ (1:length(alpha_fin)) + (2*length(alpha_fin))] = alpha_width_calib_dkw_simp
df_airf$width[ (1:length(alpha_fin)) + (3*length(alpha_fin))] = alpha_width_calib_dw

df_airf$low_width[ (1:length(alpha_fin))] = low_width_overfit
df_airf$low_width[ (1:length(alpha_fin)) + (length(alpha_fin))] = low_width_calib  
df_airf$low_width[ (1:length(alpha_fin)) + (2*length(alpha_fin))] = low_width_calib_dkw_simp
df_airf$low_width[ (1:length(alpha_fin)) + (3*length(alpha_fin))] = low_width_calib_dw

df_airf$up_width[ (1:length(alpha_fin))] = up_width_overfit
df_airf$up_width[ (1:length(alpha_fin)) + (length(alpha_fin))] = up_width_calib  
df_airf$up_width[ (1:length(alpha_fin)) + (2*length(alpha_fin))] = up_width_calib_dkw_simp
df_airf$up_width[ (1:length(alpha_fin)) + (3*length(alpha_fin))] = up_width_calib_dw

#####


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

p1
#save_plot("/Users/sidsarkar/Desktop/Conf/bike/airfoil_plots/airf_cov_int_cow.pdf",p1)

p2 = df_airf%>%
  ggplot(mapping = aes(x = alpha, y = width, color = Method)) +
  geom_line() + 
  ylab("Width") + 
  ggtitle(expression(paste("Airfoil dataset: Average width of prediction set ",sep = "")),
          subtitle = expression(paste("with ", delta," = 0.1, I = [0,0.5] ", sep = ""))) +
  theme(legend.position="bottom")

p2
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
# ggsave(p21 , file = "/Users/sidsarkar/Desktop/Conf/bike/airfoil_plots/airf_comb.pdf",  width = 11, height = 3.7, dpi = 300, units = "in" )
# save_plot("/Users/sidsarkar/Desktop/Conf/bike/airfoil_plots/airf_comb_cow.pdf",p21)



################################################################################################
################### with conf bands ##########################################
################################################################################################

p_conf_1_comb =  df_airf%>%
  ggplot(mapping = aes(x = alpha, y = diff_coverage, color = Method)) +
  geom_line( linewidth = 1) + 
  geom_ribbon(data = df_airf, 
              aes(ymin = diff_coverage - 1.64*sqrt(coverage*(1-coverage)/n_test),
                  ymax = diff_coverage + 1.64*sqrt(coverage*(1-coverage)/n_test),
                  fill = Method), alpha = 0.2, colour = NA)  + 
  geom_hline(yintercept = 0,linetype = 2, color = "grey") + 
  theme_bw() + 
  ylab(expression(paste("Coverage - (1 -", alpha, ")", sep=""))) + 
  ggtitle(expression(paste("Coverage probability",sep = "")),
          subtitle = expression(paste("with ", delta," = 0.1, I = [0,0.5] ", " (difference with (1-",alpha,"))",  sep = ""))) +
  theme(legend.position="bottom") +
  theme(text = element_text(size = 13)) +
  theme(legend.text=element_text(size=13))
  

p_conf_1_comb  
  
p_conf_2_comb =  df_airf%>%
  ggplot(mapping = aes(x = alpha, y = width, color = Method)) +
  geom_line(linewidth = 1) + 
  geom_ribbon(data = df_airf, 
              aes(ymin = low_width,
                  ymax = up_width,
                  fill = Method), alpha = 0.17, colour = NA) + 
  ylab("Width") + 
  theme_bw() + 
  ggtitle(expression(paste("Average width of prediction set",sep = "")),
          subtitle = expression(paste("with ", delta," = 0.2, I = [0,0.5] ", sep = ""))) +
  theme(legend.position="bottom") +
  theme(text = element_text(size = 13)) +
  theme(legend.text=element_text(size=13)) 
p_conf_2_comb


p_conf_21 = ggarrange(p_conf_1_comb ,
                p_conf_2_comb + theme(text = element_text(size = 13)) ,
                ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

p_conf_21 = annotate_figure(p_conf_21, top = text_grob("Airfoil dataset",size = 17))
p_conf_21

ggsave(p_conf_21,
       file = "/Users/sidsarkar/Documents/Projects/Conf_redone/bike/airfoil_plots/airf_comb_I_del.pdf",
       width = 11, height = 3.7, dpi = 300, units = "in" )

save_plot("/Users/sidsarkar/Documents/Projects/Conf_redone/bike/airfoil_plots/airf_comb_I_del_cow.pdf",
          p_conf_21)
