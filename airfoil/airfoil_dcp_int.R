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

alpha_width_overfit = numeric(length(alpha_fin))
alpha_width_calib = numeric(length(alpha_fin))
alpha_width_calib_dw = numeric(length(alpha_fin))
alpha_width_calib_dkw_simp = numeric(length(alpha_fin))
alpha_width_calib_dkw_complex = numeric(length(alpha_fin))


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
  
  ######### 
  
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

  alpha_cov_overfit[idx] = 1- mean(cov_overfit)
  alpha_cov_calib[idx] = 1- mean(cov_calib)
  alpha_cov_calib_dw[idx] = 1- mean(cov_calib_dw)
  alpha_cov_calib_dkw_simp[idx] = 1- mean(cov_calib_dkw_simp)

  width_overfit = pred_overfit[,2] - pred_overfit[,1]
  width_calib = pred_calib[,2] - pred_calib[,1] 
  width_calib_dw = pred_calib_dw[,2] - pred_calib_dw[,1]
  width_calib_dkw_simp = pred_calib_dkw_simp[,2] - pred_calib_dkw_simp[,1]

  alpha_width_overfit[idx] = mean(width_overfit)
  alpha_width_calib[idx] = mean(width_calib)
  alpha_width_calib_dw[idx] = mean(width_calib_dw)
  alpha_width_calib_dkw_simp[idx] = mean(width_calib_dkw_simp)
}

df_airf = data.frame( alpha = rep(alpha_fin,  4), 
                     Method = rep(c("training data quantile", "split dcp","DKW simul", "DW simul"), each = length(alpha_fin)),
                     coverage = rep(NA, 4*length(alpha_fin)),
                     width = rep(NA, 4*length(alpha_fin)))

df_airf$coverage[ (1:length(alpha_fin))] = alpha_cov_overfit
df_airf$coverage[ (1:length(alpha_fin)) + (length(alpha_fin))] = alpha_cov_calib  
df_airf$coverage[ (1:length(alpha_fin)) + (2*length(alpha_fin))] = alpha_cov_calib_dkw_simp
df_airf$coverage[ (1:length(alpha_fin)) + (3*length(alpha_fin))] = alpha_cov_calib_dw

df_airf$width[ (1:length(alpha_fin))] = alpha_width_overfit
df_airf$width[ (1:length(alpha_fin)) + (length(alpha_fin))] = alpha_width_calib  
df_airf$width[ (1:length(alpha_fin)) + (2*length(alpha_fin))] = alpha_width_calib_dkw_simp
df_airf$width[ (1:length(alpha_fin)) + (3*length(alpha_fin))] = alpha_width_calib_dw

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

save_plot("/Users/sidsarkar/Desktop/Conf/bike/airfoil_plots/airf_cov_int_cow.pdf",p1)

p2 = df_airf%>%
  ggplot(mapping = aes(x = alpha, y = width, color = Method)) +
  geom_line() + 
  ylab("Width") + 
  ggtitle(expression(paste("Airfoil dataset: Average width of prediction set ",sep = "")),
          subtitle = expression(paste("with ", delta," = 0.1, I = [0,0.5] ", sep = ""))) +
  theme(legend.position="bottom")

save_plot("/Users/sidsarkar/Desktop/Conf/bike/airfoil_plots/airf_width_int_cow.pdf",p2)


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
ggsave(p21 , file = "/Users/sidsarkar/Desktop/Conf/bike/airfoil_plots/airf_comb.pdf",  width = 11, height = 3.7, dpi = 300, units = "in" )
save_plot("/Users/sidsarkar/Desktop/Conf/bike/airfoil_plots/airf_comb_cow.pdf",p21)
