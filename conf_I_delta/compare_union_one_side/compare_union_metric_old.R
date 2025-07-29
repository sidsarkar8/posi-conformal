#### metric specific union bound

setwd("~/Documents/Projects/Conf_redone/conf_I_delta/compare_union_one_side")

source("~/Documents/Projects/Conf_redone/conf_I_delta/dumbgun_wellner_one_side/get_bds_dw_one_side.R")
source("~/Documents/Projects/Conf_redone/conf_I_delta/dkw_one_side/get_bds_dkw_one_side.R")

set.seed(2024)

library(dplyr)
library(quantregForest)
library(ggplot2)
library(ggpubr)

##########################################
##########################################

n_calib = 1000
delta_n = 0.1


num_k_seq = c(2,4,8,16)

alpha_close_tail = list( seq(from = 0.1, to = 0.2, length.out = num_k_seq[1]),
                         seq(from = 0.1, to = 0.2, length.out = num_k_seq[2]),
                         seq(from = 0.1, to = 0.2, length.out = num_k_seq[3]),
                         seq(from = 0.1, to = 0.2, length.out = num_k_seq[4]))


alpha_close_mid = list( seq(from = 0.45, to = 0.55, length.out = num_k_seq[1]),
                        seq(from = 0.45, to = 0.55, length.out = num_k_seq[2]),
                        seq(from = 0.45, to = 0.55, length.out = num_k_seq[3]),
                        seq(from = 0.45, to = 0.55, length.out = num_k_seq[4]))

#### father apart alpha's
alpha_far = list( seq(from = 0.02, to = 0.8, length.out = num_k_seq[1]),
                  seq(from = 0.02, to = 0.8, length.out =  num_k_seq[2]),
                  seq(from = 0.02, to = 0.8, length.out = num_k_seq[3]),
                  seq(from = 0.02, to = 0.8, length.out = num_k_seq[4]))


########################################################
############## close at tail #################
########################################################

plot1 = list()

for(num_k_idx in 1:4)
{
  print(num_k_idx)
  
  alpha_preset = alpha_close_tail[[num_k_idx]]
  
  ############ new I,delta methods #####
  
  dw_quant_val1 = dw_quant_one_side(delta = delta_n,
                                    I = alpha_preset,
                                    n = n_calib,
                                    n_runs = 50000,
                                    type = "Array" )
  
  dkw_quant_val1 = dkw_quant_one_side(delta = delta_n,
                                      I = alpha_preset,
                                      n = n_calib,
                                      n_runs = 50000,
                                      type = "Array" )
  ############ 
  temp_df = expand.grid(Method = c("DKW","DW","DKW-union", "DW-union"), alpha = alpha_preset)
  temp_df = temp_df %>% mutate(alpha.corr = NA)
  
  for( i in 1:length(alpha_preset))
  {
    alpha_n = alpha_preset[i]
    
    dw_quant_val_temp = dw_quant_one_side(delta = delta_n/length(alpha_preset),
                                          I = alpha_n,
                                          n = n_calib,
                                          n_runs = 50000,
                                          type = "Array" )
    
    dkw_quant_val_temp = dkw_quant_one_side(delta = delta_n/length(alpha_preset),
                                            I = alpha_n,
                                            n = n_calib,
                                            n_runs = 50000,
                                            type = "Array" )
    
    temp_df$alpha.corr[(i-1)*4 + 1  ] = alpha_dkw_one_side(alpha = alpha_n,
                                                           I = alpha_preset,
                                                           n = n_calib,
                                                           dkw_quant_val = dkw_quant_val1,
                                                           type = "Array")
    
    temp_df$alpha.corr[(i-1)*4 + 2  ] = alpha_dw_one_side(alpha = alpha_n,
                                                          I = alpha_preset,
                                                          n = n_calib,
                                                          dw_quant_val = dw_quant_val1,
                                                          type = "Array")
    
    # temp_df$alpha.corr[(i-1)*4 + 3 ] = 1 - (qbinom(1-(delta_n/length(alpha_preset)), 
    #                                                size = n_calib, prob = 1- alpha_n)/n_calib)
    temp_df$alpha.corr[(i-1)*4 + 3 ] =  alpha_dkw_one_side(alpha = alpha_n,
                                                           I = alpha_n,
                                                           n = n_calib,
                                                           dkw_quant_val = dkw_quant_val_temp,
                                                           type = "Array")
    
    temp_df$alpha.corr[(i-1)*4 + 4 ] =  alpha_dw_one_side(alpha = alpha_n,
                                                          I = alpha_n,
                                                          n = n_calib,
                                                          dw_quant_val = dw_quant_val_temp,
                                                          type = "Array")
    
  }
  
  p_temp = temp_df %>%
    ggplot(mapping = aes(x = alpha, y = alpha.corr - alpha , color = Method)) +
    geom_line() + geom_point() + theme_bw() + 
    geom_hline(yintercept = 0,linetype = 2) +
    theme(text = element_text(size = 10)) + 
    scale_x_continuous(breaks = round(alpha_preset,3)) + 
    labs(x = NULL, y = NULL) + 
    ylim(-0.05,0) + 
    theme(legend.position = "none")
  
  plot1[[num_k_idx]] = p_temp
  
  #    ggtitle(expression(paste("Comparing union bd vs. simul, tail [0.1,0.2], n = 5000, ",
  #                             delta," = 0.1" ,sep = "")),
  #            subtitle = expression(paste("(difference with ",alpha,")", sep = ""))) +
  # ggsave(paste("plots_compare_union_alphas_5000/tail",num_k_idx,"_5000_metric.pdf", sep = ""),
  #        width = 8, height = 6, dpi = 300, units = "in" )
}

plot1_comb = ggpubr::ggarrange(plot1[[1]], plot1[[2]],
                               plot1[[3]], plot1[[4]],
                               common.legend = TRUE, legend="bottom",
                               ncol = 1)

plot1_comb = annotate_figure(plot1_comb,
                             top = text_grob(expression(paste("Comparing union bd vs. simul, [0.1,0.2], n = 1000, ",
                                                              delta," = 0.1" ,sep = "")), 
                                             size = 12),
                             left = text_grob("alpha - alpha.corr", 
                                              rot = 90))
plot1_comb
# ggsave("plot_tail_5000_metric.pdf",
#        width = 6, height = 6, dpi = 300, units = "in" )


########################################################
############## close at middle #################
########################################################


plot2 = list()

for(num_k_idx in 1:4)
{
  print(num_k_idx)
  
  alpha_preset = alpha_close_mid[[num_k_idx]]
  ############ new I,delta methods #####
  
  dw_quant_val1 = dw_quant_one_side(delta = delta_n,
                                    I = alpha_preset,
                                    n = n_calib,
                                    n_runs = 50000,
                                    type = "Array" )
  
  dkw_quant_val1 = dkw_quant_one_side(delta = delta_n,
                                      I = alpha_preset,
                                      n = n_calib,
                                      n_runs = 50000,
                                      type = "Array" )
  ############ 
  temp_df = expand.grid(Method = c("DKW","DW","DKW-union", "DW-union"), alpha = alpha_preset)
  temp_df = temp_df %>% mutate(alpha.corr = NA)
  
  for( i in 1:length(alpha_preset))
  {
    alpha_n = alpha_preset[i]
    
    dw_quant_val_temp = dw_quant_one_side(delta = delta_n/length(alpha_preset),
                                          I = alpha_n,
                                          n = n_calib,
                                          n_runs = 50000,
                                          type = "Array" )
    
    dkw_quant_val_temp = dkw_quant_one_side(delta = delta_n/length(alpha_preset),
                                            I = alpha_n,
                                            n = n_calib,
                                            n_runs = 50000,
                                            type = "Array" )
    
    temp_df$alpha.corr[(i-1)*4 + 1  ] = alpha_dkw_one_side(alpha = alpha_n,
                                                           I = alpha_preset,
                                                           n = n_calib,
                                                           dkw_quant_val = dkw_quant_val1,
                                                           type = "Array")
    
    temp_df$alpha.corr[(i-1)*4 + 2  ] = alpha_dw_one_side(alpha = alpha_n,
                                                          I = alpha_preset,
                                                          n = n_calib,
                                                          dw_quant_val = dw_quant_val1,
                                                          type = "Array")
    
    # temp_df$alpha.corr[(i-1)*4 + 3 ] = 1 - (qbinom(1-(delta_n/length(alpha_preset)), 
    #                                                size = n_calib, prob = 1- alpha_n)/n_calib)
    temp_df$alpha.corr[(i-1)*4 + 3 ] =  alpha_dkw_one_side(alpha = alpha_n,
                                                           I = alpha_n,
                                                           n = n_calib,
                                                           dkw_quant_val = dkw_quant_val_temp,
                                                           type = "Array")
    
    temp_df$alpha.corr[(i-1)*4 + 4 ] =  alpha_dw_one_side(alpha = alpha_n,
                                                          I = alpha_n,
                                                          n = n_calib,
                                                          dw_quant_val = dw_quant_val_temp,
                                                          type = "Array")
    
  }  
  
  p_temp = temp_df %>%
    ggplot(mapping = aes(x = alpha, y = alpha.corr - alpha , color = Method)) +
    geom_line() + geom_point() + theme_bw() + 
    geom_hline(yintercept = 0,linetype = 2) +
    theme(text = element_text(size = 10)) + 
    scale_x_continuous(breaks = round(alpha_preset,3)) + 
    labs(x = NULL, y = NULL) + 
    ylim(-0.05,0) + 
    theme(legend.position = "none")
  
  plot2[[num_k_idx]] = p_temp
  
  #    ggtitle(expression(paste("Comparing union bd vs. simul, tail [0.1,0.2], n = 5000, ",
  #                             delta," = 0.1" ,sep = "")),
  #            subtitle = expression(paste("(difference with ",alpha,")", sep = ""))) +
  # ggsave(paste("plots_compare_union_alphas_5000/tail",num_k_idx,"_5000_metric.pdf", sep = ""),
  #        width = 8, height = 6, dpi = 300, units = "in" )
}

plot2_comb = ggpubr::ggarrange(plot2[[1]], plot2[[2]],
                               plot2[[3]], plot2[[4]],
                               common.legend = TRUE, legend="bottom",
                               ncol = 1)

plot2_comb = annotate_figure(plot2_comb,
                             top = text_grob(expression(paste("Comparing union bd vs. simul, [0.45,0.55], n = 5000, ",
                                                              delta," = 0.1" ,sep = "")), 
                                             size = 12),
                             left = text_grob("alpha - alpha.corr", 
                                              rot = 90))
plot2_comb
# ggsave("plot_middle_5000_metric.pdf",
#        width = 6, height = 6, dpi = 300, units = "in" )

########################################################
############## far away values #################
########################################################

plot3 = list()

for(num_k_idx in 1:4)
{
  print(num_k_idx)
  
  alpha_preset = alpha_far[[num_k_idx]]
  ############ new I,delta methods #####
  
  dw_quant_val1 = dw_quant_one_side(delta = delta_n,
                                    I = alpha_preset,
                                    n = n_calib,
                                    n_runs = 50000,
                                    type = "Array" )
  
  dkw_quant_val1 = dkw_quant_one_side(delta = delta_n,
                                      I = alpha_preset,
                                      n = n_calib,
                                      n_runs = 50000,
                                      type = "Array" )
  ############ 
  temp_df = expand.grid(Method = c("DKW","DW","DKW-union", "DW-union"), alpha = alpha_preset)
  temp_df = temp_df %>% mutate(alpha.corr = NA)
  
  for( i in 1:length(alpha_preset))
  {
    alpha_n = alpha_preset[i]
    
    dw_quant_val_temp = dw_quant_one_side(delta = delta_n/length(alpha_preset),
                                          I = alpha_n,
                                          n = n_calib,
                                          n_runs = 50000,
                                          type = "Array" )
    
    dkw_quant_val_temp = dkw_quant_one_side(delta = delta_n/length(alpha_preset),
                                            I = alpha_n,
                                            n = n_calib,
                                            n_runs = 50000,
                                            type = "Array" )
    
    temp_df$alpha.corr[(i-1)*4 + 1  ] = alpha_dkw_one_side(alpha = alpha_n,
                                                           I = alpha_preset,
                                                           n = n_calib,
                                                           dkw_quant_val = dkw_quant_val1,
                                                           type = "Array")
    
    temp_df$alpha.corr[(i-1)*4 + 2  ] = alpha_dw_one_side(alpha = alpha_n,
                                                          I = alpha_preset,
                                                          n = n_calib,
                                                          dw_quant_val = dw_quant_val1,
                                                          type = "Array")
    
    # temp_df$alpha.corr[(i-1)*4 + 3 ] = 1 - (qbinom(1-(delta_n/length(alpha_preset)), 
    #                                                size = n_calib, prob = 1- alpha_n)/n_calib)
    temp_df$alpha.corr[(i-1)*4 + 3 ] =  alpha_dkw_one_side(alpha = alpha_n,
                                                           I = alpha_n,
                                                           n = n_calib,
                                                           dkw_quant_val = dkw_quant_val_temp,
                                                           type = "Array")
    
    temp_df$alpha.corr[(i-1)*4 + 4 ] =  alpha_dw_one_side(alpha = alpha_n,
                                                          I = alpha_n,
                                                          n = n_calib,
                                                          dw_quant_val = dw_quant_val_temp,
                                                          type = "Array")
    
  }
  
  
  p_temp = temp_df %>%
    ggplot(mapping = aes(x = alpha, y = alpha.corr - alpha , color = Method)) +
    geom_line() + geom_point() + theme_bw() + 
    geom_hline(yintercept = 0,linetype = 2) +
    theme(text = element_text(size = 10)) + 
    scale_x_continuous(breaks = round(alpha_preset,3)) + 
    labs(x = NULL, y = NULL) + 
    ylim(-0.05,0) + 
    theme(legend.position = "none")
  
  plot3[[num_k_idx]] = p_temp
  
  #    ggtitle(expression(paste("Comparing union bd vs. simul, tail [0.1,0.2], n = 5000, ",
  #                             delta," = 0.1" ,sep = "")),
  #            subtitle = expression(paste("(difference with ",alpha,")", sep = ""))) +
  # ggsave(paste("plots_compare_union_alphas_5000/tail",num_k_idx,"_5000_metric.pdf", sep = ""),
  #        width = 8, height = 6, dpi = 300, units = "in" )
}

plot3_comb = ggpubr::ggarrange(plot3[[1]], plot3[[2]],
                               plot3[[3]], plot3[[4]],
                               common.legend = TRUE, legend="bottom",
                               ncol = 1)

plot3_comb = annotate_figure(plot3_comb,
                             top = text_grob(expression(paste("Comparing union bd vs. simul, [0.02,0.8], n = 5000, ",
                                                              delta," = 0.1" ,sep = "")), 
                                             size = 12),
                             left = text_grob("alpha - alpha.corr", 
                                              rot = 90))
plot3_comb
# ggsave("plot_far_5000_metric.pdf",
#        width = 6, height = 6, dpi = 300, units = "in" )

##########

# plot_all = ggarrange(
#   plot1[[1]] + ylim(c(-0.05,0)), plot1[[2]] + ylim(c(-0.05,0)), plot1[[3]] + ylim(c(-0.05,0)), 
#   plot2[[1]] + ylim(c(-0.05,0)), plot2[[2]] + ylim(c(-0.05,0)), plot2[[3]] + ylim(c(-0.05,0)), 
#   plot3[[1]] + ylim(c(-0.05,0)), plot3[[2]] + ylim(c(-0.05,0)), plot3[[3]] + ylim(c(-0.05,0))
# )
# 
# 
# 
# plot_all = annotate_figure(plot_all,
#                 top = text_grob(expression(paste("Comparing union bd vs. simul, n = 5000, ",
#                                                  delta," = 0.1" ,sep = "")), 
#                                 size = 12),
#                 left = text_grob("alpha - alpha.corr", 
#                                  rot = 90))
# 
# annotate_figure(plot_all,
#                 top = text_grob(expression(paste("Comparing union bd vs. simul, n = 5000, ",
#                                                  delta," = 0.1" ,sep = "")), 
#                                 size = 12),
#                 left = text_grob("alpha - alpha.corr \n
#                                  [0.1,0.2]                                        [0.45,0.55]                                        [0.02,0.8]", 
#                                  rot = 90))


plot_all = ggarrange(
  plot1[[1]] + ylim(c(-0.04,0)) + ggtitle("[0.1,0.2]") + theme(legend.text = element_text(size = 11),
                                                               legend.title = element_text(size = 11)), 
  plot2[[1]] + ylim(c(-0.04,0)) + ggtitle("[0.45,0.55]")+ theme(legend.text = element_text(size = 11)), 
  plot3[[1]] + ylim(c(-0.04,0)) + ggtitle("[0.02,0.8]")+ theme(legend.text = element_text(size = 11)),
  plot1[[2]] + ylim(c(-0.04,0))+ theme(legend.text = element_text(size = 11)), 
  plot2[[2]] + ylim(c(-0.04,0))+ theme(legend.text = element_text(size = 11)),
  plot3[[2]] + ylim(c(-0.04,0))+ theme(legend.text = element_text(size = 11)),
  plot1[[3]] + ylim(c(-0.04,0))+ theme(legend.text = element_text(size = 11)), 
  plot2[[3]] + ylim(c(-0.04,0))+ theme(legend.text = element_text(size = 11)), 
  plot3[[3]] + ylim(c(-0.04,0))+ theme(legend.text = element_text(size = 11)),
  common.legend = TRUE, legend="bottom")



plot_all = annotate_figure(plot_all,
                           top = text_grob(expression(paste("Comparing union bd vs. simul, n = 1000, ",
                                                            delta," = 0.1" ,sep = "")),
                                           size = 14),
                           left = text_grob("alpha - alpha.corr",
                                            rot = 90))

plot_all

ggsave("plot_all_one_side_1000_metric_old.pdf",
       width = 9, height = 7, dpi = 300, units = "in" )
##### validity

# n_sam = 5000
# cov_al_cor = matrix(nrow = n_sam, ncol = 3)
# colnames(cov_al_cor) = c("Union", "DKW", "DW")
# 
# for( idx in 1:n_sam)
# {
#   if(!(idx%%1000)){print(idx)}
#   temp_data = runif(n_calib)
#   F_n_temp = ecdf(temp_data)
#   
#   q1 = temp_df%>% filter(Method == "Union")
#   cov_al_cor[idx,1] = prod(F_n_temp(1 - q1$alpha) <= 1 - q1$alpha.corr) 
#   
#   q2 = temp_df%>% filter(Method == "dkw_one_side")
#   cov_al_cor[idx,2] = prod(F_n_temp(1 - q2$alpha) <= 1 - q2$alpha.corr) 
#   
#   
#   q3 = temp_df%>% filter(Method == "DW")
#   cov_al_cor[idx,3] = prod(F_n_temp(1 - q3$alpha) <= 1 - q3$alpha.corr) 
# }
# 
# mean(cov_al_cor[,1]) 
# mean(cov_al_cor[,2])
# mean(cov_al_cor[,3])
# delta_n

