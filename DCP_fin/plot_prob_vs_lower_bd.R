library(ggplot2)
library(dplyr)

setwd("~/Desktop/Conf/DCP_fin")

########################################################################

method_names = c("training data quantile", "split conformal","DKW simul", "DW simul")

####################################################################
source("~/Desktop/Conf/dumbgun_wellner_new/get_bds_dw1.R")
source("~/Desktop/Conf/dkw/get_bds_dkw.R")

delta_n = 0.1
n_calib = 10000
j_norm = (1:n_calib)/n_calib

#############lower bound ###########################################

alpha_n = 1:19/20
del_m_dkw = (1/n_calib + sqrt( log(2/delta_n)/(2*n_calib)))/3
bds_dkw = get_bds_dkw_simp(alpha = delta_n, n = n_calib)

theory_bd_dkw_fs1 = numeric(length(alpha_n))
theory_bd_dkw_exact = numeric(length(alpha_n))


for(i in 1:19)
{
  ## lower bd
  if(alpha_n[i]<1/2){
    alpha_temp = alpha_n[i]
    print("<")
  } else if(alpha_n[i]>1/2 + del_m_dkw){
    alpha_temp = alpha_n[i] - del_m_dkw
    print(">")
  }else{
    alpha_temp = 1/2
    print("=")
  }
  
  theory_bd_dkw_fs1[i] = 1 - (0.05)^(1/(4*alpha_temp*(1- alpha_temp))) 
  
  ## exact
  j_al = which(bds_dkw$low >= 1- alpha_n[i] ) %>% min()
  theory_bd_dkw_exact[i] = 1 - pbeta(1-alpha_n[i], shape1 = j_al, 
                                 shape2 = (n_calib - j_al + 1) )
  
}

#################################################

quant_df_val  = dw_quant(alpha = delta_n, n = n_calib)

v = 3/2

bds_dw_1 =  get_bds_dw(alpha = delta_n, n = n_calib, v = v, quant_val = quant_df_val)
del_m_dw = max(bds_dw_1$up-bds_dw_1$low)/2

theory_bd_dw_fs1 = numeric(length(alpha_n))
theory_bd_dw_exact = numeric(length(alpha_n))

for(i in 1:19)
{
  if(alpha_n[i]<1/2){
    alpha_temp = alpha_n[i]
    print("<")
  } else if(alpha_n[i]>1/2 + del_m_dw){
    alpha_temp = alpha_n[i] - del_m_dw
    print(">")
  }else{
    alpha_temp = 1/2
    print("=")
  }
  
  theory_bd_dw_fs1[i] = 1 - exp( - C_dw(1-alpha_temp) - v*D_dw(1-alpha_temp) - quant_df_val )
  
  ## exact
  j_al = which(bds_dw_1$low >= 1- alpha_n[i] ) %>% min()
  theory_bd_dw_exact[i] =1 - pbeta(1-alpha_n[i], shape1 = j_al, 
                                 shape2 = (n_calib - j_al + 1) )
}


####################################################

df_plot = data.frame(Method = rep(c("DKW simul","DW simul"), each = 38),
                     alpha = rep(alpha_n , 4),
                     delta= rep(rep(c("Exact","Lower bound"),each = 19), 2),
                     value = c( theory_bd_dkw_exact, theory_bd_dkw_fs1,
                                theory_bd_dw_exact, theory_bd_dw_fs1))

plt_fix_al = df_plot %>% 
  ggplot(mapping = aes(x = alpha, y = value, color = delta)) +
  geom_line(aes(group = delta)) +
  facet_grid(~ Method) + 
  scale_x_continuous(breaks = seq(1/20, 19/20, by = 0.2)) +
  ylab("probability") +
  ggtitle(expression(paste("P( Coverage ">="1-",alpha," ) for fixed ", alpha, " with ", delta, " = 0.1, n = 10000",  sep = ""))) +
  theme(legend.position="bottom") +
  theme(text = element_text(size = 10)) 
plt_fix_al

plt_fix_al_same = df_plot  %>%
  ggplot(mapping = aes(x = alpha, y = value, color = Method, linetype = delta)) +
  geom_line() + 
  coord_cartesian(ylim=c(0.985,1)) + 
  scale_x_continuous(breaks = seq(1/20, 19/20, by = 0.2)) +
  ylab("probability") +
  ggtitle(expression(paste("P( Coverage ">="1-",alpha," ) for fixed ", alpha,  sep = "")),
          subtitle = expression(paste("with ", delta," = 0.1, calibration set size = 10000 ", sep = ""))) +
  theme(legend.position="bottom") +
  labs(color = "", linetype = "") + 
  theme(text = element_text(size = 12)) 

plt_fix_al_same

ggsave("./fixed_al_exact_vs_lowerbd_same.pdf", width = 5, height = 3.5, dpi = 300, units = "in",
       plot = plt_fix_al_same)


