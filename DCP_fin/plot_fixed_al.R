source("/Users/sidsarkar/Desktop/Conf/dkw/get_bds_dkw.R")
source("~/Desktop/Conf/dumbgun_wellner_new/dw_quantile_gen1.R")
source("~/Desktop/Conf/dumbgun_wellner_new/get_bds_dw1.R")

library(dplyr)

low_idx_dw = numeric(19)
low_idx_dkw_simp= numeric(19)
prob_dkw = numeric(19)
prob_dw = numeric(19)

n_calib = 10000

bds_dw = get_bds_dw(n = n_calib, alpha = 0.1) %>% mutate("Empirical.CDF.values" = 1:n/n) %>% mutate(Method = "DW")
bds_dw = bds_dw %>% rename("lower.bound" = "low",
                               "upper.bound" = "up") 

bds_dkw_simp = get_bds_dkw_simp(n = n_calib, alpha = 0.1) %>% mutate("Empirical.CDF.values" = 1:n/n) %>% mutate(Method = "DKW")
bds_dkw_simp = bds_dkw_simp  %>% rename("lower.bound" = "low",
                                        "upper.bound" = "up") 

for(idx in 1:19)
{  
  #print(idx)
  alpha_n = idx/20
  
  ### Duembgun Wellner
  ### which point will give us the lower bound of F(x) >= 1- alpha
  low_idx_dw[idx] = min( which(bds_dw$low >= 1- alpha_n))

  ### DKW simple
  ### which point will give us the lower bound of F(x) >= 1- alpha
  low_idx_dkw_simp[idx] = min( which(bds_dkw_simp$low >= 1- alpha_n))

  prob_dw[idx]  =  pbeta(1 - alpha_n, shape1 = low_idx_dw[idx], 
                         shape2 = n_calib - low_idx_dw[idx] + 1,
                         lower.tail = F)
  prob_dkw[idx]  = pbeta(1 - alpha_n, shape1 = low_idx_dkw_simp[idx], 
                         shape2 = n_calib - low_idx_dkw_simp[idx] + 1,
                         lower.tail = F)
}  

df_plot = data.frame(Method = rep(c("DKW simul","DW simul"), each = 38),
                     alpha = rep(alpha_n , 4),
                     delta= rep(rep(c("True","Lower bound"),each = 19), 2),
                     value = c(prob_dkw, theory_bd_dkw_fs1,
                                prob_dw, theory_bd_dw_fs1))
