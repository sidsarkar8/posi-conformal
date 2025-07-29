###### CLOUD ######
######## plots for dcp
######## One-sided intervalss
### functions now return a correct value of alpha to take directly quantile with

set.seed(2024)

source("/home/siddhaas/Conf_redone/conf_I_delta/dumbgun_wellner_one_side/get_bds_dw_one_side.R")
source("/home/siddhaas/Conf_redone/conf_I_delta/dkw_one_side/get_bds_dkw_one_side.R")

library(dplyr)

n_calib = 10000
delta_n = 0.1

alpha_preset = c(0.5,1)

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

write.csv(dw_quant_val1, file = "/home/siddhaas/Conf_redone/conf_I_delta/dumbgun_wellner_one_side/dw_quant_10000_I.csv")
write.csv(dkw_quant_val1, file = "/home/siddhaas/Conf_redone/conf_I_delta/dkw_one_side/dkw_quant_10000_I.csv")
