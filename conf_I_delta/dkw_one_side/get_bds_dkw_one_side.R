##### bounds generator for DKW simple formulation I
source("~/Documents/Projects/Conf_redone/conf_I_delta/dkw_one_side/dkw_quantile_one_side.R")


alpha_dkw_one_side = function(alpha, I, n, dkw_quant_val = NA, type = "Interval",  v = 3/2)
{
  
  if(type == "Interval"){
    #### if supplied I is an interval, and if it is correctly ordered
    if(I[1]>I[2]){stop("Incorrectly order interval I")}
    
    #### TODO find quantile    
    
    #### using quantile return alpha correction function  
    ########################################################
    if(alpha >= I[1] & alpha<= I[2])
    {
      alpha_comp_corr = 1 - alpha + dkw_quant_val/sqrt(n)
    }else{
      stop("alpha not in range of I interval")
    }
    
    alpha_comp_corr = max(0, alpha_comp_corr)
    alpha_comp_corr = min(alpha_comp_corr, 1)
    
    return(1 - alpha_comp_corr)
    ########################################################
    
  }else if(type == "Array"){
    ##### if supplied I is an array
    
    #### TODO find quantile
    
    #### using quantile return alpha correction 
    ########################################################
    if(alpha %in% I)
    {
      alpha_comp_corr = 1 - alpha + dkw_quant_val/sqrt(n)
    }else{
      stop("alpha not in range of I interval")
    }
    
    alpha_comp_corr = max(0, alpha_comp_corr)
    alpha_comp_corr = min(alpha_comp_corr, 1)
    
    return(1 - alpha_comp_corr)
    ########################################################
    
  }else{
    #### if neither type it is an error
    stop("invalid choice of I")
  }
}
# 
# qq1 = dkw_quant_one_side(0.2, c(0.2,0.3), type = "Interval", 1000, )
# 
# delta1 = 0.1
# n1 = 1000
# alpha1 = 0.2
# qq2 = dkw_quant_one_side(delta = delta1,
#                             I =  alpha1, type = "Array", n = n1, n_runs = 10000)
# 
# alpha_dkw_one_side(alpha = alpha1,  n = n1, I = alpha1, dkw_quant_val = qq2, type = "Array")
# 
# k = 1
# while(1 - pbeta(1- alpha1, shape1 = k, shape2 = n1-k+1) < 1 - delta1){
#   k = k + 1
# }
# 1 - k/n1
# 
# (qbinom(delta1, size = n1,
#         prob = alpha1) + 1)/(n1+1)
