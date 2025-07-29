########################################################

source("~/Documents/Projects/Conf_redone/conf_I_delta/dumbgun_wellner_one_side/dw_quantile_one_side.R")

fin_func_dw = function(p,q,quant,v = 3/2,n)
{
  return(KL_dw(p,q) - ((C_v_dw(p,q,v = v)+ quant)/n)  )
}

alpha_dw_one_side  = function(alpha, I, n, dw_quant_val = NA, type = "Interval",
                        v = 3/2)
{
  
  if(type == "Interval"){
    #### if supplied I is an interval, and if it is correctly ordered
    if(I[1]>I[2]){stop("Incorrectly order interval I")}
    
    #### TODO find quantile    
    
    #### using quantile return alpha correction function  
    ########################################################
    if(alpha >= I[1] & alpha<= I[2])
    {
      #### accounting for quantile, leading to a correction factor
      if( fin_func_dw(q = 1-alpha, p = 1-10^(-10), quant = dw_quant_val, v = v, n = n) < 0 ){
        alpha_comp_corr = 1
      }else{
        alpha_comp_corr = uniroot( fin_func_dw, lower = 1-alpha, upper = 1-10^(-10),
                                   q = 1-alpha, quant = dw_quant_val, v=v, n = n, tol = 10^-10 )$root
      }
    }else{
      stop("alpha not in range of I interval")
    }
    
    alpha_comp_corr = max(0, alpha_comp_corr)
    alpha_comp_corr = min(alpha_comp_corr, 1)
    
    return(1-alpha_comp_corr)
    ########################################################
    
  }else if(type == "Array"){
    ##### if supplied I is an array
    
    #### TODO find quantile
    
    #### using quantile return alpha correction 
    ########################################################
    if(alpha %in% I)
    {
      #### accounting for quantile, leading to a correction factor
      if( fin_func_dw(q = 1-alpha, p = 1-10^(-10), quant = dw_quant_val, v = v, n = n) < 0 ){
        alpha_comp_corr = 1
      }else{
        alpha_comp_corr = uniroot( fin_func_dw, lower = 1-alpha, upper = 1-10^(-10),
                                   q = 1-alpha, quant = dw_quant_val, v=v, n = n, tol = 10^-10 )$root
      }
      
    }else{
      stop("alpha not in range of I interval")
    }
    
    alpha_comp_corr = max(0, alpha_comp_corr)
    alpha_comp_corr = min(alpha_comp_corr, 1)
    
    return(1-alpha_comp_corr)
    ########################################################
    
  }else{
    #### if neither type it is an error
    stop("invalid choice of I")
  }
}


###########

# delta1 = 0.1
# n1 = 500
# alpha1 = 0.2
# qq2 = dw_quant_one_side(delta = delta1,
#                             I =  alpha1, type = "Array", n = n1, n_runs = 50000)
# 
# n1*alpha_dw_one_side(alpha = alpha1,  n = n1, I = alpha1, dw_quant_val = qq2, type = "Array")
# 
# k = 1
# while(1 - pbeta(1- alpha1, shape1 = k, shape2 = n1-k+1) < 1 - delta1){
#   k = k + 1
# }
# 
# n1*(1 - (k-1)/n1)
# 
# (n1 - qbinom(1-delta1, size = n1, prob = 1- alpha1))


# ans = numeric(length(1:100))
# 
# for( idx in 1:100)
# {
#   n1 = idx*10
#   k = 1
#   while(1 - pbeta(1- alpha1, shape1 = k, shape2 = n1-k+1) < 1 - delta1){
#     k = k + 1
#   }
# 
#   # print(idx)
#   # print(n1*(1 - (k + 1)/n1))
#   # print(n1*(1 - (k)/n1))
#   
#   ans[idx] =  n1*(1 - (k-1)/n1) - (n1 - qbinom(1-delta1, size = n1, prob = 1- alpha1 ))
# }
# ans
# # n1*(qbinom(delta1, size = n1, 
# #         prob = alpha1) + 1)/(n1+1)

