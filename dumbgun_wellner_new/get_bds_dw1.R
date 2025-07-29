#source("~/Desktop/Conf/dumbgun_wellner_new/dw_quantile_gen1.R")
source("~/Documents/Projects/Conf/dumbgun_wellner_new/dw_quantile_gen1.R")

fin_func_dw = function(p,q,quant,v = 3/2,n)
{
  return(KL_dw(p,q) - ((C_v_dw(p,q,v = v)+ quant)/n)  )
}

get_bds_dw = function(alpha = 0.05, n, v = 3/2, quant_val = NA)
{
  if(is.na(quant_val))
  { 
    dw_quant_val = dw_quant(alpha = alpha, n = n, v = v)
  }else{
    dw_quant_val = quant_val
  }
  
  bds_dw = data.frame(low = rep(NA, n),
                      up = rep(NA, n))
  
  t_n = (1:n)/n
  #t_n_shift = (0:(n-1))/n
  
  for(i in 1:(n-1))
  {
    ### upper bound
    if( fin_func_dw(p = t_n[i], q = 1-10^(-10), quant = dw_quant_val, v = v, n = n) < 0 ){
      temp_up = 1
    }else{
      temp_up = uniroot( fin_func_dw, lower = t_n[i] , upper = 1-10^(-10),
                              p = t_n[i], quant = dw_quant_val, v=v, n = n, tol = 10^-10 )$root
    }
    
    bds_dw$up[i]  = temp_up
    bds_dw$low[n-i] = 1 - temp_up
  }
  
  #### upper for i = n
  bds_dw$up[n] = 1
  
  ### upper bound for i = 0
  if( fin_func_dw(p = 0, q = 1-10^(-10), quant = dw_quant_val, v = v, n = n) < 0 ){
    temp_up_0 = 1
  }else{
    temp_up_0 = uniroot( fin_func_dw, lower = 0 , upper = 1-10^(-10),
                       p = 0, quant = dw_quant_val, v=v, n = n, tol = 10^-10 )$root
  }
  
  ##### upper for i = 0
  bds_dw$low[n] = 1 - temp_up_0
  
  return(bds_dw)
}

#bds_dw = get_bds_dw(n=100)

################################### plotting
## plot method
# plot( y = (1:n)/n, x = 1:n, col = "green", type = "l", ylim = c(-0.5,1.5),
#        main = paste("Dumbgun Wellner, n =", n))
# points( y = bds_dw$low , x = 1:n, col = "red", pch = 4)
# points( y = bds_dw$up, x = 1:n, col = "blue", pch = 4)

## plot method vs Full essHist
# plot( y = (1:n)/n, x = 1:n, col = "green", type = "l", ylim = c(-0.5,1.5),
#       main = paste("Dumbgun Wellner vs Full essHist( bf ), n =", n))
# 
# points( y = bds_dw$low , x = 1:n, col = "red", pch = 4)
# points( y = bds_dw$up, x = 1:n, col = "blue", pch = 4)
# points( y = low_bd_1 , x = 1:n, col = "blue")
# points( y = up_bd_1, x = 1:n, col = "red")

## plot method vs LP
# plot( y = (1:n)/n, x = 1:n, col = "green", type = "l", ylim = c(-0.5,1.5),
#       main = paste("Dumbgun Wellner vs Full essHist( LP ), n =", n))
# 
# points( y = bds_dw$low , x = 1:n, col = "red", pch = 4)
# points( y = bds_dw$up, x = 1:n, col = "blue", pch = 4)
# points( y = low_bd_lp , x = 1:n, col = "blue")
# points( y = up_bd_lp, x = 1:n, col = "red")

