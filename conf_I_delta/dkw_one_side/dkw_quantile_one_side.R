#### implementation of dkw for a particular interval I

library(dplyr)


dkw_quant_one_side = function(  delta, I, n, n_runs = 10000, type = "Interval" )
{
  sam = numeric(n_runs)
  
  if(type == "Interval"){
    
    #### if supplied I is an interval, and if it is correctly ordered
    if(I[1]>I[2]){stop("Incorrectly ordered interval I")}
    
    ### set where 1- alpha lies
    I_comp = 1 - I[2:1]
    
    
    for( i in 1:n_runs)
    {
      U_sort = sort( runif(n))
      F_n_hat = ecdf(U_sort)
      
      idx_low_I = (F_n_hat(I_comp[1])*n) + 1
      idx_up_I = (F_n_hat(I_comp[2])*n)
      
      # if(U_sort[idx_up_I] == I_comp[2])
      # {
      #   U_total = c(I_comp[1], U_sort[idx_low_I:idx_up_I])
      # }else{
      #   U_total = c(I_comp[1], U_sort[idx_low_I:idx_up_I], I_comp[2])
      # }
      # 
      # n_total = length(U_total)
      # U_total_tail = U_total[2:length(U_total)]
      # U_total_head = U_total[1:(length(U_total)-1)]
      
      U_total_head =  c(I_comp[1], U_sort[idx_low_I:idx_up_I])
      U_total_tail =  c(U_sort[idx_low_I:idx_up_I], I_comp[2])
      
      
      
      sam[i] = max( sqrt(n)*pmax(F_n_hat(U_total_head)- U_total_head, 0),
                    sqrt(n)*pmax(F_n_hat(U_total_head)- U_total_tail,0))
      }
    
  }else if(type  == "Array"){
    ##### if supplied I is an array
    
    ### set where 1- alpha lies
    I_comp = 1 - I
    
    for( i in 1:n_runs)
    {
      U_sort = sort(runif(n))
      F_n_hat = ecdf(U_sort)
      
      sam[i] = max(sqrt(n)*pmax(F_n_hat(I_comp) - I_comp,0))
    }
  }
  else{
    #### if neither type it is an error
    stop("invalid choice of I")
  }
  
  n_cuttoff = ceiling(n_runs*(1-delta))
  
  return(sort(sam)[n_cuttoff])
}
