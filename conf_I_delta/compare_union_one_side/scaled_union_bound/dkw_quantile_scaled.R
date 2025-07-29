#### union-bound method
#### implementation of dkw for a particular interval I

library(dplyr)

dkw_quant_union_finite = function(  delta, I, n, n_runs = 10000)
{
  sam = matrix(nrow = length(I), ncol = n_runs)
  
  
  ##### if supplied I is an array
  
  ### set where 1- alpha lies
  I_comp = 1 - I
  
  for( j in 1:n_runs)
  {
    U_sort = sort(runif(n))
    F_n_hat = ecdf(U_sort)
    
    for( i in 1:length(I))
    {
      sam[i,j] = max(sqrt(n)*pmax(F_n_hat(I_comp[i]) - I_comp[i], 0))
    }
  }
  
  ### sort the samples
  
  sam_sort = sam
  for(i in 1:length(I)){
    sam_sort[i,] = sort(sam[i,])
  }
  
  
  n_cuttoff = ceiling(n_runs*(1-delta))
  
  return(as.numeric(sam_sort[,n_cuttoff]))
}


######


dkw_quant_joint_finite = function(  delta, I, n, n_runs = 10000)
{
  
  temp_quant_finite = dkw_quant_union_finite( delta/length(I), I, n, n_runs)
  
  sam = numeric(n_runs)
  
  ##### if supplied I is an array
  
  ### set where 1-alpha lies
  I_comp = 1 - I
  
  for( j in 1:n_runs)
  {
    U_sort = sort(runif(n))
    F_n_hat = ecdf(U_sort)
    
    sam[j] = max(sqrt(n)*pmax(((F_n_hat(I_comp) - I_comp)/temp_quant_finite),
                               0))
  }
  
  n_cuttoff = ceiling(n_runs*(1-delta))
  
  quant_across_methods = data.frame(alpha = I, union = temp_quant_finite,
                                    joint =  temp_quant_finite* (sort(sam)[n_cuttoff]))
  return(quant_across_methods)
}


#dkw_quant_union_finite(0.1, I = c(0.1,0.2,0.3), n = 100)
dkw_quant_joint_finite(0.1, I = c(0.1, 0.11, 0.12,0.13), n = 100, n_runs = 500000)

