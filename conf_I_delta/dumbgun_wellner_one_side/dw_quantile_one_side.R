#### implementation of dumbgun wellner for a particular interval I
library(dplyr)
library(dplyr)

#########################
KL_dw = function(p,q){
  if(length(p)!= length(q))
  {
    print("p and q length not same")
    stop()
  }
  
  val = numeric(length(p))
  
  for( i in 1:length(p))
  {
    val[i] = p[i]*log(p[i]/q[i]) + (1-p[i])*log((1-p[i])/(1-q[i]))
    if( p[i] == 0 ){val[i] = log(1/(1-q[i]))}
    if( p[i] == 1 ){val[i] = log(1/q[i]) }
    if( q[i] == 0 ){val[i] = Inf }
    if( q[i] == 1 ){val[i] = Inf}
    
    if( p[i] == 0 & q[i] == 0 ){val[i] = 0}
    if( p[i] == 1 & q[i] == 0 ){val[i] = Inf }
    if( p[i] == 0 & q[i] == 1 ){val[i] = Inf }
    if( p[i] == 1 & q[i] == 1 ){val[i] = 0}
  }
  
  return(val)
}

#########################
C_dw = function(t)
{
  return( log( 1 - log(2*t) - log(2*(1-t)) ) )
}

#########################
D_dw = function(t)
{
  return( log(1 + (C_dw(t))^2) )
}

C_v_dw = function(u,t, v)
{
  len_t = length(t)
  if(length(u)!= len_t)
  {
    print("not same length")
    stop
  }
  
  temp = numeric(len_t)
  for( i in 1:len_t)
  {
    if(max(u[i],t[i]) < 1/2){
      temp[i] = C_dw(max(u[i],t[i])) + (v*D_dw(max(u[i],t[i])))
    } else if(min(u[i],t[i]) > 1/2){
      temp[i] = C_dw(min(u[i],t[i])) + (v*D_dw(min(u[i],t[i])))
    } else {
      temp[i] = 0
    }
  }
  return(temp)
}

#########################

dw_quant_one_side = function(  delta, I, n, n_runs = 10000, v = 3/2,
                               type = "Interval" )
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
      
      ### case where there is a data-point at top intervaln
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
      
      ### 2 cases to be computed for [X_(i),X_(i+1)) 
      ### one is d(i,i) and other is d(i,i+1)
      
      U_total_head =  c(I_comp[1], U_sort[idx_low_I:idx_up_I])
      U_total_tail =  c(U_sort[idx_low_I:idx_up_I], I_comp[2])
      
      #### only compute when difference is positive
      idx_pos_1 = (F_n_hat(U_total_head) - U_total_head) > 0 
      idx_pos_2 = (F_n_hat(U_total_head) - U_total_tail) > 0
      
      #### if no terms match the condition above, then return 0
      if(sum(idx_pos_1) > 0)
      {
        dw_stat_full_1 = n*KL_dw( F_n_hat(U_total_head), U_total_head) - 
          C_v_dw( F_n_hat(U_total_head), U_total_head, v = v)
        dw_stat_pos_1 = dw_stat_full_1[idx_pos_1]
      }else{
        dw_stat_pos_1 = -Inf
      }
      
      if(sum(idx_pos_2) > 0)
      {
        dw_stat_full_2 = n*KL_dw( F_n_hat(U_total_head), U_total_tail) - 
          C_v_dw( F_n_hat(U_total_head), U_total_tail, v = v)                  
        dw_stat_pos_2 = dw_stat_full_2[idx_pos_2]
      }else{
        dw_stat_pos_2 = -Inf
      }
      
      
      sam[i] = max(dw_stat_pos_1, 
                   dw_stat_pos_2)
    }
    
  }else if(type  == "Array"){
    ##### if supplied I is an array
    
    ### set where 1- alpha lies
    I_comp = 1 - I
    
    for( i in 1:n_runs)
    {
      U_sort = sort(runif(n))
      F_n_hat = ecdf(U_sort)
      
      idx_pos = F_n_hat(I_comp) - I_comp > 0
      
      if(sum(idx_pos) > 0)
      {      
        dw_stat_full = n*KL_dw( F_n_hat(I_comp), I_comp) - 
        C_v_dw( F_n_hat(I_comp), I_comp, v = v)
        dw_stat_pos = dw_stat_full[idx_pos]
      }else{
        dw_stat_pos = -Inf
      }      
      sam[i] = max(dw_stat_pos)
    }
  }
  else{
    #### if neither type it is an error
    stop("invalid choice of I")
  }
  
  n_cuttoff = ceiling(n_runs*(1-delta))
  
  return(sort(sam)[n_cuttoff])
}






