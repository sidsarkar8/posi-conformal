set.seed(2022)
library(dplyr)
library(ranger)
library(quantregForest)
library(tidyverse)
##########################################
##########################################

##### bounds generator for DKW simple formulation


get_bds_dkw_simp = function(alpha = 0.05, n)
{
  bds_dkw_simp = data.frame(low = rep(NA, n),
                            up = rep(NA, n))
  
  bds_dkw_simp$low = (1:n)/n - sqrt(( log(2) - log(alpha))/(2*n))
  bds_dkw_simp$up = (1:n)/n + sqrt(( log(2) - log(alpha))/(2*n))
  
  bds_dkw_simp$low[ bds_dkw_simp$low < 0] = 0
  bds_dkw_simp$up[ bds_dkw_simp$up > 1] = 1
  
  return(bds_dkw_simp)
}

lambda_uppbd_func = function(lambda, n, alpha )
{
  x = lambda
  mu = 0.4345
  temp =  -1 + (x + 1/(4*x) + (3*mu)/x + (3*mu)/(2*x^2) )/sqrt(n)
  temp = temp - (mu/2)*(4 +  (1/x^2) )/n 
  temp = temp + (mu/2)*(4*x + (1/x) )*n^(-1.5)
  temp = (temp*2*x)/(3*sqrt(n)) + 1
  temp = temp*(exp(-2*x^2))
  temp = temp - alpha/2
  return(temp)
}

get_bds_dkw_complex = function(alpha = 0.05, n)
{
  lam_alpha =  uniroot( lambda_uppbd_func, lower =  1 , upper = 10, 
                        n = n , alpha = alpha, tol = 10^-10  )$root
  
  bds_dkw_complex = data.frame(low = rep(NA, n),
                               up = rep(NA, n))
  
  bds_dkw_complex$low = (1:n)/n - lam_alpha/sqrt(n)
  bds_dkw_complex$up = (1:n)/n  + lam_alpha/sqrt(n)
  
  bds_dkw_complex$low[ bds_dkw_complex$low < 0] = 0
  bds_dkw_complex$up[ bds_dkw_complex$up > 1] = 1
  
  return(bds_dkw_complex)
}

##########################################
##########################################

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
dw_quant = function( v = 3/2, alpha = 0.05, n_runs = 10000, n )
{
  sam = numeric(n_runs)
  t_n = (1:n)/(n) 
  t_n_shift = (0:(n-1))/n
  for( i in 1:n_runs)
  {
    U_sort = sort( runif(n))
    
    sam[i] = max( n*KL_dw(t_n_shift, U_sort) - C_v_dw(t_n_shift, U_sort,v = v),  
                  n*KL_dw(t_n, U_sort) - C_v_dw(t_n, U_sort,  v = v) ) 
  }
  
  n_cuttoff = ceiling(n_runs*(1-alpha))
  
  return(sort(sam)[n_cuttoff])
}



##########################################
##########################################
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

#####################

airf = read.table("/Users/sidsarkar/Desktop/Conf/bike/airfoil_self_noise.dat")
colnames(airf)[6] = "y"
n_train = 700
n_calib = 300
n_test = 503


shuffle_idx = sample(nrow(airf),rep = F)
airf_train = airf[shuffle_idx[1:n_train],]
airf_calib = airf[shuffle_idx[(n_train + 1):(n_train + n_calib)],]
airf_test = airf[shuffle_idx[- (1:(n_train + n_calib))],]

delta_n = 0.1
print("generating bounds")
bds_dw = get_bds_dw(alpha = delta_n, n = n_calib) 
bds_dkw_simp = get_bds_dkw_simp(alpha = delta_n, n = n_calib) 
bds_dkw_complex = get_bds_dkw_complex(alpha = delta_n, n = n_calib) 

airf_rqfit <- quantregForest( y = airf_train$y, x = dplyr::select(airf_train, -y))


################# quant on data calib
quant_ecdf = predict(airf_rqfit, newdata = dplyr::select(airf_calib, -y), 
                     what = ecdf)

score_n = numeric(n_calib)

for(i in 1:n_calib)
{
  score_n[i] = abs(quant_ecdf[[i]](airf_calib$y[i]) - 0.5)
}

################# quant on data train  
quant_ecdf_train = predict(airf_rqfit, newdata = dplyr::select(airf_train, -y) , 
                           what = ecdf)
score_n_overfit = numeric(n_train)

for(i in 1:n_train)
{
  score_n_overfit[i] = abs(quant_ecdf_train[[i]](airf_train$y[i]) - 0.5)
}

##### cutoff for CDFs
print("cutoffs and test data")
##### computing cutoffs

alpha_cov_overfit = numeric(19)
alpha_cov_calib = numeric(19)
alpha_cov_calib_dw = numeric(19)
alpha_cov_calib_dkw_simp = numeric(19)
alpha_cov_calib_dkw_complex = numeric(19)

alpha_width_overfit = numeric(19)
alpha_width_calib = numeric(19)
alpha_width_calib_dw = numeric(19)
alpha_width_calib_dkw_simp = numeric(19)
alpha_width_calib_dkw_complex = numeric(19)


for(idx in 2:19)
{  
  print(idx)
  alpha_n = idx/20
  
  ### Duembgun Wellner
  ### which point will give us the lower bound of F(x) >= 1- alpha
  low_idx_dw = min( which(bds_dw$low >= 1- alpha_n))
  F_s_cutoff_dw =  low_idx_dw/n_calib
  F_s_upperbd_dw = bds_dw$up[low_idx_dw] 
  
  s_cutoff_dw = quantile(score_n, probs = F_s_cutoff_dw)
  
  ### DKW simple
  ### which point will give us the lower bound of F(x) >= 1- alpha
  low_idx_dkw_simp = min( which(bds_dkw_simp$low >= 1- alpha_n))
  F_s_cutoff_dkw_simp =  low_idx_dkw_simp/n_calib
  F_s_upperbd_dkw_simp = bds_dkw_simp$up[low_idx_dkw_simp] 
  
  s_cutoff_dkw_simp = quantile(score_n, probs = F_s_cutoff_dkw_simp)
  
  ### DKW complex
  ### which point will give us the lower bound of F(x) >= 1- alpha
  low_idx_dkw_complex = min( which(bds_dkw_complex$low >= 1- alpha_n))
  F_s_cutoff_dkw_complex =  low_idx_dkw_complex/n_calib
  F_s_upperbd_dkw_complex = bds_dkw_complex$up[low_idx_dkw_complex] 
  
  s_cutoff_dkw_complex = quantile(score_n, probs = F_s_cutoff_dkw_complex)
  
  ### standard split conf
  s_cutoff = quantile(score_n, probs = 1- alpha_n)
  
  ### cutoff based on training data
  s_cutoff_overfit = quantile(score_n_overfit, probs = 1 - alpha_n)
  
  ######### 
  
  pred_overfit = predict(airf_rqfit, newdata = dplyr::select(airf_test, -y) ,
                         what = c( 0.5 - s_cutoff_overfit  , 0.5 + s_cutoff_overfit))
  pred_calib = predict(airf_rqfit, newdata = dplyr::select(airf_test, -y) ,
                       what = c( 0.5 - s_cutoff, 0.5 + s_cutoff ))
  pred_calib_dw = predict(airf_rqfit, newdata = dplyr::select(airf_test, -y) ,
                          what = c( 0.5 - s_cutoff_dw, 0.5 + s_cutoff_dw ))
  pred_calib_dkw_simp = predict(airf_rqfit, newdata = dplyr::select(airf_test, -y) ,
                                what = c( 0.5 - s_cutoff_dkw_simp, 0.5 + s_cutoff_dkw_simp ))
  pred_calib_dkw_complex = predict(airf_rqfit, newdata = dplyr::select(airf_test, -y) ,
                                   what = c( 0.5 - s_cutoff_dkw_complex, 0.5 + s_cutoff_dkw_complex ))
  
  cov_overfit = airf_test$y < pred_overfit[,1] | airf_test$y > pred_overfit[,2]
  cov_calib = airf_test$y < pred_calib[,1]  | airf_test$y > pred_calib[,2]
  cov_calib_dw = airf_test$y < pred_calib_dw[,1]  | airf_test$y > pred_calib_dw[,2]
  cov_calib_dkw_simp = airf_test$y < pred_calib_dkw_simp[,1]  | airf_test$y > pred_calib_dkw_simp[,2]
  cov_calib_dkw_complex = airf_test$y < pred_calib_dkw_complex[,1]| airf_test$y > pred_calib_dkw_complex[,2] 
  
  alpha_cov_overfit[idx] = 1- mean(cov_overfit)
  alpha_cov_calib[idx] = 1- mean(cov_calib)
  alpha_cov_calib_dw[idx] = 1- mean(cov_calib_dw)
  alpha_cov_calib_dkw_simp[idx] = 1- mean(cov_calib_dkw_simp)
  alpha_cov_calib_dkw_complex[idx] = 1- mean(cov_calib_dkw_complex)
  
  width_overfit = pred_overfit[,2] - pred_overfit[,1]
  width_calib = pred_calib[,2] - pred_calib[,1] 
  width_calib_dw = pred_calib_dw[,2] - pred_calib_dw[,1]
  width_calib_dkw_simp = pred_calib_dkw_simp[,2] - pred_calib_dkw_simp[,1]
  width_calib_dkw_complex = pred_calib_dkw_complex[,2] - pred_calib_dkw_complex[,1]
  
  alpha_width_overfit[idx] = mean(width_overfit)
  alpha_width_calib[idx] = mean(width_calib)
  alpha_width_calib_dw[idx] = mean(width_calib_dw)
  alpha_width_calib_dkw_simp[idx] = mean(width_calib_dkw_simp)
  alpha_width_calib_dkw_complex[idx] = mean(width_calib_dkw_complex)
}

df_dcp = data.frame( alpha = rep(1:19/20,  4), 
                     Method = rep(c("training data quantile", "split dcp","DKW univ", "DW univ"), each = 19),
                     coverage = rep(NA, 4*19),
                     width = rep(NA, 4*19),
                     base_coverage = rep(NA,4*19),
                     base_width = rep(NA,4*19))


df_dcp$coverage[ (1:19)] = alpha_cov_overfit
df_dcp$coverage[ (1:19) + (19)] = alpha_cov_calib  
df_dcp$coverage[ (1:19) + (38)] = alpha_cov_calib_dkw_simp
df_dcp$coverage[ (1:19) + (57)] = alpha_cov_calib_dw

df_dcp$width[ (1:19)] = alpha_width_overfit
df_dcp$width[ (1:19) + (19)] = alpha_width_calib  
df_dcp$width[ (1:19) + (38)] = alpha_width_calib_dkw_simp
df_dcp$width[ (1:19) + (57)] = alpha_width_calib_dw

df_dcp$base_width = rep(alpha_width_overfit,4)
df_dcp$base_coverage = rep(alpha_cov_overfit,4)


df_dcp = df_dcp %>% mutate(diff_coverage = coverage + alpha - 1)
df_dcp = df_dcp %>% mutate(diff_base_coverage = coverage - base_coverage)
df_dcp = df_dcp %>% mutate(diff_base_width = width - base_width)

df_dcp%>% filter(alpha> 1/20) %>%
  ggplot(mapping = aes(x = alpha, y = diff_coverage, color = Method)) +
  geom_line() + 
  geom_hline(yintercept = 0,linetype = 2, color = "grey") + 
  ylab(expression(paste("Coverage - (1 -", alpha, ")", sep=""))) + 
  ggtitle(expression(paste("Airfoil dataset: Coverage of methods with ",
                           delta," = 0.1" ,sep = "")),
          subtitle = expression(paste("(difference with (1-",alpha,"))", sep = ""))) +
  theme(legend.position="bottom")


df_dcp%>% filter(alpha> 1/20) %>%
  ggplot(mapping = aes(x = alpha, y = width, color = Method)) +
  geom_line() + 
  ylab("Width") + 
  ggtitle(expression(paste("Airfoil dataset: Width with ",
                           delta," = 0.1" ,sep = "")),
          subtitle = expression(paste("(difference with (1-",alpha,"))", sep = ""))) +
  theme(legend.position="bottom")

