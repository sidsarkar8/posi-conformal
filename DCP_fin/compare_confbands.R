library(ggplot2)
library(dplyr)
library(tidyverse)
library(ggpubr)
#library(Rmpfr)

###################### DKW ######################
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

###################### DW ######################


source("/Users/sidsarkar/Documents/Projects/Conf_redone/dumbgun_wellner_new/dw_quantile_gen1.R")
source("/Users/sidsarkar/Documents/Projects/Conf_redone/dumbgun_wellner_new/get_bds_dw1.R")

###################### Asymptotic Eicker ######################
log2 = function(x)
{return(log(log(x)))}
log3 = function(x)
{return(log(log(log(x))))}


get_bds_eic = function(alpha = 0.05, n)
{
  bds_eic = data.frame(low = rep(NA, n),
                       up = rep(NA, n))
  
  t_eic = log( -2/(sqrt(pi)*log(1-alpha)) )
  
  t_eic_n = (2*log2(n))^(1/2)*( 1 + ((log3(n) + 2*t_eic)/(4*log2(n))) )
  
  bds_eic$low = (1:n)/n -  t_eic_n*sqrt((1:n)*(n - (1:n)))/n^(1.5)
  bds_eic$up = (1:n)/n + t_eic_n*sqrt((1:n)*(n - (1:n)))/n^(1.5)
  
  bds_eic$low[1] = 0
  bds_eic$low[n] = bds_eic$low[n-1]
  
  bds_eic$up[1] = bds_eic$up[2]
  bds_eic$up[n] = 1
  
  
  bds_eic$low[ bds_eic$low < 0] = 0
  bds_eic$up[ bds_eic$up > 1] = 1
  
  return(bds_eic)
}

###################### EIC with all indices ######################

get_bds_eic1 = function(alpha = 0.05, n)
{
  bds_eic = data.frame(low = rep(NA, n),
                       up = rep(NA, n))
  
  t_eic = log( -2/(sqrt(pi)*log(1-alpha)) )
  
  t_eic_n = (2*log2(n))^(1/2)*( 1 + ((log3(n) + 2*t_eic)/(4*log2(n))) )
  
  bds_eic$low = (1:n)/(n+1) -  t_eic_n*sqrt( (1:n)*(n + 1 - (1:n)))*(n+1)^(-1)*(n+2)^(-0.5)
  bds_eic$up = (1:n)/(n+1) + t_eic_n*sqrt( (1:n)*(n + 1 - (1:n)))*(n+1)^(-1)*(n+2)^(-0.5)
  
  bds_eic$low[ bds_eic$low < 0] = 0
  bds_eic$up[ bds_eic$up > 1] = 1
  
  return(bds_eic)
}

###################### Asymptotic AD ######################
fin_func_ad = function(p,q,quant,n)
{
  if(p==q){ 
    temp = -quant
  }else{
    temp = abs( sqrt(n)*(p-q)/(sqrt(q)*sqrt(1-q))) - quant 
  }
  return(temp)
}

get_bds_ad = function(alpha = 0.05, n)
{
  bds_ad = data.frame(low = rep(NA, n),
                      up = rep(NA, n))
  
  t_ad = log( -2/(sqrt(pi)*log(1-alpha)) )
  
  t_ad_n = (2*log2(n))^(1/2)*( 1 + ((log3(n) + 2*t_ad)/(4*log2(n))) )
  
  t_n = (1:n)/(n) 
  
  for( i in 2:(n-1))
  {
    ### lower bound
    if( fin_func_ad(p = t_n[i], q = 10^(-10), quant = t_ad_n, n = n) < 0 ){
      bds_ad$low[i] = 0
    }else{
      bds_ad$low[i] = uniroot( fin_func_ad, lower = 10^(-10) , upper = t_n[i], 
                               p = t_n[i], quant = t_ad_n, n = n, tol = 10^-10  )$root
    }
    
    ### upper bound
    if( fin_func_ad(p = t_n[i], q = 1-10^(-10), quant = t_ad_n, n = n) < 0 ){
      bds_ad$up[i] = 1
    }else{
      bds_ad$up[i] = uniroot( fin_func_ad, lower = t_n[i] , upper = 1-10^(-10),
                              p = t_n[i], quant = t_ad_n, n = n, tol = 10^-10 )$root
    }
  }
  
  bds_ad$low[1] = 0
  bds_ad$low[n] = bds_ad$low[n-1]
  
  bds_ad$up[1] = bds_ad$up[2]
  bds_ad$up[n] = 1
  
  bds_ad$low[ bds_ad$low < 0] = 0
  bds_ad$up[ bds_ad$up > 1] = 1
  
  return(bds_ad)
}

###################### Simes ######################

get_bds_simes = function(alpha = 0.05, n , k = NA)
{
  if(is.na(k)){k = n/2}
  #if(is.na(k_up)){k_up = n/2}
  
  bds_simes = data.frame(low = rep(NA, n),
                         up = rep(NA, n))
  
  for (i in 1:n)
  {
    #print(i)
    if(i < k){
      bds_simes$up[n - i + 1] = 1
      bds_simes$low[i] = 0
    }else{
      simes_seq = (i:(i - k + 1)) / (n:(n - k + 1))
      temp_log_simes = ( sum(log(simes_seq))/k + log(alpha/2)/k )
      
      #bds_simes$low[i] = (alpha/2)^(1/k)*(prod(simes_seq))^(1/k)
      bds_simes$low[i] = as.double(exp(temp_log_simes))  
      
      #bds_simes$up[n - i + 1] = 1 - (alpha / 2) ^ (1 / k) * (prod(simes_seq))^(1/k)
      bds_simes$up[n - i + 1] = 1 - as.double( exp(temp_log_simes))
    }
  }
  
  bds_simes$low[bds_simes$low<0] = 0
  bds_simes$up[bds_simes$up>1] = 1
  
  return(bds_simes)
}
################################## bates

source("/Users/sidsarkar/Documents/Projects/Conf_redone/dumbgun_wellner_new/get_bds_bates.R")


##################################
n = 500

get_bds_dw_n = get_bds_dw(n = n, alpha = 0.1) %>% mutate("Empirical.CDF.values" = 1:n/n) %>% mutate(Method = "DW")
get_bds_dw_n = get_bds_dw_n  %>% rename("lower.bound" = "low",
                                        "upper.bound" = "up") 

get_bds_dkw_n = get_bds_dkw_simp(n = n, alpha = 0.1) %>% mutate("Empirical.CDF.values" = 1:n/n) %>% mutate(Method = "DKW")
get_bds_dkw_n = get_bds_dkw_n  %>% rename("lower.bound" = "low",
                                          "upper.bound" = "up") 

get_bds_eic_n = get_bds_eic1(n = n, alpha = 0.1) %>% mutate("Empirical.CDF.values" = 1:n/n) %>% mutate(Method = "Eicker (Asymp)")
get_bds_eic_n = get_bds_eic_n  %>% rename("lower.bound" = "low",
                                          "upper.bound" = "up") 

get_bds_ad_n = get_bds_ad(n = n, alpha = 0.1) %>% mutate("Empirical.CDF.values" = 1:n/n) %>% mutate(Method = "AD (Asymp)")
get_bds_ad_n = get_bds_ad_n  %>% rename("lower.bound" = "low",
                                        "upper.bound" = "up") 

get_bds_simes_n = get_bds_simes(n = n, alpha = 0.1, k = n/2) %>% mutate("Empirical.CDF.values" = 1:n/n) %>% mutate(Method = "Simes")
get_bds_simes_n = get_bds_simes_n  %>% rename("lower.bound" = "low",
                                              "upper.bound" = "up") 

get_bds_bates_n = get_bds_bates(n = n, alpha = 0.1) %>% mutate("Empirical.CDF.values" = 1:n/n) %>% mutate(Method = "Bates")
get_bds_bates_n = get_bds_bates_n  %>% rename("lower.bound" = "low",
                                        "upper.bound" = "up") 


get_bds_fin = rbind(get_bds_dw_n, get_bds_dkw_n, get_bds_eic_n, get_bds_ad_n, get_bds_simes_n, get_bds_bates_n) %>% as.data.frame()
#get_bds_fin = rbind(get_bds_dw_n, get_bds_dkw_n, get_bds_eic_n,get_bds_ad_n) %>% as.data.frame()


get_bds_fin = get_bds_fin %>% mutate(lower.centred = lower.bound -  Empirical.CDF.values)
get_bds_fin = get_bds_fin %>% mutate(upper.centred = upper.bound -  Empirical.CDF.values)

ggplot(data = get_bds_fin, aes(x = Empirical.CDF.values, y = lower.bound, colour = Method)) +
  geom_line() +theme_bw() + 
  geom_line(data = get_bds_fin, aes(x = Empirical.CDF.values , y = upper.bound, colour = Method)) +
  geom_point(aes(y =Empirical.CDF.values, x =Empirical.CDF.values), col = "black") +
  ylab("") +
  theme(legend.position="bottom")+
  ggtitle(expression( paste("Comparison of methods with ", delta," = 0.1, n = ",100, sep = "")))+
  theme(text = element_text(size = 15))

#ggsave("/Users/sidsarkar/Desktop/Conf/DCP_fin/compareMethod.pdf") 

get_bds_fin = get_bds_fin %>% filter(Method!= "Simes")

p1 = ggplot(data = get_bds_fin , aes(x = Empirical.CDF.values, y = lower.centred, colour = Method)) +
  geom_line(linewidth = 1) +
  theme_bw() + 
  geom_line(data = get_bds_fin , aes(x = Empirical.CDF.values , y = upper.centred, colour = Method),
            linewidth = 1) +
  ylab("centred") +
  theme(legend.position="bottom", legend.text = element_text(size=14)) +
  # ggtitle(expression( paste("Comparison of methods with ", delta," = 0.1, n = ",1000, sep = "")))+
  # theme(text = element_text(size = 15))+ 
  theme(text = element_text(size = 15))
#ggsave("/Users/sidsarkar/Desktop/Conf/DCP_fin/compareMethod_1000.pdf", width = 8, height = 6, dpi = 300, units = "in" )


p2 = ggplot(data = get_bds_fin %>% filter(Empirical.CDF.values < 0.15), aes(x = Empirical.CDF.values, y = lower.centred, colour = Method)) +
  geom_line(linewidth = 1) +
  theme_bw() + 
  geom_line(data = get_bds_fin %>% filter(Empirical.CDF.values < 0.15), aes(x = Empirical.CDF.values , y = upper.centred, colour = Method),
            linewidth = 1) +
  ylab("centred") +
  theme(legend.position="bottom", legend.text = element_text(size=14)) +
  # ggtitle(expression( paste("Comparison of methods with ", delta," = 0.1, n = ",1000, sep = "")),
  #         subtitle = "Zoomed in on values < 0.2")+
  theme(text = element_text(size = 15))
#ggsave("/Users/sidsarkar/Desktop/Conf/DCP_fin/compareMethod_zoom_1000.pdf", width = 8, height = 6, dpi = 300, units = "in" )

p21 = ggarrange(p1, 
                p2, 
                ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

p21 = annotate_figure(p21, top = text_grob(paste("Comparison of confidence bands (centred) for m = ", n),size = 18))

p21
ggsave(paste("/Users/sidsarkar/Documents/Projects/Conf_redone/DCP_fin/compareMethod_comb_",n,"_new.pdf",
             sep= ""), width = 12, height = 6, dpi = 300, units = "in" )


########################### 
#########comparison across different n#########
###########################
n_list = c(200,500,1000)
for(n1 in n_list)
{
  n = n1
  print(n1)
  get_bds_dw_n = get_bds_dw(n = n, alpha = 0.1) %>% mutate("Empirical.CDF.values" = 1:n/n) %>% mutate(Method = "DW")
  get_bds_dw_n = get_bds_dw_n  %>% rename("lower.bound" = "low",
                                          "upper.bound" = "up") 
  
  get_bds_dkw_n = get_bds_dkw_simp(n = n, alpha = 0.1) %>% mutate("Empirical.CDF.values" = 1:n/n) %>% mutate(Method = "DKW")
  get_bds_dkw_n = get_bds_dkw_n  %>% rename("lower.bound" = "low",
                                            "upper.bound" = "up") 
  
  get_bds_eic_n = get_bds_eic1(n = n, alpha = 0.1) %>% mutate("Empirical.CDF.values" = 1:n/n) %>% mutate(Method = "Eicker (Asymp)")
  get_bds_eic_n = get_bds_eic_n  %>% rename("lower.bound" = "low",
                                            "upper.bound" = "up") 
  
  get_bds_ad_n = get_bds_ad(n = n, alpha = 0.1) %>% mutate("Empirical.CDF.values" = 1:n/n) %>% mutate(Method = "AD (Asymp)")
  get_bds_ad_n = get_bds_ad_n  %>% rename("lower.bound" = "low",
                                          "upper.bound" = "up") 
  
  get_bds_simes_n = get_bds_simes(n = n, alpha = 0.1, k = n/2) %>% mutate("Empirical.CDF.values" = 1:n/n) %>% mutate(Method = "Simes")
  get_bds_simes_n = get_bds_simes_n  %>% rename("lower.bound" = "low",
                                                "upper.bound" = "up") 
  
  get_bds_bates_n = get_bds_bates(n = n, alpha = 0.1) %>% mutate("Empirical.CDF.values" = 1:n/n) %>% mutate(Method = "Bates")
  get_bds_bates_n = get_bds_bates_n  %>% rename("lower.bound" = "low",
                                                "upper.bound" = "up") 
  
  get_bds_fin = rbind(get_bds_dw_n, get_bds_dkw_n, get_bds_eic_n, 
                      get_bds_ad_n, get_bds_simes_n, get_bds_bates_n) %>% as.data.frame()
  #get_bds_fin = rbind(get_bds_dw_n, get_bds_dkw_n, get_bds_eic_n,get_bds_ad_n) %>% as.data.frame()
  
  
  get_bds_fin = get_bds_fin %>% mutate(lower.centred = lower.bound -  Empirical.CDF.values)
  get_bds_fin = get_bds_fin %>% mutate(upper.centred = upper.bound -  Empirical.CDF.values)
  
  if(n1 == 200)
  {
    p200 = ggplot(data = get_bds_fin %>% filter(Empirical.CDF.values < 0.2), aes(x = Empirical.CDF.values, y = lower.centred, colour = Method)) +
      geom_line() +
      geom_line(data = get_bds_fin %>% filter(Empirical.CDF.values < 0.2), aes(x = Empirical.CDF.values , y = upper.centred, colour = Method)) +
      ylab("centred") +
      theme(legend.position="bottom")+
      # ggtitle(expression( paste("Comparison of methods with ", delta," = 0.1, n = ",1000, sep = "")),
      #         subtitle = "Zoomed in on values < 0.2")+
      theme(text = element_text(size = 15))}
  
  if(n1 == 500){
    
    p500 = ggplot(data = get_bds_fin %>% filter(Empirical.CDF.values < 0.2), aes(x = Empirical.CDF.values, y = lower.centred, colour = Method)) +
      geom_line() +
      geom_line(data = get_bds_fin %>% filter(Empirical.CDF.values < 0.2), aes(x = Empirical.CDF.values , y = upper.centred, colour = Method)) +
      ylab("centred") +
      theme(legend.position="bottom")+
      # ggtitle(expression( paste("Comparison of methods with ", delta," = 0.1, n = ",1000, sep = "")),
      #         subtitle = "Zoomed in on values < 0.2")+
      theme(text = element_text(size = 15))
  }
  
  if(n1 == 1000){
    
    p1000 = ggplot(data = get_bds_fin %>% filter(Empirical.CDF.values < 0.2), aes(x = Empirical.CDF.values, y = lower.centred, colour = Method)) +
      geom_line() +
      geom_line(data = get_bds_fin %>% filter(Empirical.CDF.values < 0.2), aes(x = Empirical.CDF.values , y = upper.centred, colour = Method)) +
      ylab("centred") +
      theme(legend.position="bottom")+
      # ggtitle(expression( paste("Comparison of methods with ", delta," = 0.1, n = ",1000, sep = "")),
      #         subtitle = "Zoomed in on values < 0.2")+
      theme(text = element_text(size = 15))
  }
}

p_n_list = ggarrange(p200,p500, p1000 , 
                ncol=3, nrow=1, common.legend = TRUE, legend="bottom")
p_n_list = annotate_figure(p_n_list, top = text_grob(paste("Comparison of confidence bands (centred) for n = 200, 500, 1000"),size = 17))
p_n_list
#ggsave("/Users/sidsarkar/Desktop/Conf/DCP_fin/compareMethod_across_n.pdf", width = 12, height = 4, dpi = 300, units = "in" )


#ggsave("/Users/sidsarkar/Desktop/Conf/DCP_fin/compareMethod.pdf") 
# #################
# n_runs = 10000
# meth_flag = data.frame(method = get_bds_fin$Method %>% unique(), prop = rep(0,5))
# 
# pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
#                      max = n_runs, # Maximum value of the progress bar
#                      style = 3,    # Progress bar style (also available style = 1 and style = 2)
#                      width = 50,   # Progress bar width. Defaults to getOption("width")
#                      char = "=")
# for( i in 1:n_runs)
# {
#   U_samp = sort(runif(n))
#   for(j in 1:5)
#   {
#     temp_low = prod(U_samp > get_bds_fin %>% filter(Method == meth_flag$method[j]) %>% select(lower.bound))
#     temp_up = prod(U_samp < get_bds_fin %>% filter(Method == meth_flag$method[j]) %>% select(upper.bound))
#     meth_flag[j,"prop"] = meth_flag[j,"prop"] + (temp_up*temp_low) 
#   }
#   setTxtProgressBar(pb,i)
# }
# meth_flag$prop = meth_flag$prop/n_runs
# meth_flag
