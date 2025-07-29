library(ggplot2)
library(dplyr)

setwd("~/Desktop/Conf/DCP_fin")

df_dcp_1 = NULL
for( p_idx in 1:10)
{
df_dcp_temp = read.csv(paste("./dcp_files/df_dcp_",p_idx,".csv",sep=""))[,-1] %>% as.data.frame()
df_dcp_temp = df_dcp_temp %>% mutate(simul_id = simul_id + (p_idx-1)*100)
df_dcp_1 = rbind(df_dcp_1, df_dcp_temp)
}

########################################################################

df_dcp_1 %>% 
  ggplot(mapping = aes(x = alpha, y = diff_coverage, color = Method)) +
  geom_line(aes(group = simul_id),alpha = 0.3, 
            show.legend = FALSE) +
  facet_grid(~ Method) + 
  scale_x_continuous(breaks = seq(1/20, 19/20, by = 0.2)) +
  geom_hline(yintercept = 0,linetype = 2) + 
  ylab(expression(paste("Coverage - (1 -", alpha, ")", sep=""))) + 
  ggtitle(expression(paste("Coverage of methods with ",
                           delta," = 0.1" ,sep = "")),
          subtitle = expression(paste("(difference with (1-",alpha,"))", sep = ""))) +
  theme(text = element_text(size = 13)) 

#ggsave("./dcp_diff_coverage_multiple_1.pdf", width = 12, height = 3.5, dpi = 300, units = "in" )

####################################################################

df_dcp_1 %>% 
  ggplot(mapping = aes(x = alpha, y = ratio_width, color = Method)) +
  geom_line(aes(group = simul_id), alpha = 0.3, 
            show.legend = FALSE) +
  facet_grid(~ Method) + 
  geom_hline(yintercept = 1,linetype = 2) + 
  scale_x_continuous(breaks = seq(1/20, 19/20, by = 0.2)) +
  ylab("Width / Width(true quantiles)") + 
  ggtitle(expression(paste("Width of prediction set across methods with ",
                           delta," = 0.1", sep = "" )),
          subtitle = "(ratio to true quantile widths)") +
  theme(text = element_text(size = 13)) 

#ggsave("./dcp_ratio_width_multiple_1.pdf", width = 12, height = 3.5, dpi = 300, units = "in" )

####################################################################
####################################################################

M = 1000

method_names = df_dcp_1$Method %>% unique()

df_prop_cross = data.frame(Method = rep(method_names, each = M), simul_id = rep(1:M, 4))
df_prop_cross = cbind(df_prop_cross, matrix(nrow = 4*M, ncol = 19))

df_prop_cross_full =  data.frame(Method = rep(method_names, each = M), simul_id = rep(1:M, 4),
                                 flag_if_cross = rep(NA,4*M))

for( idx in 1:(4*M))
{
  temp_method_name = df_prop_cross$Method[idx]
  temp_simul_id = df_prop_cross$simul_id[idx]
  
  temp_df_diff_cov = df_dcp_1 %>% filter(Method == temp_method_name, simul_id == temp_simul_id ) %>% .$diff_coverage
  df_prop_cross[idx,3:21] = as.numeric( temp_df_diff_cov >= 0)
  df_prop_cross_full$flag_if_cross[idx] = prod(temp_df_diff_cov>=0)
}

df_prop_cross_full =  df_prop_cross_full %>% group_by(Method) %>% summarise( prop_if_cross_full := 1- mean(flag_if_cross) )

####################################################

colnames(df_prop_cross)[3:21] = paste("flag_al",1:19, sep="") 
df_prop_cross_al = df_prop_cross %>% group_by(Method) %>% summarise( prop_if_cross1 := 1- mean(flag_al1) )
df_prop_cross_sd = df_prop_cross %>% group_by(Method) %>% summarise( sd_if_cross1 := sd(flag_al1) )
for(idx in 2:19)
{
  temp_name = paste("prop_if_cross", idx, sep="")
  temp_al = paste("flag_al", idx ,sep = "")
  df_prop_cross_al = cbind(df_prop_cross_al, 
                           df_prop_cross %>% group_by(Method) %>% summarise( {{temp_name}} := 1- mean(.data[[temp_al]]) )%>% .[,-1])
  
  temp_name = paste("sd_if_cross", idx, sep="")
  df_prop_cross_sd = cbind(df_prop_cross_sd, 
                           df_prop_cross %>% group_by(Method) %>% summarise( {{temp_name}} := sd(.data[[temp_al]]) )%>% .[,-1])
  
}



####################################################################
source("~/Desktop/Conf/dumbgun_wellner/get_bds_dw.R")

####################################################################
alpha_n = 1:19/20
theory_bd_dkw = pnorm( sqrt(log(2/0.1)/(2*(alpha_n)*(1- alpha_n))) )
theory_bd_dkw_fs = 1 - (0.05)^(1/(4*alpha_n*(1- alpha_n)))
del_m_dkw = (1/10000 + sqrt( log(20)/(20000)))/3
theory_bd_dkw_fs1 = numeric(length(alpha_n))
for(i in 1:19)
{
  if(alpha_n[i]<1/2){
    alpha_temp = alpha_n[i]
    print("<")
  } else if(alpha_n[i]>1/2 + del_m_dkw){
    alpha_temp = alpha_n[i] - del_m_dkw
    print(">")
  }else{
    alpha_temp = 1/2
    print("=")
  }
  
  theory_bd_dkw_fs1[i] = 1 - (0.05)^(1/(4*alpha_temp*(1- alpha_temp))) 
}



plot( y = 1- df_prop_cross_al[1,-1] , x = alpha_n, ylim = c(0.95,1),
      ylab = "prop", xlab = "alpha", main = "DKW, delta = 0.1", type = "l",lty = 2)
points( y = 1- df_prop_cross_al[1,-1] , x = alpha_n, ylim = c(0.95,1))
points(y = theory_bd_dkw, x = alpha_n, col = "red", type = "l") 
points(y = theory_bd_dkw_fs1, x = alpha_n, col = "green", type = "l")

legend("bottomright", legend=c("asymp theory", "empirical", "fin sample using beta"),
       col=c("red", "black", "green"), lty=c(1,2,1), cex=0.8)


#################################################

quant_df = read.csv(file = "~/Desktop/Conf/dumbgun_wellner/quant_df_one_10000.csv")
quant_df_val  = as.data.frame(quant_df)[2] %>% as.numeric()

v = 3/2
theory_bd_dw = pnorm(sqrt( 2*( C_dw(1-alpha_n) + v*D_dw(1-alpha_n) + quant_df_val)))
theory_bd_dw_fs = 1 - exp( - C_dw(1-alpha_n) - v*D_dw(1-alpha_n) - quant_df_val )

bds_dw_1 =  get_bds_dw(alpha = 0.1, n = 10000, v = v, quant_val = quant_df_val)
del_m_dw = max(bds_dw_1$up-bds_dw_1$low)/2
theory_bd_dw_fs1 = numeric(length(alpha_n))

for(i in 1:19)
{
  if(alpha_n[i]<1/2){
    alpha_temp = alpha_n[i]
    print("<")
  } else if(alpha_n[i]>1/2 + del_m_dw){
    alpha_temp = alpha_n[i] - del_m_dw
    print(">")
  }else{
    alpha_temp = 1/2
    print("=")
  }
  
  theory_bd_dw_fs1[i] = 1 - exp( - C_dw(1-alpha_temp) - v*D_dw(1-alpha_temp) - quant_df_val )
}


plot( y = 1- df_prop_cross_al[2,-1] , x = alpha_n, ylim = c(0.95,1),
      ylab = "prop", xlab = "alpha", main = "DW, delta = 0.1", type = "l",lty = 2)
points( y = 1- df_prop_cross_al[2,-1] , x = alpha_n, ylim = c(0.95,1))
points(y = theory_bd_dw, x = alpha_n, col = "red", type = "l")
points(y = theory_bd_dw_fs1, x = alpha_n, col = "green", type = "l")
points(y = theory_bd_dw_fs, x = alpha_n, col = "yellow", type = "l")

legend("bottomright", legend=c("asymp theory", "empirical", "fin sample using beta"),
       col=c("red", "black", "green"), lty=c(1,2,1), cex=0.8)

# plot( y = 1- df_prop_cross_al[3,-1] , x = alpha_n)
# plot( y = 1- df_prop_cross_al[4,-1] , x = alpha_n)

###############

df_plot = data.frame(Method = rep(c("DKW univ","DW univ"), each = 38),
                     alpha = rep(alpha_n , 4),
                     delta= rep(rep(c("Empirical","Theoretical lower bound"),each = 19), 2),
                     value = c( 1 - ((df_prop_cross_al[1,-1]) %>% as.numeric()), theory_bd_dkw_fs1,
                                    1 - ((df_prop_cross_al[2,-1]) %>% as.numeric()), theory_bd_dw_fs1)
                     )

plt_fix_al = df_plot %>% 
  ggplot(mapping = aes(x = alpha, y = value, color = delta)) +
  geom_line(aes(group = delta)) +
  facet_grid(~ Method) + 
  scale_x_continuous(breaks = seq(1/20, 19/20, by = 0.2)) +
  ylab("probability") +
  ggtitle(expression(paste("P( Coverage ">="1-",alpha," ) for fixed ", alpha, " with ", delta, " = 0.1 ",  sep = "")),
          subtitle = expression(paste(" Number of simulations = 1000", sep ="" )) ) +
  theme(legend.position="bottom") +
  theme(text = element_text(size = 10)) 
plt_fix_al
#ggsave("./fixed_al_theory_bd.pdf", width = 6, height = 3.5, dpi = 300, units = "in",
#       plot = plt_fix_al)

#################################
#################################
f1 = function(x, y)
{
  if(x == "Empirical")
  {del = y}else{
    del = 0
  }
  return(del*(1-del))
}
  
df_plot_sd = data.frame(Method = rep(c("DKW univ","DW univ"), each = 38),
                         alpha = rep(alpha_n , 4),
                         delta= rep(rep(c("Empirical","Theoretical lower bound"),each = 19), 2),
                         value = c( 1 - ((df_prop_cross_al[1,-1]) %>% as.numeric()), theory_bd_dkw_fs1,
                                    1 - ((df_prop_cross_al[2,-1]) %>% as.numeric()), theory_bd_dw_fs1))
df_plot_sd = df_plot_sd %>% mutate(up = qbinom(0.9, prob = value, size = 1000)/1000)
df_plot_sd = df_plot_sd %>% mutate(low = qbinom(0.1, prob = value, size = 1000)/1000)
df_plot_sd = df_plot_sd %>% mutate(sd = sqrt(value*(1-value)/1000))


df_plot_sd[(df_plot_sd$delta != "Empirical"), "up"] = df_plot_sd[(df_plot_sd$delta != "Empirical"), "value"]
df_plot_sd[(df_plot_sd$delta != "Empirical"), "low"] = df_plot_sd[(df_plot_sd$delta != "Empirical"), "value"]
df_plot_sd[(df_plot_sd$delta != "Empirical"), "sd"] = 0
  ggplot() + 
    geom_line(df_plot_sd, 
              mapping = aes(x = alpha, y = up, group = delta), linetype = "dashed") +
    geom_line(df_plot_sd, 
              mapping = aes(x = alpha, y = value - sd, group = delta), linetype = "dashed") +
    geom_line(df_plot_sd, 
              mapping = aes(x = alpha, y = value, group = delta, color = delta)) +
    facet_grid(~ Method) + 
    scale_x_continuous(breaks = seq(1/20, 19/20, by = 0.2)) +
    ylab("probability") +
    ggtitle(expression(paste("P( Coverage ">="1-",alpha," ) for fixed ", alpha, " with ", delta, " = 0.1 ",  sep = "")),
            subtitle = expression(paste(" Number of simulations = 1000", sep ="" )) ) +
    theme(legend.position="bottom") +
    theme(text = element_text(size = 10)) 


