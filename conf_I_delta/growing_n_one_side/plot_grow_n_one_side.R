library(ggplot2)
library(dplyr)

setwd("~/Documents/Projects/Conf/conf_I_delta/growing_n_one_side")
df_dcp_100 = read.csv("server_output_10000/df_dcp_I_delta_100.csv")[,-1]
df_dcp_500 = read.csv("server_output_10000/df_dcp_I_delta_500.csv")[,-1]
df_dcp_1000 = read.csv("server_output_10000/df_dcp_I_delta_1000.csv")[,-1]
df_dcp_2000 = read.csv("server_output_10000/df_dcp_I_delta_2000.csv")[,-1]
df_dcp_5000 = read.csv("server_output_10000/df_dcp_I_delta_5000.csv")[,-1]
df_dcp_10000 = read.csv("server_output_10000/df_dcp_I_delta_10000.csv")[,-1]

df_dcp_100 = df_dcp_100 %>% filter(Method != "DKW") %>% filter(Method != "DW")
df_dcp_500 = df_dcp_500 %>% filter(Method != "DKW") %>% filter(Method != "DW")
df_dcp_1000 = df_dcp_1000 %>% filter(Method != "DKW") %>% filter(Method != "DW")
df_dcp_2000 = df_dcp_2000 %>% filter(Method != "DKW") %>% filter(Method != "DW")
df_dcp_5000 = df_dcp_5000 %>% filter(Method != "DKW") %>% filter(Method != "DW")
df_dcp_10000 = df_dcp_10000 %>% filter(Method != "DKW") %>% filter(Method != "DW")
########################################

M = max(df_dcp_100$simul_id)
alpha_fin = df_dcp_100$alpha %>% unique()
method_names = df_dcp_100$Method %>% unique() 

####################### avg comparisons

df_grow_n = df_dcp_100 %>% group_by(alpha,Method) %>%
  summarise(avg_cov = mean(coverage),
            avg_ratio_width = mean(ratio_width)) %>% mutate(n = 100)

df_grow_n = rbind(df_grow_n,  df_dcp_500 %>% group_by(Method,alpha) %>%
                    summarise(avg_cov = mean(coverage),
                              avg_ratio_width = mean(ratio_width)) %>% mutate(n = 500))

df_grow_n = rbind(df_grow_n,  df_dcp_1000 %>% group_by(Method,alpha) %>% 
                    summarise(avg_cov = mean(coverage),
                              avg_ratio_width = mean(ratio_width)) %>% mutate(n = 1000))

df_grow_n = rbind(df_grow_n,  df_dcp_2000 %>% group_by(Method,alpha) %>% 
                    summarise(avg_cov = mean(coverage),
                              avg_ratio_width = mean(ratio_width)) %>% mutate(n = 2000))

df_grow_n = rbind(df_grow_n,  df_dcp_5000 %>% group_by(Method,alpha) %>% 
                    summarise(avg_cov = mean(coverage),
                              avg_ratio_width = mean(ratio_width)) %>% mutate(n = 5000))

df_grow_n = rbind(df_grow_n,  df_dcp_10000 %>% group_by(Method,alpha) %>% 
                    summarise(avg_cov = mean(coverage),
                              avg_ratio_width = mean(ratio_width)) %>% mutate(n = 10000))


ggplot(df_grow_n) + geom_line(aes(x = alpha, y = avg_ratio_width, color = Method))+
    facet_wrap(~ n, nrow=2) +
  ylab("Width / Width(true quantiles)") + 
  ggtitle(expression(paste("Avg.Width of prediction set with growing calibration size, ",
                           delta," = 0.1", sep = "" )),
          subtitle = "(ratio to true quantile widths)") +
  geom_hline(yintercept = 1,linetype = 2) + 
  theme(text = element_text(size = 13)) + 
  theme(legend.position = "bottom")

ggsave("./grow_n_width_one_side.pdf", width = 12, height = 6, dpi = 300, units = "in" )


ggplot(df_grow_n) + geom_line(aes(x = alpha, y = avg_cov - (1-alpha), color = Method))+
  facet_wrap(~ n, nrow=2) + 
  geom_hline(yintercept = 0,linetype = 2) + 
  ylab(expression(paste("Avg.Coverage - (1 -", alpha, ")", sep=""))) + 
  ggtitle(expression(paste("Avg.Coverage of methods with growing calibration size, ",
                           delta," = 0.1" ,sep = "")),
          subtitle = expression(paste("(difference with (1-",alpha,"))", sep = ""))) +
  theme(text = element_text(size = 13)) + 
  theme(legend.position = "bottom")

ggsave("./grow_n_cov_one_side.pdf", width = 12, height = 6, dpi = 300, units = "in" )

################################################################################################
#################### total behaviour ###############################################################
################################################################################################
df_dcp = df_dcp_1000

df_prop_cross = data.frame(Method = rep(method_names, each = M),
                           simul_id = rep(1:M, 4),
                           flag_if_cross = rep(NA, 4*M))
df_prop_cross = cbind(df_prop_cross, matrix(nrow = 4*M, 
                                            ncol = length(alpha_fin)))

df_prop_cross_full =  data.frame(Method = rep(method_names, each = M),
                                 simul_id = rep(1:M, 4),
                                 flag_if_cross = rep(NA,4*M))

for( idx in 1:(4*M))
{
  temp_method_name = df_prop_cross$Method[idx]
  temp_simul_id = df_prop_cross$simul_id[idx]

  temp_df_diff_cov = df_dcp %>% filter(Method == temp_method_name,
                                       simul_id == temp_simul_id ) %>%
    dplyr::select(diff_coverage)

  df_prop_cross[idx,3 + (1:length(alpha_fin))] = as.numeric( temp_df_diff_cov >= 0)
  df_prop_cross_full$flag_if_cross[idx] = prod(temp_df_diff_cov>=0)
}

df_prop_cross_full =  df_prop_cross_full %>%
  group_by(Method) %>% summarise( prop_if_cross_full := 1- mean(flag_if_cross) )

####################################################

colnames(df_prop_cross)[-(1:3)] = paste("flag_al",1:length(alpha_fin), sep="")

df_prop_cross_al = df_prop_cross %>%
  group_by(Method) %>% summarise( prop_if_cross1 := 1- mean(flag_al1) )

for(idx in (1:length(alpha_fin)) )
{
  temp_name = paste("prop_if_cross", idx, sep="")
  temp_al = paste("flag_al", idx ,sep = "")
  df_prop_cross_al = cbind(df_prop_cross_al,
                           df_prop_cross %>% group_by(Method) %>% summarise( {{temp_name}} := 1- mean(.data[[temp_al]]) )%>% .[,-1])
}


df_prop_cross_full = df_prop_cross_full %>% mutate(diff_coverage = -0.025, alpha = 0.25)

########################################################################
# ann_text = paste("# of curves >= 1-alpha =" ,
#                  (1- df_prop_cross_full$prop_if_cross_full))
ann_text = paste("Prop. of curves >= 1-alpha :" ,
                 (1- df_prop_cross_full$prop_if_cross_full))


###################################################################
###################################################################

df_dcp %>%
  ggplot(mapping = aes(x = alpha, y = diff_coverage, color = Method)) +
  geom_line(aes(group = simul_id),alpha = 0.3,
             show.legend = FALSE) +
  facet_grid(~ Method) +
  scale_x_continuous(breaks = alpha_fin[(1:6)*8]) +
  geom_hline(yintercept = 0,linetype = 2) +
  ylab(expression(paste("Coverage - (1 -", alpha, ")", sep=""))) +
  ggtitle(expression(paste("Coverage of methods with ",
                           delta," = 0.1, n = 1000" ,sep = "")),
          subtitle = expression(paste("(difference with (1-",alpha,"))", sep = ""))) +
  theme(text = element_text(size = 13)) +
  geom_text(data = df_prop_cross_full, label =  ann_text,
            col = "black")

ggsave("./simul_cov_1000.pdf",
            width = 12, height = 3.9, dpi = 300, units = "in" )

###################################################################

# df_dcp %>%
#   ggplot(mapping = aes(x = alpha, y = ratio_width, color = Method)) +
#   geom_line(aes(group = simul_id), alpha = 0.3,
#             show.legend = FALSE) +
#   facet_grid(~ Method) +
#   geom_hline(yintercept = 1,linetype = 2) +
#   scale_x_continuous(breaks = alpha_fin[(1:6)*8]) +
#   ylab("width / width(true quantiles)") +
#   ggtitle(expression(paste("Width of prediction set across methods ",
#                            delta,"= 0.1", sep = "" )),
#           subtitle = "(ratio to true quantile widths)") +
#   theme(text = element_text(size = 13))

#ggsave("/home/siddhaarth/DCP_calib_fin/ratio_width2/dcp_ratio_width_1.pdf",
#     ratio_width = 12, height = 3.5, dpi = 300, units = "in" )


###################################################################
