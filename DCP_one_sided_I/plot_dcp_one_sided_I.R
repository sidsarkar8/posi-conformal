### for I,delta

library(ggplot2)
library(dplyr)

setwd("/Users/sidsarkar/Documents/Projects/Conf_redone/DCP_one_sided_I")

########################################################################

df_dcp_1 = 
  read.csv("/Users/sidsarkar/Documents/Projects/Conf_redone/DCP_one_sided_I/dcp_files/df_dcp_I_large.csv")[,-1] %>% as.data.frame()

alpha_seq = seq(from = 0.02, to = 0.5, length.out = 19)

####################################################################

M = 500

method_names = df_dcp_1$Method %>% unique()


df_prop_cross = data.frame(Method = rep(method_names, each = M), simul_id = rep(1:M, 4))
df_prop_cross = cbind(df_prop_cross, matrix(nrow = 4*M, ncol = 19))

df_prop_cross_full =  data.frame(Method = rep(method_names, each = M), simul_id = rep(1:M, 4),
                                 flag_if_cross = rep(NA,4*M))

for( idx in 1:(4*M))
{
  if(idx%%100 == 0){cat(' ',idx)}
  temp_method_name = df_prop_cross$Method[idx]
  temp_simul_id = df_prop_cross$simul_id[idx]
  
  temp_df_diff_cov = df_dcp_1 %>% filter(Method == temp_method_name, 
                                         simul_id == temp_simul_id ) %>% 
    dplyr::select(diff_coverage)
  df_prop_cross[idx,3:21] = as.numeric(temp_df_diff_cov >= 0)
  df_prop_cross_full$flag_if_cross[idx] = prod(temp_df_diff_cov>=0)
}

df_prop_cross_full =  df_prop_cross_full %>% group_by(Method) %>% 
  summarise( prop_if_cross_full := 1- mean(flag_if_cross) )

####################################################

colnames(df_prop_cross)[3:21] = paste("flag_al",1:19, sep="") 

df_prop_cross_al = df_prop_cross %>% group_by(Method) %>% 
  summarise( prop_if_cross1 := 1- mean(flag_al1) )

for(idx in 2:19)
{
  temp_name = paste("prop_if_cross", idx, sep="")
  temp_al = paste("flag_al", idx ,sep = "")
  df_prop_cross_al = cbind(df_prop_cross_al, 
                           df_prop_cross %>% group_by(Method) %>% summarise( {{temp_name}} := 1- mean(.data[[temp_al]]) )%>% .[,-1])
}


df_prop_cross_full = df_prop_cross_full %>% mutate(diff_coverage = -0.025, alpha = 0.2)

########################################################################
# ann_text = paste("# of curves above 1-alpha =" ,
#                  (1- df_prop_cross_full$prop_if_cross_full)*M)

ann_text = paste("prop. of curves above 1-alpha =" ,
                 (1- df_prop_cross_full$prop_if_cross_full))

df_dcp_1$Method[df_dcp_1$Method == "DKW univ"] = "DKW simul"
df_dcp_1$Method[df_dcp_1$Method == "DW univ"] = "DW simul"

df_prop_cross_full$Method[df_prop_cross_full$Method == "DKW univ"] = "DKW simul"
df_prop_cross_full$Method[df_prop_cross_full$Method == "DW univ"] = "DW simul"


df_dcp_1 %>% 
  ggplot(mapping = aes(x = alpha, y = diff_coverage, color = Method)) +
  geom_line(aes(group = simul_id),alpha = 0.3, 
            show.legend = FALSE) +
  facet_grid(~ Method) + 
  ylim(-0.03,0.03) +
  scale_x_continuous(breaks = alpha_seq) +
  geom_hline(yintercept = 0,linetype = 2) + 
  ylab(expression(paste("Coverage - (1 -", alpha, ")", sep=""))) + 
  ggtitle(expression(paste("Coverage of methods with ",
                           delta," = 0.1" ,sep = "")),
          subtitle = expression(paste("(difference with (1-",alpha,"))", sep = ""))) +
  theme(text = element_text(size = 13)) +
  geom_text(data = df_prop_cross_full, label =  ann_text, 
            col = "black")

# ggsave("./dcp_diff_coverage_multiple_1.pdf", width = 12, height = 3.5, 
#        dpi = 300, units = "in" )

####################################################################

df_dcp_1 %>% 
  ggplot(mapping = aes(x = alpha, y = ratio_width, color = Method)) +
  geom_line(aes(group = simul_id), alpha = 0.3, 
            show.legend = FALSE) +
  facet_grid(~ Method) + 
  geom_hline(yintercept = 1,linetype = 2) + 
  scale_x_continuous(breaks = alpha_seq) +
  ylab("Width / Width(true quantiles)") + 
  ggtitle(expression(paste("Width of prediction set across methods with ",
                           delta," = 0.1", sep = "" )),
          subtitle = "(ratio to true quantile widths)") +
  theme(text = element_text(size = 13)) 

# ggsave("./dcp_diff_coverage_multiple_one_sided.pdf", width = 12, height = 3.5, 
#        dpi = 300, units = "in" )

####################################################################

df_dcp_1 %>%
  ggplot(mapping = aes(x = alpha, y = ratio_width, color = Method)) +
  geom_line(aes(group = simul_id), alpha = 0.3,
            show.legend = FALSE) +
  facet_grid(~ Method) +
  geom_hline(yintercept = 1,linetype = 2) +
  scale_x_continuous(breaks = alpha_seq) +
  ylab("Width / Width(true quantiles)") +
  ggtitle(expression(paste("Width of prediction set across methods with ",
                           delta," = 0.1", sep = "" )),
          subtitle = "(ratio to true quantile widths)") +
  theme(text = element_text(size = 13))

# ggsave("./dcp_ratio_width_multiple_one_sided.pdf", width = 12, height = 3.5, 
#        dpi = 300, units = "in" )

####################################################################
