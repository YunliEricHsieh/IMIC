topDir <- "~/IMIC/table/predicted_growth/"

######## consensus model ########
data <- read.table(paste0(topDir, 'consensus_results_table.csv'),
                   header = T, sep = ',')

library(ggplot2)
library(ggpubr)

a <- ggscatter(data, x =  'MAG_ab', y = 'IMIC',
               xlab = '',
               ylab = 'Predicted growth rate [1/h]',
               add = "reg.line",                               
               conf.int = TRUE,                               
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 'a. consensus')+
  font("title", size = 15, face = 'bold')+
  font("y.text", size = 15)+
  font('ylab', size = 15)+
  stat_cor(method = "spearman", label.x = 20, label.y = 40, size= 5.5, cor.coef.name = 'rho')+
  theme(axis.text.x = element_blank())

a

######## carveme model ########

data <- read.table(paste0(topDir, 'carveme_results_table.csv'),
                   header = T, sep = ',')

b <- ggscatter(data, x =  'MAG_ab', y = 'IMIC',
               xlab = '',
               ylab = '',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 'b. CarveMe')+
  font("title", size = 15, face = 'bold')+
  font("y.text", size = 15)+
  stat_cor(method = "spearman", label.x = 40, label.y = 35, size= 5.5, cor.coef.name = 'rho')+
  theme(axis.text.x = element_blank())

b

######## gapseq models ########
data <- read.table(paste0(topDir, 'gapseq_results_table.csv'),
                   header = T, sep = ',')

c <- ggscatter(data, x =  'MAG_ab', y = 'IMIC',
               xlab = 'Measured relative abundance (%)',
               ylab = 'Predicted growth rate [1/h]',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 'c. gapseq')+
  font("title", size = 15, face = 'bold')+
  font("y.text", size = 15)+
  font("x.text", size = 15)+
  font("xlab", size = 15)+
  font("ylab", size = 15)+
  stat_cor(method = "spearman", label.x = 30, label.y = 30, size= 5.5, cor.coef.name = 'rho')

c

######## kbase models ########
data <- read.table(paste0(topDir, 'kbase_results_table.csv'),
                   header = T, sep = ',')

d <- ggscatter(data, x =  'MAG_ab', y = 'IMIC',
               xlab = 'Measured relative abundance (%)',
               ylab = '',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 'd. KBase')+
  font("xlab", size = 15)+
  font("title", size = 15, face = 'bold')+
  font("x.text", size = 15)+
  font("y.text", size = 15)+
  stat_cor(method = "spearman", label.x = 40, label.y = 30, size= 5.5, cor.coef.name = 'rho')

d

library(patchwork)
a+b+c+d+plot_layout(nrow = 2)


ggsave(filename = "~/IMIC/Figure/Fig S1.svg",
       width = 25,height = 20,units = "cm")
