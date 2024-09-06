topDir <- "~/IMIC/study_case/table/predicted_growth/"
e_data <- read.table(paste0(topDir, 'equal_results_table.csv'),
                     header = T, sep = ',')
l_data <- read.table(paste0(topDir, 'less_results_table.csv'),
                     header = T, sep = ',')
m_data <- read.table(paste0(topDir, 'more_results_table.csv'),
                     header = T, sep = ',')

topDir <- "~/IMIC/study_case/table/"
qPCR <- read.table(paste0(topDir, 'qPCR-data.csv'),
                   header = T, sep = ',')

##### Organize data from qPCR data #####
equal <- qPCR[which(qPCR$condition == 'equal'),]
less <- qPCR[which(qPCR$condition == 'less'),]
more <- qPCR[which(qPCR$condition == 'more'),]

time <- c(0,4,8,24)
N_equal <- data.frame()
N_less <- data.frame()
N_more <- data.frame()

for (i in 1:length(time)) {
  subset_data <- equal$quantity[equal$species == 'EC' & equal$time == time[i]]
  N_equal[i,1] <- mean(log10(subset_data))
  N_equal[i,2] <- sd(log10(subset_data))
  subset_data_1 <- equal$quantity[equal$species == 'PP' & equal$time == time[i]]
  N_equal[4+i,1] <- mean(log10(subset_data_1))
  N_equal[4+i,2] <- sd(log10(subset_data_1))
  
  subset_data <- less$quantity[less$species == 'EC' & less$time == time[i]]
  N_less[i,1] <- mean(log10(subset_data))
  N_less[i,2] <- sd(log10(subset_data))
  subset_data_1 <- less$quantity[less$species == 'PP' & less$time == time[i]]
  N_less[4+i,1] <- mean(log10(subset_data_1))
  N_less[4+i,2] <- sd(log10(subset_data_1))
  
  subset_data <- more$quantity[more$species == 'EC' & more$time == time[i]]
  N_more[i,1] <- mean(log10(subset_data))
  N_more[i,2] <- sd(log10(subset_data))
  subset_data_1 <- more$quantity[more$species == 'PP' & more$time == time[i]]
  N_more[4+i,1] <- mean(log10(subset_data_1))
  N_more[4+i,2] <- sd(log10(subset_data_1))
}

merge_measured <- rbind(N_equal,N_less,N_more)
colnames(merge_measured) <- c('quantity','SD')
merge_measured$Species <- rep(c(rep('EC',4),rep('PP',4)),3)
merge_measured$Time <- as.factor(rep(c(rep(c('0h','4h','8h','24h'),2)),3))

###### organize the data of predicted growth rate ######
time <- c('0h','4h','8h','24h')
N_e_data <- data.frame()
N_l_data <- data.frame()
N_m_data <- data.frame()

for (i in 1:length(time)) {
  subset_data <- e_data$Growth_rate[e_data$Species == 'EC' & e_data$Timepoint == time[i]]
  N_e_data[i,1] <- mean(subset_data)
  N_e_data[i,2] <- sd(subset_data)
  subset_data_1 <- e_data$Growth_rate[e_data$Species == 'PP' & e_data$Timepoint == time[i]]
  N_e_data[4+i,1] <- mean(subset_data_1)
  N_e_data[4+i,2] <- sd(subset_data_1)
  
  subset_data <- l_data$Growth_rate[l_data$Species == 'EC' & l_data$Timepoint == time[i]]
  N_l_data[i,1] <- mean(subset_data)
  N_l_data[i,2] <- sd(subset_data)
  subset_data_1 <- l_data$Growth_rate[l_data$Species == 'PP' & l_data$Timepoint == time[i]]
  N_l_data[4+i,1] <- mean(subset_data_1)
  N_l_data[4+i,2] <- sd(subset_data_1)
  
  subset_data <- m_data$Growth_rate[m_data$Species == 'EC' & m_data$Timepoint == time[i]]
  N_m_data[i,1] <- mean(subset_data)
  N_m_data[i,2] <- sd(subset_data)
  subset_data_1 <- m_data$Growth_rate[m_data$Species == 'PP' & m_data$Timepoint == time[i]]
  N_m_data[4+i,1] <- mean(subset_data_1)
  N_m_data[4+i,2] <- sd(subset_data_1)
}

merge_data <- rbind(N_e_data,N_l_data,N_m_data)
colnames(merge_data) <- c('Mean', 'SD_1')
merge_data$Species <- rep(c(rep('EC',4),rep('PP',4)),3)

EC <- cbind(merge_measured[which(merge_measured$Species == 'EC'),], merge_data[which(merge_data$Species == 'EC'),1:2])
PP <- cbind(merge_measured[which(merge_measured$Species == 'PP'),], merge_data[which(merge_data$Species == 'PP'),1:2])

library(ggplot2)
library(ggpubr)

a <- ggscatter(data = EC, x =  'quantity', y = 'Mean',
               ylab = 'Growth Rate [1/h]',
               xlab = 'log(cell quantity)',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               color = '#F8766D',
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  font("ylab", size = 15)+
  font("xlab", size = 15)+
  labs(title = 'a.')+
  font("title", size = 15, face = 'bold')+
  font("x.text", size = 14)+
  font("y.text", size = 14)+
  stat_cor(method = "spearman", label.x = 6, label.y = 11, size= 6)+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA))+
  geom_errorbar(aes(ymin = Mean-SD_1, ymax= Mean+SD_1), color = '#F8766D', width=.1)+
  geom_errorbar(aes(xmin = quantity-SD, xmax= quantity+SD), color = '#F8766D', width=.1)+
  annotate('text', x = 5.5, y = 12, label = 'E. coli', size = 7, fontface = 'italic')

a

b<- ggscatter(data = PP, x =  'quantity', y = 'Mean',
              ylab = '',
              xlab = 'log(cell quantity)',
              add = "reg.line",                                 # Add regression line
              conf.int = TRUE,                                  # Add confidence interval
              color = '#00BFC4',
              add.params = list(color = "blue",
                                fill = "lightgray"))+
  font("xlab", size = 15)+
  labs(title = 'b.')+
  font("title", size = 15, face = 'bold')+
  font("x.text", size = 14)+
  font("y.text", size = 14)+
  stat_cor(method = "spearman", label.x = 6, label.y = 18, size= 6)+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA))+
  geom_errorbar(aes(ymin = Mean-SD_1, ymax= Mean+SD_1), color = '#00BFC4', width=.1)+
  geom_errorbar(aes(xmin = quantity-SD, xmax= quantity+SD), color = '#00BFC4', width=.1)+
  annotate('text', x = 5.8, y = 22, label = 'P. putida', size = 7, fontface = 'italic')

b

library(patchwork)
a+b+ plot_layout(nrow = 1)

ggsave(filename = '~/IMIC/Figure/Fig S5.svg',
       width = 22,height = 15,units = "cm")
