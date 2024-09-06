library(ggplot2)
library(forcats)
library(ggpubr)
library(scales)

topDir <- "~/IMIC/table/sensitivity_result/"

######## Sensitivity analysis ########
data <- read.table(paste0(topDir, 'consensus_sensitivity_analysis.csv'),
                   header = T, sep = ',')

lambda <- c('L = 0.1','L = 0.5','L = 1','L = 2','L = 3','L = 4','L = 5',
            'L = 6','L = 7','L = 8','L = 9','L = 10','L = 12','L = 14',
            'L = 15','L = 16','L = 18','L = 20','L = 22','L = 24',
            'L = 25','L = 50','L = 75','L = 100')

new_data <- matrix(nrow = length(lambda), ncol = 2)
colnames(new_data) <- c('Mean', 'SD')

for (i in 1:length(lambda)) {
  subset_data <- data$Solutions[data$Lambda == lambda[i]]
  new_data[i,1] <- mean(subset_data)
  new_data[i,2] <- sd(subset_data)
}

new_data <- data.frame(new_data)
new_data$Lambda <- as.numeric(as.character(c('0.1','0.5','1','2','3','4','5',
                                             '6','7','8','9','10','12','14',
                                             '15','16','18','20','22','24',
                                             '25','50','75','100')))


######### IMIC parameter testing #############
testtopDir <- "~/IMIC/table/parameter_test/"

consensus <- read.table(paste0(testtopDir, 'consensus_lambda_test_results.csv'),
                        header = T, sep = ',')
carveme <- read.table(paste0(testtopDir, 'carveme_lambda_test_results.csv'),
                      header = T, sep = ',')
gapseq <- read.table(paste0(testtopDir, 'gapseq_lambda_test_results.csv'),
                     header = T, sep = ',')
kbase <- read.table(paste0(testtopDir, 'kbase_lambda_test_results.csv'),
                    header = T, sep = ',')

cor_result <- data.frame()

for (i in 4:length(consensus)){
  cor_r <- cor.test(consensus[,3],consensus[,i], method = 'spearman',exact = F)
  cor_result[i-3,1] <- cor_r$estimate
  cor_result[i-3,2] <- cor_r$p.value
}

num <- nrow(cor_result)


for (i in 4:length(carveme)){
  cor_r <- cor.test(carveme[,3],carveme[,i], method = 'spearman',exact = F)
  cor_result[num+i-3,1] <- cor_r$estimate
  cor_result[num+i-3,2] <- cor_r$p.value
}

num <- nrow(cor_result)


for (i in 4:length(gapseq)){
  cor_r <- cor.test(gapseq[,3],gapseq[,i], method = 'spearman',exact = F)
  cor_result[num+i-3,1] <- cor_r$estimate
  cor_result[num+i-3,2] <- cor_r$p.value
}

num <- nrow(cor_result)


for (i in 4:length(kbase)){
  cor_r <- cor.test(kbase[,3],kbase[,i], method = 'spearman',exact = F)
  cor_result[num+i-3,1] <- cor_r$estimate
  cor_result[num+i-3,2] <- cor_r$p.value
}

colnames(cor_result) <- c('cor_ab', 'p_ab')

cor_result$Lambda <- as.numeric(as.character(rep(c('0.1','0.5','1','2','3','4','5',
                                                   '6','7','8','9','10','12','14',
                                                   '15','16','18','20','22','24',
                                                   '25','50','75','100'),4)))

cor_result$Methods <- as.factor(c(rep('Consensus',24),rep('CarveMe',24),rep('gapseq',24),rep('KBase',24)))

# Define the scale of the secondary axis to match primary axis scale
primary_range <- range(new_data$Mean, na.rm = TRUE)
secondary_range <- c(-0.7, 1)

# Combine plots
a<- ggplot()+
  # First plot: Sensitivity analysis
  geom_line(data = new_data,aes(x=Lambda, y=Mean, group = 1), linetype = 'dashed')+
  geom_point(data = new_data,aes(x=Lambda, y=Mean, group = 1))+
  geom_errorbar(data = new_data,aes(x=Lambda, y=Mean, 
                                    ymin=Mean-SD, ymax=Mean+SD), width=.2,position=position_dodge(0.05))+
  # Second plot: Lambda testing
  geom_line(data = cor_result, aes(x = Lambda, y = rescale(cor_ab, from = secondary_range, to = primary_range), 
                                   color = Methods, group = Methods))+
  geom_point(data = cor_result, aes(x = Lambda, y = rescale(cor_ab, from = secondary_range, to = primary_range), 
                                    color = Methods, group = Methods))+
  labs(x = "", y = "Value of objective function") +
  scale_y_continuous(
    sec.axis = sec_axis(~rescale(. , from = primary_range, to = secondary_range), name = "Spearman correlation coefficient")
  ) +
  # Adding a vertical line at x = 12
  geom_vline(xintercept = 12, color = "darkred", size = 1, linetype = "solid") +
  labs(title = 'a.')+
  font("ylab", size = 16)+
  font("title", size = 16, face = 'bold')+
  font("y.text", size = 15)+
  font("legend.text", size = 12)+
  theme(legend.position = c(0.9,0.25),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = 'black', size = 1, fill = NA))

a

############## study case ###################
topDir <- "~/IMIC/study_case/table/parameter_test/"
e_data <- read.table(paste0(topDir, 'equal_lambda_testing.csv'),
                     header = T, sep = ',')
l_data <- read.table(paste0(topDir, 'less_lambda_testing.csv'),
                     header = T, sep = ',')
m_data <- read.table(paste0(topDir, 'more_lambda_testing.csv'),
                     header = T, sep = ',')
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
  subset_data_1 <- equal$quantity[equal$species == 'PP' & equal$time == time[i]]
  N_equal[4+i,1] <- mean(log10(subset_data_1))
  
  subset_data <- less$quantity[less$species == 'EC' & less$time == time[i]]
  N_less[i,1] <- mean(log10(subset_data))
  subset_data_1 <- less$quantity[less$species == 'PP' & less$time == time[i]]
  N_less[4+i,1] <- mean(log10(subset_data_1))
  
  subset_data <- more$quantity[more$species == 'EC' & more$time == time[i]]
  N_more[i,1] <- mean(log10(subset_data))
  subset_data_1 <- more$quantity[more$species == 'PP' & more$time == time[i]]
  N_more[4+i,1] <- mean(log10(subset_data_1))
}

colnames(N_equal) <- 'log-quantity'
N_equal$Species <- c(rep('EC',4),rep('PP',4))
N_equal$Time <- rep(c('0h','4h','8h','24h'),2)

colnames(N_less) <- 'log-quantity'
N_less$Species <- c(rep('EC',4),rep('PP',4))
N_less$Time <- rep(c('0h','4h','8h','24h'),2)

colnames(N_more) <- 'log-quantity'
N_more$Species <- c(rep('EC',4),rep('PP',4))
N_more$Time <- rep(c('0h','4h','8h','24h'),2)


##### Organize data from preidction data #####
lambda <- c('L = 0.1','L = 0.5','L = 1','L = 2','L = 3','L = 4','L = 5',
            'L = 6','L = 7','L = 8','L = 9','L = 10','L = 12','L = 14',
            'L = 15','L = 16','L = 18','L = 20','L = 22','L = 24',
            'L = 25','L = 50','L = 75','L = 100')

time <- c('0h','4h','8h','24h')

N_e_data <- data.frame()
N_l_data <- data.frame()
N_m_data <- data.frame()

for (i in 1:length(lambda)) {
  for (j in 1:length(time)){
    subset_data <- e_data$Growth_rate[e_data$Species == 'EC' & e_data$Time == time[j] & e_data$Lambda == lambda[i]]
    N_e_data[j,i] <- mean(subset_data)
    subset_data_1 <- e_data$Growth_rate[e_data$Species == 'PP' & e_data$Time == time[j] & e_data$Lambda == lambda[i]]
    N_e_data[4+j,i] <- mean(subset_data_1) 
    
    subset_data <- l_data$Growth_rate[l_data$Species == 'EC' & l_data$Time == time[j] & l_data$Lambda == lambda[i]]
    N_l_data[j,i] <- mean(subset_data)
    subset_data_1 <- l_data$Growth_rate[l_data$Species == 'PP' & l_data$Time == time[j] & l_data$Lambda == lambda[i]]
    N_l_data[4+j,i] <- mean(subset_data_1) 
    
    subset_data <- m_data$Growth_rate[m_data$Species == 'EC' & m_data$Time == time[j] & m_data$Lambda == lambda[i]]
    N_m_data[j,i] <- mean(subset_data)
    subset_data_1 <- m_data$Growth_rate[m_data$Species == 'PP' & m_data$Time == time[j] & m_data$Lambda == lambda[i]]
    N_m_data[4+j,i] <- mean(subset_data_1) 
  }
}


cor_result_1 <- data.frame()

for (i in 1:length(lambda)){
  cor_r <- cor.test(N_e_data[,i],N_equal[,1], method = 'spearman',exact = F)
  cor_result_1[i,1] <- cor_r$estimate
  cor_result_1[i,2] <- cor_r$p.value
}

num <- nrow(cor_result_1)
for (i in 1:length(lambda)){
  cor_r <- cor.test(N_l_data[,i],N_less[,1], method = 'spearman',exact = F)
  cor_result_1[num+i,1] <- cor_r$estimate
  cor_result_1[num+i,2] <- cor_r$p.value
}

num <- nrow(cor_result_1)
for (i in 1:length(lambda)){
  cor_r <- cor.test(N_m_data[,i],N_more[,1], method = 'spearman',exact = F)
  cor_result_1[num+i,1] <- cor_r$estimate
  cor_result_1[num+i,2] <- cor_r$p.value
}

colnames(cor_result_1) <- c('Cor', 'p_value')
cor_result_1$Lambda <-  as.numeric(as.character(rep(c('0.1','0.5','1','2','3','4','5',
                                                      '6','7','8','9','10','12','14',
                                                      '15','16','18','20','22','24',
                                                      '25','50','75','100'),3)))

cor_result_1$Ratio <- as.factor(c(rep('1:1',24),rep('1:1000',24),rep('1000:1',24)))

###### Sensitivity analysis ######
topDir <- "~/IMIC/study_case/table/sensitivity_result/"

equal <- read.table(paste0(topDir, 'equal_sensitivity_analysis.csv'),
                    header = T, sep = ',')
more <- read.table(paste0(topDir, 'more_sensitivity_analysis.csv'),
                   header = T, sep = ',')
less <- read.table(paste0(topDir, 'less_sensitivity_analysis.csv'),
                   header = T, sep = ',')

lambda <- c('L = 0.1','L = 0.5','L = 1','L = 2','L = 3','L = 4','L = 5',
            'L = 6','L = 7','L = 8','L = 9','L = 10','L = 12','L = 14',
            'L = 15','L = 16','L = 18','L = 20','L = 22','L = 24',
            'L = 25','L = 50','L = 75','L = 100')

new_equal <- matrix(nrow = length(lambda), ncol = 2)
colnames(new_equal) <- c('Mean', 'SD')

for (i in 1:length(lambda)) {
  subset_data <- equal$Solutions[equal$Lambda == lambda[i]]
  new_equal[i,1] <- mean(subset_data)
  new_equal[i,2] <- sd(subset_data)
}

new_more <- matrix(nrow = length(lambda), ncol = 2)
colnames(new_more) <- c('Mean', 'SD')

for (i in 1:length(lambda)) {
  subset_data <- more$Solutions[more$Lambda == lambda[i]]
  new_more[i,1] <- mean(subset_data)
  new_more[i,2] <- sd(subset_data)
}

new_less <- matrix(nrow = length(lambda), ncol = 2)
colnames(new_less) <- c('Mean', 'SD')

for (i in 1:length(lambda)) {
  subset_data <- less$Solutions[less$Lambda == lambda[i]]
  new_less[i,1] <- mean(subset_data)
  new_less[i,2] <- sd(subset_data)
}

new_data_1 <- data.frame(rbind(new_equal,new_less,new_more))
new_data_1$Lambda <-as.numeric(as.character(rep(c('0.1','0.5','1','2','3','4','5',
                                                  '6','7','8','9','10','12','14',
                                                  '15','16','18','20','22','24',
                                                  '25','50','75','100'),3)))

new_data_1$Ratio <- as.factor(c(rep('1:1',24),rep('1:1000',24),rep('1000:1',24)))

# Define the scale of the secondary axis to match primary axis scale
primary_range_1 <- range(new_data_1$Mean, na.rm = TRUE)
secondary_range_1 <- c(-0.5, 1)

# Combine plots
b <- ggplot()+
  # First plot: Sensitivity analysis
  geom_line(data = new_data_1,aes(x=Lambda, y = Mean, color = Ratio, group = Ratio), linetype = 'dashed')+
  geom_point(data = new_data_1,aes(x=Lambda, y = Mean, color = Ratio, group = Ratio))+
  geom_errorbar(data = new_data_1,aes(x=Lambda, y=Mean, 
                                      ymin=Mean-SD, ymax=Mean+SD, color = Ratio), width=.2,position=position_dodge(0.05))+
  # Second plot: Lambda testing
  geom_line(data = cor_result_1, aes(x = Lambda, y = rescale(Cor, from = secondary_range_1, to = primary_range_1), 
                                     color = Ratio, group = Ratio))+
  geom_point(data = cor_result_1, aes(x = Lambda, y = rescale(Cor, from = secondary_range_1, to = primary_range_1), 
                                      color = Ratio, group = Ratio))+
  labs(x = "Balancing factor (λ)", y = "Value of objective function") +
  scale_y_continuous(
    sec.axis = sec_axis(~rescale(. , from = primary_range_1, to = secondary_range_1), name = "Spearman correlation coefficient")
  ) +
  # Adding a vertical line at x = 12
  geom_vline(xintercept = 6, color = "darkred", size = 1, linetype = "solid") +
  labs(title = 'b.')+
  font("ylab", size = 16)+
  font("xlab", size = 16)+
  font("title", size = 16, face = 'bold')+
  font("y.text", size = 15)+
  font("x.text", size = 15)+
  font("legend.text", size = 12)+
  theme(legend.position = c(0.9,0.2),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.title = element_blank(),
        panel.background = element_blank(),
        panel.border = element_rect(color = 'black', size = 1, fill = NA))+
  scale_color_manual(values = c("#999999", "#E69F00", "lightblue2"))

b

library(patchwork)
a+b+ plot_layout(nrow = 2)

ggsave(filename = '~/IMIC/Figure/Fig 3.svg',
       width = 20,height = 20,units = "cm")
