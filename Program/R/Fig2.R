topDir <- "~/IMIC/table/predicted_growth/"
micom <- read.table(paste0(topDir, 'consensus_MICOM_results_table.csv'),
                    header = T, sep = ',')

library(ggplot2)
library(ggpubr)
micom <- micom %>% mutate_if(is.numeric, ~round(., 6))
ab1 <- ggscatter(micom, x =  'MAG_ab', y = 'MICOM_ab1',
                 xlab = 'Measured relative abundance (%)',
                 ylab = '')+
  labs(title = 'b.')+
  font("title", size = 15, face = 'bold')+
  font("xlab", size = 15)+
  font("x.text", size = 13)+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA))

ab1

ab <- ggscatter(micom, x =  'MAG_ab', y = 'MICOM',
                xlab = 'Measured relative abundance (%)',
                ylab = 'Predicted growth rate [1/h]',
                add = "reg.line",                                 # Add regression line
                conf.int = TRUE,                                  # Add confidence interval
                add.params = list(color = "blue",
                                  fill = "lightgray"))+
  ylim(0,150)+
  labs(title = 'a.')+
  font("ylab", size = 15)+
  font("xlab", size = 15)+
  font("title", size = 15, face = 'bold')+
  font("x.text", size = 13)+
  font("y.text", size = 13)+
  stat_cor(method = "spearman", label.x = 5, label.y = 140, size= 5.5, cor.coef.name = 'rho')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA))

ab

##########################################################################
library(forcats)
library(scales)

topDir <- "~/IMIC/table/parameter_test/"

##### coco test #####
data <- read.table(paste0(topDir, 'consensus_coco_test_alpha_0.3.csv'),
                   header = T, sep = ',')
data <- data %>% mutate_if(is.numeric, ~round(., 6))

##### coco test with abundance = 1 results #####
data_ab <- read.table(paste0(topDir, 'consensus_coco_test_with_abundance_1_alpha_0.3.csv'),
                      header = T, sep = ',')
data_ab <- data_ab %>% mutate_if(is.numeric, ~round(., 6))

cor_result <- data.frame()

# calculate the correlation
for (i in 4:length(data)){
  cor_r <- cor.test(data[,2],data[,i], method = 'spearman',exact = F)
  cor_result[i-3,1] <- cor_r$estimate
}

for (i in 4:length(data_ab)){
  if (sum(is.na(data_ab[, i])) != 70){
    cor_r <- cor.test(data_ab[,2],data_ab[,i], method = 'spearman',exact = F)
    cor_result[i-3,2] <- cor_r$estimate
  }
  else {
    cor_result[i-3,2] <- NA
  }
}

colnames(cor_result) <- c('cor', 'cor_ab')

# find the correlation coefficience value of corresponding delta
data_1 <- data[,-c(1:3)]
#colnames(data_1)[which(cor_result$cor == max(cor_result$cor))]

a <- cor_result[grepl('D1.G', colnames(data_1)),1]
a <- c(a,cor_result[grepl('D2.G', colnames(data_1)),1])
a <- c(a,cor_result[grepl('D3.G', colnames(data_1)),1])
a <- c(a,cor_result[grepl('D4.G', colnames(data_1)),1])    
a <- c(a,cor_result[grepl('D5.G', colnames(data_1)),1])
a <- c(a,cor_result[grepl('D6.G', colnames(data_1)),1])
a <- c(a,cor_result[grepl('D7.G', colnames(data_1)),1])
a <- c(a,cor_result[grepl('D8.G', colnames(data_1)),1])
a <- c(a,cor_result[grepl('D9.G', colnames(data_1)),1])
a <- c(a,cor_result[grepl('D10.G', colnames(data_1)),1])

data_2 <- data_ab[,-c(1:3)]
a <- c(a,cor_result[grepl('D1.G', colnames(data_2)),2])
a <- c(a,cor_result[grepl('D2.G', colnames(data_2)),2])
a <- c(a,cor_result[grepl('D3.G', colnames(data_2)),2])
a <- c(a,cor_result[grepl('D4.G', colnames(data_2)),2])    
a <- c(a,cor_result[grepl('D5.G', colnames(data_2)),2])
a <- c(a,cor_result[grepl('D6.G', colnames(data_2)),2])
a <- c(a,cor_result[grepl('D7.G', colnames(data_2)),2])
a <- c(a,cor_result[grepl('D8.G', colnames(data_2)),2])
a <- c(a,cor_result[grepl('D9.G', colnames(data_2)),2])
a <- c(a,cor_result[grepl('D10.G', colnames(data_2)),2])

delt <- data.frame(a)
colnames(delt) <- 'cor_value'

delt$delta  <- as.factor(rep(c(rep('10', 10),rep('20', 10),rep('30', 10),rep('40', 10),rep('50', 10),
                               rep('60', 10),rep('70', 10),
                               rep('80', 10),rep('90', 10),rep('100', 10)),2))
delt$abundance <- c(rep('With abundance',100),rep('Without abundance',100))

delta_p <- ggplot(delt, aes(x= fct_inorder(delta), y = cor_value, color = factor(abundance)))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(), alpha = 0.5)+
  theme_linedraw()+
  theme(legend.position = 'none')+
  ylim(0,1)+
  xlab('𝛿: scaling factor')+
  ylab('Spearman correlation coefficient')+
  labs(title = 'c.')+
  font("title", size = 15, face = 'bold')+
  font("ylab", size = 15)+ 
  font("xlab", size = 15)+
  font("y.text", size = 13)+
  font("x.text", size = 13)

delta_p

# find the correlation coefficience value of corresponding gamma
b <- cor_result[grepl('G1.p', colnames(data_1)),1]
b <- c(b,cor_result[grepl('G2.p', colnames(data_1)),1])
b <- c(b,cor_result[grepl('G3.p', colnames(data_1)),1])
b <- c(b,cor_result[grepl('G4.p', colnames(data_1)),1])    
b <- c(b,cor_result[grepl('G5.p', colnames(data_1)),1])
b <- c(b,cor_result[grepl('G6.p', colnames(data_1)),1])
b <- c(b,cor_result[grepl('G7.p', colnames(data_1)),1])
b <- c(b,cor_result[grepl('G8.p', colnames(data_1)),1])
b <- c(b,cor_result[grepl('G9.p', colnames(data_1)),1])
b <- c(b,cor_result[grepl('G10.p', colnames(data_1)),1])

b <- c(b,cor_result[grepl('G1.p', colnames(data_2)),2])
b <- c(b,cor_result[grepl('G2.p', colnames(data_2)),2])
b <- c(b,cor_result[grepl('G3.p', colnames(data_2)),2])
b <- c(b,cor_result[grepl('G4.p', colnames(data_2)),2])    
b <- c(b,cor_result[grepl('G5.p', colnames(data_2)),2])
b <- c(b,cor_result[grepl('G6.p', colnames(data_2)),2])
b <- c(b,cor_result[grepl('G7.p', colnames(data_2)),2])
b <- c(b,cor_result[grepl('G8.p', colnames(data_2)),2])
b <- c(b,cor_result[grepl('G9.p', colnames(data_2)),2])
b <- c(b,cor_result[grepl('G10.p', colnames(data_2)),2])


gamm <- data.frame(b)
colnames(gamm) <- 'cor_value'

gamm$gamma  <- as.factor(rep(c(rep('10', 10),rep('20', 10),rep('30', 10),rep('40', 10),rep('50', 10),
                               rep('60', 10),rep('70', 10),
                               rep('80', 10),rep('90', 10),rep('100', 10)),2))
gamm$abundance <- c(rep('With abundance',100),rep('Without abundance',100))

gamma_p <- ggplot(gamm, aes(x= fct_inorder(gamma), y = cor_value, color = factor(abundance)))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(), alpha = 0.5)+
  theme_linedraw()+
  theme(axis.text.y = element_blank(),
        legend.title = element_blank(),
        legend.position= 'inside')+
  ylim(0,1)+
  xlab('𝛾: regulating factor')+
  ylab('')+
  labs(title = 'd.')+
  font("xlab", size = 15)+
  font("title", size = 15, face = 'bold')+
  font("legend.text", size = 12)+
  font("x.text", size = 13)

gamma_p

############################################################################
consensus <- read.table(paste0(topDir, 'consensus_lambda_test_results.csv'),
                        header = T, sep = ',')
consensus <- consensus %>% mutate_if(is.numeric, ~round(., 6))

cor_result <- data.frame()

for (i in 4:length(consensus)){
  cor_r <- cor.test(consensus[,2],consensus[,i], method = 'spearman',exact = F)
  cor_result[i-3,1] <- cor_r$estimate
  cor_result[i-3,2] <- cor_r$p.value
}


colnames(cor_result) <- c('cor_ab', 'p_ab')

cor_result$Lambda <- as.factor(c('0.1','0.5','1','2','3','4','5',
                                 '6','7','8','9','10','12','14',
                                 '15','16','18','20','22','24',
                                 '25','50','75','100'))

library(forcats)
library(ggplot2)

e <- ggplot(cor_result, aes(x = fct_inorder(Lambda), y = cor_ab))+
  geom_col() +
  labs(x = "Balancing factor (λ)", y = "Spearman correlation coefficient") +
  theme_linedraw()+
  labs(title = 'e.')+
  font("title", size = 15, face = 'bold')+
  theme(axis.text.y = element_text(size = 15),
        axis.text.x = element_text(size = 15, angle = 90, hjust = 1),
        axis.title = element_text(size = 15))

e

######## consensus model ########
topDir <- "~/IMIC/table/predicted_growth/"
data <- read.table(paste0(topDir, 'consensus_results_table.csv'),
                   header = T, sep = ',')

library(ggplot2)
library(ggpubr)

f <- ggscatter(data, x =  'MAG_ab', y = 'IMIC',
               xlab = 'Measured relative abundance (%)',
               ylab = 'Predicted growth rate [1/h]',
               add = "reg.line",                               
               conf.int = TRUE,                               
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 'f.')+
  font("title", size = 15, face = 'bold')+
  font("y.text", size = 15)+
  font('ylab', size = 15)+
  font("x.text", size = 15)+
  font('xlab', size = 15)+
  stat_cor(method = "spearman", label.x = 20, label.y = 40, size= 5.5, cor.coef.name = 'rho')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),)

f

library(patchwork)
ab + ab1 + delta_p+gamma_p+e+f+ plot_layout(nrow = 3)

ggsave(filename = '~/IMIC/Figure/Fig 2.svg',
       width = 25,height = 30,units = "cm")

