topDir <- "~/IMIC/table/"
consensus <- read.table(paste0(topDir, 'predicted_growth/consensus_results_table.csv'),
                   header = T, sep = ',')
micom.consensus <- read.table(paste0(topDir, 'predicted_growth/consensus_MICOM_results_table.csv'),
                              header = T, sep = ',')
coco.consensus <- read.table(paste0(topDir, 'parameter_test/consensus_coco_test_alpha_0.3.csv'),
                        header = T, sep = ',')
coco.consensus.ab1 <- read.table(paste0(topDir, 'parameter_test/consensus_coco_test_with_abundance_1_alpha_0.3.csv'),
                             header = T, sep = ',')
carveme <- read.table(paste0(topDir, 'predicted_growth/carveme_results_table.csv'),
                        header = T, sep = ',')
micom.carveme <- read.table(paste0(topDir, 'predicted_growth/carveme_MICOM_results_table.csv'),
                              header = T, sep = ',')
coco.carveme <- read.table(paste0(topDir, 'parameter_test/carveme_coco_test_alpha_0.5.csv'),
                             header = T, sep = ',')
coco.carveme.ab1 <- read.table(paste0(topDir, 'parameter_test/carveme_coco_test_with_abundance_1_alpha_0.5.csv'),
                                 header = T, sep = ',')
gapseq <- read.table(paste0(topDir, 'predicted_growth/gapseq_results_table.csv'),
                        header = T, sep = ',')
micom.gapseq <- read.table(paste0(topDir, 'predicted_growth/gapseq_MICOM_results_table.csv'),
                              header = T, sep = ',')
coco.gapseq <- read.table(paste0(topDir, 'parameter_test/gapseq_coco_test_alpha_0.5.csv'),
                           header = T, sep = ',')
coco.gapseq.ab1 <- read.table(paste0(topDir, 'parameter_test/gapseq_coco_test_with_abundance_1_alpha_0.5.csv'),
                               header = T, sep = ',')
kbase <- read.table(paste0(topDir, 'predicted_growth/kbase_results_table.csv'),
                        header = T, sep = ',')
micom.kbase <- read.table(paste0(topDir, 'predicted_growth/kbase_MICOM_results_table.csv'),
                              header = T, sep = ',')
coco.kbase <- read.table(paste0(topDir, 'parameter_test/kbase_coco_test_alpha_0.5.csv'),
                          header = T, sep = ',')
coco.kbase.ab1 <- read.table(paste0(topDir, 'parameter_test/kbase_coco_test_with_abundance_1_alpha_0.5.csv'),
                              header = T, sep = ',')

library(ggplot2)
library(ggpubr)
library(dplyr)

# IMIC compare with replication rate
# consensus - IMIC
# Round only numeric columns to 8 decimal places
consensus <- consensus %>% mutate_if(is.numeric, ~round(., 6))
a <- ggscatter(consensus, x =  'Replication_rate', y = 'IMIC',
                xlab = '',
                ylab = 'Predicted growth rate [1/h]',
                add = "reg.line",                                 # Add regression line
                conf.int = TRUE,                                  # Add confidence interval
                add.params = list(color = "blue",
                                  fill = "lightgray"))+
  labs(title = 'a')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = max(consensus$IMIC), size= 7, cor.coef.name = 'rho')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        axis.text.x = element_blank())

a
# consensus - MICOM-ab1
micom.consensus <- micom.consensus %>% mutate_if(is.numeric, ~round(., 6))
b <- ggscatter(micom.consensus, x =  'Replication_rate', y = 'MICOM_ab1',
               xlab = '',
               ylab = '')+
  labs(title = 'b')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        axis.text.x = element_blank())+
  scale_y_continuous(
    limits = c(15.6, 15.9),
    breaks = seq(15.6, 15.9, by = 0.1),
    labels = number_format(accuracy = 0.1))

b
# consensus - coco-ab1
# calculate the correlation for different parameter values
coco.consensus.ab1 <- coco.consensus.ab1 %>% mutate_if(is.numeric, ~round(., 6))
cor_result <- data.frame()

for (i in 4:length(coco.consensus.ab1)){
  if (sum(is.na(coco.consensus.ab1[, i])) < 70){
    cor_r <- cor.test(coco.consensus.ab1[,3],coco.consensus.ab1[,i], method = 'spearman',exact = F)
    cor_result[i-3,1] <- cor_r$estimate
  }
  else {
    cor_result[i-3,1] <- NA
  }
}

colnames(cor_result) <- c('cor')

# calculate the mean correlation value and find the one most close to the mean value
closest_value <- cor_result$cor[which.min(abs(cor_result$cor - mean(cor_result$cor, na.rm = TRUE)))]

na <- colnames(coco.consensus.ab1)[which(cor_result$cor == closest_value)+3]

library(grid)

c <- ggscatter(coco.consensus.ab1, x =  'Replication_rate', y = na,
         xlab = '',
         ylab = '',
         add = "reg.line",                                 # Add regression line
         conf.int = TRUE,                                  # Add confidence interval
         add.params = list(color = "blue",
                 fill = "lightgray"))+
  labs(title = 'c')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 15.93, size= 7, cor.coef.name = 'rho')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
    axis.text.x = element_blank())+
  scale_y_continuous(
    limits = c(15.7, 16.0),
    breaks = seq(15.7, 16.0, by = 0.1),
    labels = number_format(accuracy = 0.1))

c

# Distribution of correlation
di <- ggplot(cor_result, aes(x = cor)) +
  geom_histogram(binwidth = 0.005, fill = "gray", color = "black") +
  labs(x = "", y = "") +
  font("y.text", size = 13)+
  font("x.text", size = 13)+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
    panel.background = element_blank(),
    panel.grid = element_blank())

# Combine the plots
c1 <- c + annotation_custom(ggplotGrob(di), xmin = 0.1, xmax = 0.8, ymin = 15.85, ymax = 16.0)

c1

# consensus - MICOM
micom.consensus <- micom.consensus %>% mutate_if(is.numeric, ~round(., 6))
d <- ggscatter(micom.consensus, x =  'Replication_rate', y = 'MICOM',
               xlab = '',
               ylab = '',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 'd')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = max(micom.consensus$MICOM), size= 7, cor.coef.name = 'rho')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        axis.text.x = element_blank())

d

# consensus - coco
# calculate the correlation for different parameter values
coco.consensus <- coco.consensus %>% mutate_if(is.numeric, ~round(., 6))
cor_result <- data.frame()

for (i in 4:length(coco.consensus)){
  if (sum(is.na(coco.consensus[, i])) < 70){
    cor_r <- cor.test(coco.consensus[,3],coco.consensus[,i], method = 'spearman',exact = F)
    cor_result[i-3,1] <- cor_r$estimate
  }
  else {
    cor_result[i-3,1] <- NA
  }
}

colnames(cor_result) <- c('cor')
# calculate the mean correlation value and find the one most close to the mean value
closest_value <- cor_result$cor[which.min(abs(cor_result$cor - mean(cor_result$cor, na.rm = TRUE)))]

na <- colnames(coco.consensus)[which(cor_result$cor == closest_value)+3]

e <- ggscatter(coco.consensus, x =  'Replication_rate', y = na,
               xlab = '',
               ylab = '',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 'e')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 80, size= 7, cor.coef.name = 'rho')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        axis.text.x = element_blank())

e

# Distribution of correlation
di_e <- ggplot(cor_result, aes(x = cor)) +
  geom_histogram(binwidth = 0.005, fill = "gray", color = "black") +
  labs(x = "", y = "") +
  font("y.text", size = 13)+
  font("x.text", size = 13)+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        panel.background = element_blank(),
        panel.grid = element_blank())

# Combine the plots
e1 <- e + annotation_custom(ggplotGrob(di_e), xmin = 0.1, xmax = 0.8, ymin = 36, ymax = 80)

e1

# CarveMe - IMIC
carveme <- carveme %>% mutate_if(is.numeric, ~round(., 6))
f <- ggscatter(carveme, x =  'Replication_rate', y = 'IMIC',
               xlab = '',
               ylab = 'Predicted growth rate [1/h]',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 'f')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = max(carveme$IMIC), size= 7, cor.coef.name = 'rho')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        axis.text.x = element_blank())

f
# CarveMe - MICOM-ab1
micom.carveme <- micom.carveme %>% mutate_if(is.numeric, ~round(., 6))
g <- ggscatter(micom.carveme, x =  'Replication_rate', y = 'MICOM_ab1',
               xlab = '',
               ylab = '')+
  labs(title = 'g')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        axis.text.x = element_blank())+
  scale_y_continuous(
    limits = c(36.9, 37.2),
    breaks = seq(36.9, 37.2, by = 0.1),
    labels = number_format(accuracy = 0.1))

g

# CarveMe - coco-ab1
# calculate the correlation for different parameter values
coco.carveme <- coco.carveme %>% mutate_if(is.numeric, ~round(., 6))
cor_result <- data.frame()

for (i in 4:length(coco.carveme.ab1)){
  if (sum(is.na(coco.carveme.ab1[, i])) < 70){
    cor_r <- cor.test(coco.carveme.ab1[,3],coco.carveme.ab1[,i], method = 'spearman',exact = F)
    cor_result[i-3,1] <- cor_r$estimate
  }
  else {
    cor_result[i-3,1] <- NA
  }
}

colnames(cor_result) <- c('cor')
# calculate the mean correlation value and find the one most close to the mean value
closest_value <- cor_result$cor[which.min(abs(cor_result$cor - mean(cor_result$cor, na.rm = TRUE)))]

na <- colnames(coco.carveme.ab1)[which(cor_result$cor == closest_value)+3]

h <- ggscatter(coco.carveme.ab1, x =  'Replication_rate', y = na,
               xlab = '',
               ylab = '',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 'h')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 50, size= 7, cor.coef.name = 'rho')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        axis.text.x = element_blank())

h

# Distribution of correlation
di_h <- ggplot(cor_result, aes(x = cor)) +
  geom_histogram(binwidth = 0.005, fill = "gray", color = "black") +
  labs(x = "", y = "") +
  font("y.text", size = 13)+
  font("x.text", size = 13)+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        panel.background = element_blank(),
        panel.grid = element_blank())

# Combine the plots
h1 <- h + annotation_custom(ggplotGrob(di_h), xmin = 0.1, xmax = 0.8, ymin = 5, ymax = 25)

h1

# CarveMe - MICOM
micom.carveme <- micom.carveme %>% mutate_if(is.numeric, ~round(., 6))
u <- ggscatter(micom.carveme, x =  'Replication_rate', y = 'MICOM',
               xlab = '',
               ylab = '',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 'i')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 80, size= 7, cor.coef.name = 'rho')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        axis.text.x = element_blank())

u

# CarveMe - coco
# calculate the correlation for different parameter values
coco.carveme <- coco.carveme %>% mutate_if(is.numeric, ~round(., 6))
cor_result <- data.frame()

for (i in 4:length(coco.carveme)){
  if (sum(is.na(coco.carveme[, i])) < 70){
    cor_r <- cor.test(coco.carveme[,3],coco.carveme[,i], method = 'spearman',exact = F)
    cor_result[i-3,1] <- cor_r$estimate
  }
  else {
    cor_result[i-3,1] <- NA
  }
}

colnames(cor_result) <- c('cor')
# calculate the mean correlation value and find the one most close to the mean value
closest_value <- cor_result$cor[which.min(abs(cor_result$cor - mean(cor_result$cor, na.rm = TRUE)))]

na <- colnames(coco.carveme)[which(cor_result$cor == closest_value)+3]

j <- ggscatter(coco.carveme, x =  'Replication_rate', y = na,
               xlab = '',
               ylab = '',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 'j')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 120, size= 7, cor.coef.name = 'rho')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        axis.text.x = element_blank())

j

# Distribution of correlation
di_j <- ggplot(cor_result, aes(x = cor)) +
  geom_histogram(binwidth = 0.005, fill = "gray", color = "black") +
  labs(x = "", y = "") +
  font("y.text", size = 13)+
  font("x.text", size = 13)+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        panel.background = element_blank(),
        panel.grid = element_blank())

# Combine the plots
j1 <- j + annotation_custom(ggplotGrob(di_j), xmin = 0.1, xmax = 0.8, ymin = 55, ymax = 120)

j1

# gapseq - IMIC
gapseq <- gapseq %>% mutate_if(is.numeric, ~round(., 6))
k <- ggscatter(gapseq, x =  'Replication_rate', y = 'IMIC',
               xlab = '',
               ylab = 'Predicted growth rate [1/h]',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 'k')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 30, size= 7, cor.coef.name = 'rho')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        axis.text.x = element_blank())

k
# gapseq - MICOM-ab1
micom.gapseq <- micom.gapseq %>% mutate_if(is.numeric, ~round(., 6))
l <- ggscatter(micom.gapseq, x =  'Replication_rate', y = 'MICOM_ab1',
               xlab = '',
               ylab = '')+
  labs(title = 'l')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        axis.text.x = element_blank())+
  scale_y_continuous(
    limits = c(31.9, 32.2),
    breaks = seq(31.9, 32.2, by = 0.1),
    labels = number_format(accuracy = 0.1))

l
# gapseq - coco-ab1
# calculate the correlation for different parameter values
coco.gapseq.ab1 <- coco.gapseq.ab1 %>% mutate_if(is.numeric, ~round(., 5))
cor_result <- data.frame()

for (i in 4:length(coco.gapseq.ab1)){
  if (sum(is.na(coco.gapseq.ab1[, i])) < 70){
    cor_r <- cor.test(coco.gapseq.ab1[,3],coco.gapseq.ab1[,i], method = 'spearman',exact = F)
    cor_result[i-3,1] <- cor_r$estimate
  }
  else {
    cor_result[i-3,1] <- NA
  }
}

colnames(cor_result) <- c('cor')
# calculate the mean correlation value and find the one most close to the mean value
closest_value <- cor_result$cor[which.min(abs(cor_result$cor - mean(cor_result$cor, na.rm = TRUE)))]

na <- colnames(coco.gapseq.ab1)[which(cor_result$cor == closest_value)+3]

m <- ggscatter(coco.gapseq.ab1, x =  'Replication_rate', y = na,
               xlab = '',
               ylab = '',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 'm')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 43, size= 7, cor.coef.name = 'rho')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        axis.text.x = element_blank())

m

# Distribution of correlation
di_m <- ggplot(cor_result, aes(x = cor)) +
  geom_histogram(binwidth = 0.005, fill = "gray", color = "black") +
  labs(x = "", y = "") +
  font("y.text", size = 13)+
  font("x.text", size = 13)+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        panel.background = element_blank(),
        panel.grid = element_blank())

# Combine the plots
m1 <- m + annotation_custom(ggplotGrob(di_m), xmin = 0.1, xmax = 0.8, ymin = 5, ymax = 25)

m1

# gapseq - MICOM
micom.gapseq <- micom.gapseq %>% mutate_if(is.numeric, ~round(., 6))
n <- ggscatter(micom.gapseq, x =  'Replication_rate', y = 'MICOM',
               xlab = '',
               ylab = '',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 'n')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 140, size= 7, cor.coef.name = 'rho')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        axis.text.x = element_blank())

n

# gapseq - coco
# calculate the correlation for different parameter values
coco.gapseq <- coco.gapseq %>% mutate_if(is.numeric, ~round(., 6))
cor_result <- data.frame()

for (i in 4:length(coco.gapseq)){
  if (sum(is.na(coco.gapseq[, i])) < 70){
    cor_r <- cor.test(coco.gapseq[,3],coco.gapseq[,i], method = 'spearman',exact = F)
    cor_result[i-3,1] <- cor_r$estimate
  }
  else {
    cor_result[i-3,1] <- NA
  }
}

colnames(cor_result) <- c('cor')
# calculate the mean correlation value and find the one most close to the mean value
closest_value <- cor_result$cor[which.min(abs(cor_result$cor - mean(cor_result$cor, na.rm = TRUE)))]

na <- colnames(coco.gapseq)[which(cor_result$cor == closest_value)+3]

o <- ggscatter(coco.gapseq, x =  'Replication_rate', y = na[1],
               xlab = '',
               ylab = '',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 'o')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 150, size= 7, cor.coef.name = 'rho')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        axis.text.x = element_blank())

o
# Distribution of correlation
di_o <- ggplot(cor_result, aes(x = cor)) +
  geom_histogram(binwidth = 0.005, fill = "gray", color = "black") +
  labs(x = "",y = "") +
  font("y.text", size = 13)+
  font("x.text", size = 13)+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        panel.background = element_blank(),
        panel.grid = element_blank())

# Combine the plots
o1 <- o + annotation_custom(ggplotGrob(di_o), xmin = 0.1, xmax = 0.8, ymin = 70, ymax = 150)

o1

# kbase - IMIC
kbase <- kbase %>% mutate_if(is.numeric, ~round(., 6))
p <- ggscatter(kbase, x =  'Replication_rate', y = 'IMIC',
               xlab = 'CoPTR-log2(PTR)',
               ylab = 'Predicted growth rate [1/h]',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 'p')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  font("xlab", size = 18)+
  font("x.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 40, size= 7, cor.coef.name = 'rho')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA))

p

# kbase - MICOM-ab1
micom.kbase <- micom.kbase %>% mutate_if(is.numeric, ~round(., 6))
q <- ggscatter(micom.kbase, x =  'Replication_rate', y = 'MICOM_ab1',
               xlab = 'CoPTR-log2(PTR)',
               ylab = '',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 'q')+
  font("title", size = 18, face = 'bold')+
  font("y.text", size = 16)+
  font("xlab", size = 18)+
  font("x.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 47, size= 7, cor.coef.name = 'rho')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA))

q

# kbase - coco-ab1
# calculate the correlation for different parameter values
coco.kbase.ab1 <- coco.kbase.ab1 %>% mutate_if(is.numeric, ~round(., 6))
cor_result <- data.frame()

for (i in 4:length(coco.kbase.ab1)){
  if (sum(is.na(coco.kbase.ab1[, i])) < 70){
    cor_r <- cor.test(coco.kbase.ab1[,3],coco.kbase.ab1[,i], method = 'spearman',exact = F)
    cor_result[i-3,1] <- cor_r$estimate
  }
  else {
    cor_result[i-3,1] <- NA
  }
}

colnames(cor_result) <- c('cor')
# calculate the mean correlation value and find the one most close to the mean value
closest_value <- cor_result$cor[which.min(abs(cor_result$cor - mean(cor_result$cor, na.rm = TRUE)))]

na <- colnames(coco.kbase.ab1)[which(cor_result$cor == closest_value)+3]

r <- ggscatter(coco.kbase.ab1, x =  'Replication_rate', y = na,
               xlab = 'CoPTR-log2(PTR)',
               ylab = '',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 'r')+
  font("title", size = 18, face = 'bold')+
  font("y.text", size = 16)+
  font("xlab", size = 18)+
  font("x.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 50, size= 7, cor.coef.name = 'rho')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA))

r
# Distribution of correlation
di_r <- ggplot(cor_result, aes(x = cor)) +
  geom_histogram(binwidth = 0.005, fill = "gray", color = "black") +
  labs(x = "",y = "") +
  font("y.text", size = 13)+
  font("x.text", size = 13)+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        panel.background = element_blank(),
        panel.grid = element_blank())

# Combine the plots
r1 <- r + annotation_custom(ggplotGrob(di_r), xmin = 0.1, xmax = 0.8, ymin = 10, ymax = 30)

r1

# kbase - MICOM
micom.kbase <- micom.kbase %>% mutate_if(is.numeric, ~round(., 6))
s <- ggscatter(micom.kbase, x =  'Replication_rate', y = 'MICOM',
               xlab = 'CoPTR-log2(PTR)',
               ylab = '',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 's')+
  font("title", size = 18, face = 'bold')+
  font("y.text", size = 16)+
  font("xlab", size = 18)+
  font("x.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 100, size= 7, cor.coef.name = 'rho')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA))

s

# kbase - coco
# calculate the correlation for different parameter values
coco.kbase <- coco.kbase %>% mutate_if(is.numeric, ~round(., 6))
cor_result <- data.frame()

for (i in 4:length(coco.kbase)){
  if (sum(is.na(coco.kbase[, i])) < 70){
    cor_r <- cor.test(coco.kbase[,3],coco.kbase[,i], method = 'spearman',exact = F)
    cor_result[i-3,1] <- cor_r$estimate
  }
  else {
    cor_result[i-3,1] <- NA
  }
}

colnames(cor_result) <- c('cor')
# calculate the mean correlation value and find the one most close to the mean value
closest_value <- cor_result$cor[which.min(abs(cor_result$cor - mean(cor_result$cor, na.rm = TRUE)))]

na <- colnames(coco.kbase)[which(cor_result$cor == closest_value)+3]

t <- ggscatter(coco.kbase, x =  'Replication_rate', y = na[1],
               xlab = 'CoPTR-log2(PTR)',
               ylab = '',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 't')+
  font("title", size = 18, face = 'bold')+
  font("y.text", size = 16)+
  font("xlab", size = 18)+
  font("x.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 180, size= 7, cor.coef.name = 'rho')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA))

t
# Distribution of correlation
di_t <- ggplot(cor_result, aes(x = cor)) +
  geom_histogram(binwidth = 0.005, fill = "gray", color = "black") +
  labs(x = "",y = "") +
  font("y.text", size = 13)+
  font("x.text", size = 13)+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        panel.background = element_blank(),
        panel.grid = element_blank())

# Combine the plots
t1 <- t + annotation_custom(ggplotGrob(di_t), xmin = 0.1, xmax = 0.8, ymin = 70, ymax = 180)

t1

library(patchwork)
a + b + c1 + d + e1 + f + g + h1 + u + j1 + k + l + m1 + n + o1 + p + q + r1 + s + t1 + plot_layout(nrow = 4)

ggsave(filename = '~/IMIC/Figure/Fig 4.svg',
       width = 64,height = 40,units = "cm")
