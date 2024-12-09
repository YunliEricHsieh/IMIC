topDir <- "~/IMIC/table/"
consensus <- read.table(paste0(topDir, 'predicted_growth/consensus_results_table.csv'),
                   header = T, sep = ',')
coco.consensus <- read.table(paste0(topDir, 'parameter_test/consensus_coco_test_alpha_0.5.csv'),
                        header = T, sep = ',')
coco.consensus.ab1 <- read.table(paste0(topDir, 'parameter_test/consensus_coco_test_with_abundance_1.csv'),
                             header = T, sep = ',')
carveme <- read.table(paste0(topDir, 'predicted_growth/carveme_results_table.csv'),
                        header = T, sep = ',')
coco.carveme <- read.table(paste0(topDir, 'parameter_test/carveme_coco_test_alpha_0.5.csv'),
                             header = T, sep = ',')
coco.carveme.ab1 <- read.table(paste0(topDir, 'parameter_test/carveme_coco_test_with_abundance_1.csv'),
                                 header = T, sep = ',')
gapseq <- read.table(paste0(topDir, 'predicted_growth/gapseq_results_table.csv'),
                        header = T, sep = ',')
coco.gapseq <- read.table(paste0(topDir, 'parameter_test/gapseq_coco_test_alpha_0.5.csv'),
                           header = T, sep = ',')
coco.gapseq.ab1 <- read.table(paste0(topDir, 'parameter_test/gapseq_coco_test_with_abundance_1.csv'),
                               header = T, sep = ',')
kbase <- read.table(paste0(topDir, 'predicted_growth/kbase_results_table.csv'),
                        header = T, sep = ',')
coco.kbase <- read.table(paste0(topDir, 'parameter_test/kbase_coco_test_alpha_0.5.csv'),
                          header = T, sep = ',')
coco.kbase.ab1 <- read.table(paste0(topDir, 'parameter_test/kbase_coco_test_with_abundance_1.csv'),
                              header = T, sep = ',')

library(ggplot2)
library(ggpubr)

# IMIC compare with replication rate
# consensus - IMIC
a <- ggscatter(consensus, x =  'Replication_rate', y = 'IMIC',
                xlab = '',
                ylab = 'Predicted growth rate [1/h]',
                add = "reg.line",                                 # Add regression line
                conf.int = TRUE,                                  # Add confidence interval
                add.params = list(color = "blue",
                                  fill = "lightgray"))+
  labs(title = 'a.')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 40, size= 7, cor.coef.name = 'r')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        axis.text.x = element_blank())

a
# consensus - MICOM-ab1
b <- ggscatter(consensus, x =  'Replication_rate', y = 'MICOM_ab_1',
               xlab = '',
               ylab = '',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 'b.')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 140, size= 7, cor.coef.name = 'r')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        axis.text.x = element_blank())

b
# consensus - coco-ab1
# calculate the correlation for different parameter values
cor_result <- data.frame()

for (i in 4:length(coco.consensus.ab1)){
  if (sum(is.na(coco.consensus.ab1[, i])) < 70){
    cor_r <- cor.test(coco.consensus.ab1[,2],coco.consensus.ab1[,i], method = 'spearman',exact = F)
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
  ylim(53,56)+
  labs(title = 'c.')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 54.7, size= 7, cor.coef.name = 'r')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
    axis.text.x = element_blank())

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
c1 <- c + annotation_custom(ggplotGrob(di), xmin = 0.1, xmax = 0.8, ymin = 54.5, ymax = 56)

c1

# consensus - MICOM
d <- ggscatter(consensus, x =  'Replication_rate', y = 'MICOM',
               xlab = '',
               ylab = '',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 'd.')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 140, size= 7, cor.coef.name = 'r')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        axis.text.x = element_blank())

d

# consensus - coco
# calculate the correlation for different parameter values
cor_result <- data.frame()

for (i in 4:length(coco.consensus)){
  if (sum(is.na(coco.consensus[, i])) < 70){
    cor_r <- cor.test(coco.consensus[,2],coco.consensus[,i], method = 'spearman',exact = F)
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

e <- ggscatter(coco.consensus, x =  'Replication_rate', y = na[1],
               xlab = '',
               ylab = '',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  ylim(-20, 170)+
  labs(title = 'e.')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 100, size= 7, cor.coef.name = 'r')+
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
e1 <- e + annotation_custom(ggplotGrob(di_e), xmin = 0.1, xmax = 0.8, ymin = 100, ymax = 180)

e1

# CarveMe - IMIC
f <- ggscatter(carveme, x =  'Replication_rate', y = 'IMIC',
               xlab = '',
               ylab = 'Predicted growth rate [1/h]',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 'f.')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 40, size= 7, cor.coef.name = 'r')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        axis.text.x = element_blank())

f
# CarveMe - MICOM-ab1
g <- ggscatter(carveme, x =  'Replication_rate', y = 'MICOM_ab_1',
               xlab = '',
               ylab = '',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  ylim(0,170)+
  labs(title = 'g.')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 150, size= 7, cor.coef.name = 'r')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        axis.text.x = element_blank())

g

# CarveMe - coco-ab1
# calculate the correlation for different parameter values
cor_result <- data.frame()

for (i in 4:length(coco.carveme.ab1)){
  if (sum(is.na(coco.carveme.ab1[, i])) < 70){
    cor_r <- cor.test(coco.carveme.ab1[,2],coco.carveme.ab1[,i], method = 'spearman',exact = F)
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
  ylim(20, 170)+
  labs(title = 'h.')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 130, size= 7, cor.coef.name = 'r')+
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
h1 <- h + annotation_custom(ggplotGrob(di_h), xmin = 0.1, xmax = 0.8, ymin = 110, ymax = 175)

h1

# CarveMe - MICOM
u <- ggscatter(carveme, x =  'Replication_rate', y = 'MICOM',
               xlab = '',
               ylab = '',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 'i.')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 80, size= 7, cor.coef.name = 'r')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        axis.text.x = element_blank())

u

# CarveMe - coco
# calculate the correlation for different parameter values
cor_result <- data.frame()

for (i in 4:length(coco.carveme)){
  if (sum(is.na(coco.carveme[, i])) < 70){
    cor_r <- cor.test(coco.carveme[,2],coco.carveme[,i], method = 'spearman',exact = F)
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

j <- ggscatter(coco.carveme, x =  'Replication_rate', y = na[1],
               xlab = '',
               ylab = '',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  ylim(-10, 100)+
  labs(title = 'j.')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 80, size= 7, cor.coef.name = 'r')+
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
j1 <- j + annotation_custom(ggplotGrob(di_j), xmin = 0.1, xmax = 0.8, ymin = 55, ymax = 105)

j1

# gapseq - IMIC
k <- ggscatter(gapseq, x =  'Replication_rate', y = 'IMIC',
               xlab = '',
               ylab = 'Predicted growth rate [1/h]',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 'k.')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 30, size= 7, cor.coef.name = 'r')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        axis.text.x = element_blank())

k
# gapseq - MICOM-ab1
l <- ggscatter(gapseq, x =  'Replication_rate', y = 'MICOM_ab_1',
               xlab = '',
               ylab = '',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 'l.')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 150, size= 7, cor.coef.name = 'r')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        axis.text.x = element_blank())

l
# gapseq - coco-ab1
# calculate the correlation for different parameter values
cor_result <- data.frame()

for (i in 4:length(coco.gapseq.ab1)){
  if (sum(is.na(coco.gapseq.ab1[, i])) < 70){
    cor_r <- cor.test(coco.gapseq.ab1[,2],coco.gapseq.ab1[,i], method = 'spearman',exact = F)
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
  labs(title = 'm.')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 90, size= 7, cor.coef.name = 'r')+
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
m1 <- m + annotation_custom(ggplotGrob(di_m), xmin = 0.1, xmax = 0.8, ymin = 10, ymax = 50)

m1

# gapseq - MICOM
n <- ggscatter(gapseq, x =  'Replication_rate', y = 'MICOM',
               xlab = '',
               ylab = '',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 'n.')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 140, size= 7, cor.coef.name = 'r')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA),
        axis.text.x = element_blank())

n

# gapseq - coco
# calculate the correlation for different parameter values
cor_result <- data.frame()

for (i in 4:length(coco.gapseq)){
  if (sum(is.na(coco.gapseq[, i])) < 70){
    cor_r <- cor.test(coco.gapseq[,2],coco.gapseq[,i], method = 'spearman',exact = F)
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

o <- ggscatter(coco.gapseq, x =  'Replication_rate', y = na,
               xlab = '',
               ylab = '',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  ylim(-10, 170)+
  labs(title = 'o.')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 120, size= 7, cor.coef.name = 'r')+
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
o1 <- o + annotation_custom(ggplotGrob(di_o), xmin = 0.1, xmax = 0.8, ymin = 100, ymax = 180)

o1

# kbase - IMIC
p <- ggscatter(kbase, x =  'Replication_rate', y = 'IMIC',
               xlab = 'CoPTR-log2(PTR)',
               ylab = 'Predicted growth rate [1/h]',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 'p.')+
  font("title", size = 18, face = 'bold')+
  font("ylab", size = 18)+
  font("y.text", size = 16)+
  font("xlab", size = 18)+
  font("x.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 40, size= 7, cor.coef.name = 'r')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA))

p

# kbase - MICOM-ab1
q <- ggscatter(kbase, x =  'Replication_rate', y = 'MICOM_ab_1',
               xlab = 'CoPTR-log2(PTR)',
               ylab = '',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 'q.')+
  font("title", size = 18, face = 'bold')+
  font("y.text", size = 16)+
  font("xlab", size = 18)+
  font("x.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 150, size= 7, cor.coef.name = 'r')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA))

q

# kbase - coco-ab1
# calculate the correlation for different parameter values
cor_result <- data.frame()

for (i in 4:length(coco.kbase.ab1)){
  if (sum(is.na(coco.kbase.ab1[, i])) < 70){
    cor_r <- cor.test(coco.kbase.ab1[,2],coco.kbase.ab1[,i], method = 'spearman',exact = F)
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
  labs(title = 'r.')+
  font("title", size = 18, face = 'bold')+
  font("y.text", size = 16)+
  font("xlab", size = 18)+
  font("x.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 170, size= 7, cor.coef.name = 'r')+
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
r1 <- r + annotation_custom(ggplotGrob(di_r), xmin = 0.1, xmax = 0.8, ymin = -20, ymax = 70)

r1

# kbase - MICOM
s <- ggscatter(kbase, x =  'Replication_rate', y = 'MICOM',
               xlab = 'CoPTR-log2(PTR)',
               ylab = '',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 's.')+
  font("title", size = 18, face = 'bold')+
  font("y.text", size = 16)+
  font("xlab", size = 18)+
  font("x.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 120, size= 7, cor.coef.name = 'r')+
  theme(panel.border = element_rect(color = 'black', size = 1, fill = NA))

s

# kbase - coco
# calculate the correlation for different parameter values
cor_result <- data.frame()

for (i in 4:length(coco.kbase)){
  if (sum(is.na(coco.kbase[, i])) < 70){
    cor_r <- cor.test(coco.kbase[,2],coco.kbase[,i], method = 'spearman',exact = F)
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

t <- ggscatter(coco.kbase, x =  'Replication_rate', y = na,
               xlab = 'CoPTR-log2(PTR)',
               ylab = '',
               add = "reg.line",                                 # Add regression line
               conf.int = TRUE,                                  # Add confidence interval
               add.params = list(color = "blue",
                                 fill = "lightgray"))+
  labs(title = 't.')+
  font("title", size = 18, face = 'bold')+
  font("y.text", size = 16)+
  font("xlab", size = 18)+
  font("x.text", size = 16)+
  stat_cor(method = "spearman", label.x = 0.8, label.y = 170, size= 7, cor.coef.name = 'r')+
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
t1 <- t + annotation_custom(ggplotGrob(di_t), xmin = 0.1, xmax = 0.8, ymin = 120, ymax = 220)

t1

library(patchwork)
a + b + c1 + d + e1 + f + g + h1 + u + j1 + k + l + m1 + n + o1 + p + q + r1 + s + t1 + plot_layout(nrow = 4)

ggsave(filename = '~/IMIC/Figure/Fig S2.svg',
       width = 64,height = 40,units = "cm")
