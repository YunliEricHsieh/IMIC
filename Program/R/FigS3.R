topDir <- "~/IMIC/table/variability_analysis/"
cons <- read.table(paste0(topDir, 'consensus_variability_analysis.csv'),
                   header = T, sep = ',')

cons$ratio <- abs(cons$Max_growth_rate-cons$Min_growth_rate)/cons$Max_growth_rate

for (i in 1:nrow(cons)){
  id <- strsplit(cons$ModelID,'_')
  cons$ModelID[i] <- id[[i]][2]
}

cons$ModelID <- factor(cons$ModelID, levels = c('1','2','3','4','5','6','7','8',
                                                '9','10','11','12','13','14'))

consensus <- ggplot(cons, aes(x= ModelID, y = ratio))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position='jitter', alpha = 1)+
  #  ylim(0,1)+
  labs(title = 'a. consensus')+
  xlab('')+
  theme(axis.text = element_blank())+
  ylab('Relative variability\nof growth rate')+
  font("ylab", size = 14)+
  font("title", size = 14, face = 'bold')+
  font("y.text", size = 13)

consensus

carv <- read.table(paste0(topDir, 'carveme_variability_analysis.csv'),
                   header = T, sep = ',')

carv$ratio <- abs(carv$Max_growth_rate-carv$Min_growth_rate)/abs(carv$Max_growth_rate)

for (i in 1:nrow(carv)){
  id <- strsplit(carv$ModelID,'_')
  carv$ModelID[i] <- id[[i]][2]
}

carv$ModelID <- factor(carv$ModelID, levels = c('1','2','3','4','5','6','7','8',
                                                '9','10','11','12','13','14'))

carveme <- ggplot(carv, aes(x= ModelID, y = ratio))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position='jitter', alpha = 1)+
  #  ylim(0,1)+
  labs(title = 'b. CarveMe')+
  xlab('')+
  ylab('')+
  theme(axis.text = element_blank())+
  font("title", size = 14, face = 'bold')
carveme

gaps <- read.table(paste0(topDir, 'gapseq_variability_analysis.csv'),
                   header = T, sep = ',')

gaps$ratio <- abs(gaps$Max_growth_rate-gaps$Min_growth_rate)/abs(gaps$Max_growth_rate)

for (i in 1:nrow(gaps)){
  id <- strsplit(gaps$ModelID,'_')
  gaps$ModelID[i] <- id[[i]][2]
}

gaps$ModelID <- factor(gaps$ModelID, levels = c('1','2','3','4','5','6','7','8',
                                                '9','10','11','12','13','14'))

gapseq <- ggplot(gaps, aes(x= ModelID, y = ratio))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position='jitter', alpha = 1)+
  #  ylim(0,1)+
  labs(title = 'c. gapseq')+
  xlab('MAG ID')+
  ylab('Relative variability\nof growth rate')+
  font("xlab", size = 14)+
  font("ylab", size = 14)+
  font("title", size = 14, face = 'bold')+
  font("x.text", size = 13)+
  font("y.text", size = 13)

gapseq

kbas <- read.table(paste0(topDir, 'kbase_variability_analysis.csv'),
                   header = T, sep = ',')

kbas$ratio <- abs(kbas$Max_growth_rate-kbas$Min_growth_rate)/abs(kbas$Max_growth_rate)

for (i in 1:nrow(kbas)){
  id <- strsplit(kbas$ModelID,'_')
  kbas$ModelID[i] <- id[[i]][2]
}

kbas$ModelID <- factor(kbas$ModelID, levels = c('1','2','3','4','5','6','7','8',
                                                '9','10','11','12','13','14'))

kbase <- ggplot(kbas, aes(x= ModelID, y = ratio))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(position='jitter', alpha = 1)+
  #  ylim(0,1)+
  labs(title = 'd. KBase')+
  xlab('MAG ID')+
  ylab('')+
  theme(axis.text = element_blank())+
  font("xlab", size = 14)+
  font("title", size = 14, face = 'bold')+
  font("x.text", size = 13)

library(patchwork)
consensus+carveme+gapseq+kbase+ plot_layout(nrow = 2)

ggsave(filename = "~/IMIC/Figure/Fig S3.svg",
       width = 20,height = 15,units = "cm")
