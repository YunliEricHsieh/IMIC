topDir <- "~/IMIC/table/key_reaction/"

data <- read.table(paste0(topDir, 'flux_value_of_key_rxns.csv'),
                   header = T, sep = ',')

rxn.list <- unique(data$RxnID)

new.data <- data.frame()
x <- 0

for (i in 1:length(rxn.list)){
  
  tmp.data <- data[which(data$RxnID == rxn.list[i]),]
  y <- 2
  z <- 7
  h <- 12
  
  for (j in 1:5){
    
    ab_cor_test <- cor.test(tmp.data[,y], tmp.data[,z], method = "spearman",exact=FALSE)
    exp_cor_test <- cor.test(tmp.data[,y], tmp.data[,h], method = "spearman",exact=FALSE)
    
    # reaction ID
    new.data[x+j,1] <- rxn.list[i]
    
    # correlation coeficient and p-value
    new.data[x+j,2] <- ab_cor_test$estimate
    new.data[x+j,3] <- ab_cor_test$p.value
    
    # correlation coeficient and p-value
    new.data[x+j,4] <- exp_cor_test$estimate
    new.data[x+j,5] <- exp_cor_test$p.value
    
    y <- y+1
    z <- z+1
    h <- h+1
    
  }
  x = x+5
}

# correct p-value for multiple comparisons
new.data$V3 <- p.adjust(new.data$V3, method="BH")
new.data$V5 <- p.adjust(new.data$V5, method="BH")

# change the column name
colnames(new.data) <- c('RxnID', 'ab_cor', 'ab_adj_p', 'exp_cor', 'exp_adj_p')

# add time point information
new.data$time <- as.factor(rep(c('20d', '40d', '60d', '90d', '180d'), 11))

library(ggplot2)
library(ggpubr)
library(forcats)

a <- ggplot(new.data, aes(x = RxnID, y = ab_cor)) + 
  geom_boxplot(aes(group = RxnID, fill = fct_inorder(time)), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color = fct_inorder(time)), position = position_jitter(width = 0.2), alpha = 0.6, size = 3) +
  labs(x = '', y = 'Spearman correlation coefficient') +
  scale_fill_brewer(palette = "Set1", name = "Time") +  # Ensure 'fill' uses time
  scale_color_brewer(palette = "Set1", name = "Time") +  # Ensure 'color' uses time
  theme_linedraw()+
  ylim(0,1)+
  labs(title = 'a.')+
  font("title", size = 15, face = 'bold')+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 13),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        legend.position = "none")
a

b <- ggplot(new.data, aes(x = RxnID, y = exp_cor)) + 
  geom_boxplot(aes(group = RxnID, fill = fct_inorder(time)), alpha = 0.5, outlier.shape = NA) +
  geom_point(aes(color = fct_inorder(time)), position = position_jitter(width = 0.2), alpha = 0.6, size = 3) +
  labs(x = '', y = '') +
  scale_fill_brewer(palette = "Set1", name = "Time") +  # Ensure 'fill' uses time
  scale_color_brewer(palette = "Set1", name = "Time") +  # Ensure 'color' uses time
  theme_linedraw()+
  ylim(0,1)+
  labs(title = 'b.')+
  font("title", size = 15, face = 'bold')+
  font("legend.text", size = 13)+
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 13),
        axis.text.y = element_blank(),
        axis.title.y = element_text(size = 15),
        legend.position = "right",
        legend.title = element_blank())
b

library(patchwork)
a+b+plot_layout(nrow = 1)

ggsave(filename = "~/IMIC/Figure/Fig 5.svg",
       width = 25, height = 15,units = "cm")
