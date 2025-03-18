topDir <- "~/IMIC/table/flux_sum_analysis/"
IMIC <- read.table(paste0(topDir, '14_MAG_flux_sum.csv'),
                   header = T, sep = ',')
MICOM <- read.table(paste0(topDir, 'MICOM_flux_sum.csv'),
                    header = T, sep = ',')

# Define the threshold
threshold <- 10^-5

# Apply the condition to numeric columns only
numeric_I <- sapply(IMIC, is.numeric)
numeric_M <- sapply(MICOM, is.numeric) 

# Replace values less than the threshold with zero
IMIC[, numeric_I] <- lapply(IMIC[, numeric_I], function(x) {
  x[x < threshold] <- 0  
  return(x)
})

MICOM[, numeric_M] <- lapply(MICOM[, numeric_M], function(x) {
  x[x < threshold] <- 0 
  return(x)
})


colnames(IMIC) <- c('ID','20d','40d','60d','90d','180d') 
colnames(MICOM) <- c('ID','20d','40d','60d','90d','180d') 

met_number <- data.frame()
met_number <- data.frame(c(sum(IMIC$'20d' > 0,na.rm = T),sum(IMIC$'40d' > 0,na.rm = T),
                           sum(IMIC$'60d' > 0,na.rm = T),sum(IMIC$'90d' > 0,na.rm = T),
                           sum(IMIC$'180d' > 0,na.rm = T),sum(MICOM$'20d' > 0,na.rm = T),
                           sum(MICOM$'40d' > 0,na.rm = T),sum(MICOM$'60d' > 0,na.rm = T),
                           sum(MICOM$'90d' > 0,na.rm = T),sum(MICOM$'180d' > 0,na.rm = T)))

colnames(met_number) <- 'Number'

met_number$Methods <- c(rep('IMIC',5),rep('MICOM',5))
met_number$Timepoint <- as.factor(rep(c('20d','40d','60d','90d','180d'),2))

library(ggplot2)
library(ggpubr)
library(forcats)

ggplot(met_number, aes(x = fct_inorder(Timepoint), y = Number, fill = Methods))+
  geom_bar(stat="identity",position='dodge')+
  geom_text(aes(label=Number), vjust=-0.3, size=4, position = position_dodge(0.9))+
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  theme_linedraw()+
  labs(y = 'Imported metabolites')+
  font("x.text", size = 12)+
  font("y.text", size = 12)+
  font("y.title", size = 14)+
  font("legend.text", size = 12)+
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'top')

ggsave(filename = "~/IMIC/Figure/Fig 6-1.svg",
       width = 12, height = 10,units = "cm")

# define the different type of imported metaoblites
# find the time independent imported metabolite
imic_time_indep <- as.data.frame(IMIC$ID[which(IMIC$'20d'>0 & IMIC$'40d'>0 & IMIC$'60d'>0 &IMIC$'90d'>0 & IMIC$'180d'>0)])
colnames(imic_time_indep) = 'ID'
micom_time_indep <- as.data.frame(MICOM$ID[which(MICOM$'20d'>0 & MICOM$'40d'>0 & MICOM$'60d'>0 &MICOM$'90d'>0 & MICOM$'180d'>0)])
colnames(micom_time_indep) = 'ID'

# find all of the imported metabolite
imic_import <- as.data.frame(IMIC$ID[rowSums(IMIC[, -1], na.rm = T) != 0])
colnames(imic_import) = 'ID'
micom_import <- as.data.frame(MICOM$ID[rowSums(MICOM[, -1], na.rm = T) != 0])
colnames(micom_import) = 'ID'

# find the time dependent imported metabolite
imic_time_dep <- as.data.frame(imic_import$ID[!(imic_import$ID %in% imic_time_indep$ID)])
colnames(imic_time_dep) = 'ID'
micom_time_dep <- as.data.frame(micom_import$ID[!(micom_import$ID %in% micom_time_indep$ID)])
colnames(micom_time_dep) = 'ID'

# create a table for upset plot
met_table = as.data.frame(IMIC$ID)
colnames(met_table) = 'ID'
#met_table$'IMIC all imported mets' <- as.integer(as.logical(IMIC$ID %in% imic_import$ID))
met_table$'IMIC time independent' <- as.integer(as.logical(IMIC$ID %in% imic_time_indep$ID))
met_table$'IMIC time dependent' <- as.integer(as.logical(IMIC$ID %in% imic_time_dep$ID))
#met_table$'MICOM all imported mets' <- as.integer(as.logical(IMIC$ID %in% micom_import$ID))
met_table$'MICOM time independent' <- as.integer(as.logical(IMIC$ID %in% micom_time_indep$ID))
met_table$'MICOM time dependent' <- as.integer(as.logical(IMIC$ID %in% micom_time_dep$ID))

# show the metabolites which are time dependent in IMIC but time independent in MICOM
met_table$ID[which(met_table$`IMIC time dependent` == 1 & met_table$`MICOM time independent` == 1)]
# show the metabolites which are time independent in IMIC but time dependent in MICOM
met_table$ID[which(met_table$`IMIC time independent` == 1 & met_table$`MICOM time dependent` == 1)]

library(UpSetR)
library(svglite)

svglite("~/IMIC/Figure/Fig 6-2.svg", width = 7, height = 5)

upset(met_table, nsets = 6, point.size = 3.5, line.size = 2, 
      mainbar.y.label = "Intersection size", sets.x.label = "Metabolites", 
      order.by = "freq", text.scale = c(2,2,2,2,2,2)) 

dev.off()
