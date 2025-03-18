topDir <- "~/IMIC/study_case/table/flux_sum_analysis/"
e_data <- read.table(paste0(topDir, 'equal_flux_sum.csv'),
                     header = T, sep = ',')
l_data <- read.table(paste0(topDir, 'less_flux_sum.csv'),
                     header = T, sep = ',')
m_data <- read.table(paste0(topDir, 'more_flux_sum.csv'),
                     header = T, sep = ',')

# Define the threshold
threshold <- 10^-5

# Apply the condition to numeric columns only
numeric_e <- sapply(e_data, is.numeric)
numeric_l <- sapply(l_data, is.numeric) 
numeric_m <- sapply(m_data, is.numeric) 

e_data[, numeric_e] <- lapply(e_data[, numeric_e], function(x) {
  x[x < threshold] <- 0  # Replace values less than the threshold with zero
  return(x)
})

l_data[, numeric_l] <- lapply(l_data[, numeric_l], function(x) {
  x[x < threshold] <- 0  # Replace values less than the threshold with zero
  return(x)
})

m_data[, numeric_m] <- lapply(m_data[, numeric_m], function(x) {
  x[x < threshold] <- 0  # Replace values less than the threshold with zero
  return(x)
})

new_e_data <- e_data[,-1]
new_l_data <- l_data[,-1]
new_m_data <- m_data[,-1]

# calculate the mean value for the replicate sample at each time point
row_means_e<- t(apply(new_e_data, 1, function(row) {
  # Reshape the row to have 3 columns
  row_matrix <- matrix(row, ncol = 3, byrow = TRUE)
  # Calculate the mean of each row in the reshaped matrix
  row_means_e <- apply(row_matrix, 1, mean)
  return(row_means_e)
}))

row_means_l <- t(apply(new_l_data, 1, function(row) {
  # Reshape the row to have 3 columns
  row_matrix <- matrix(row, ncol = 3, byrow = TRUE)
  # Calculate the mean of each row in the reshaped matrix
  row_means_l <- apply(row_matrix, 1, mean)
  return(row_means_l)
}))

row_means_m <- t(apply(new_m_data, 1, function(row) {
  # Reshape the row to have 3 columns
  row_matrix <- matrix(row, ncol = 3, byrow = TRUE)
  # Calculate the mean of each row in the reshaped matrix
  row_means_m <- apply(row_matrix, 1, mean)
  return(row_means_m)
}))

row_means_e <- as.data.frame(row_means_e)
row_means_l <- as.data.frame(row_means_l)
row_means_m <- as.data.frame(row_means_m)

metabolite_id <- data.frame(e_data[,1])

met_e <- data.frame(c(metabolite_id, row_means_e))
met_l <- data.frame(c(metabolite_id, row_means_l))
met_m <- data.frame(c(metabolite_id, row_means_m))

colnames(met_e) <- c('ID','h0','h4','h8','h24') 
colnames(met_l) <- c('ID','h0','h4','h8','h24') 
colnames(met_m) <- c('ID','h0','h4','h8','h24') 

# filter out the metabolites haveing 0 flux sum values
#met_e <- met_e[which(rowSums(met_e[,-1]) != 0),]
#met_l <- met_l[which(rowSums(met_l[,-1]) != 0),]
#met_m <- met_m[which(rowSums(met_m[,-1]) != 0),]

# remove the proton
met_e <- met_e[-which(met_e$ID == 'MNXM1[e]'),]
met_l <- met_l[-which(met_l$ID == 'MNXM1[e]'),]
met_m <- met_m[-which(met_m$ID == 'MNXM1[e]'),]

# find the top 20 import metabolites over all time points
top_20_e <- met_e$ID[order(rowSums(met_e[,-1]), decreasing = T)[1:20]]
top_20_l <- met_l$ID[order(rowSums(met_l[,-1]), decreasing = T)[1:20]]
top_20_m <- met_m$ID[order(rowSums(met_m[,-1]), decreasing = T)[1:20]]

top_mets <- as.data.frame(unique(c(top_20_e,top_20_l,top_20_m)))
colnames(top_mets) <- 'ID'

#write.csv(top_mets, '~/IMIC/study_case/table/top_import_mets.csv')

topDir <- "~/IMIC/study_case/table/"

top_mets <- m_data <- read.table(paste0(topDir, 'top_import_mets.csv'),
                                 header = T, sep = ',')

# heatmap
top_mets_flux <- data.frame()

# 0h
for (i in 1:nrow(top_mets)){
  top_mets_flux[i,1] <- met_e$h0[which(met_e$ID == top_mets[i,1])]
}

num = nrow(top_mets_flux)
for (i in 1:nrow(top_mets)){
  top_mets_flux[i+num,1] <- met_l$h0[which(met_l$ID == top_mets[i,1])]
}

num = nrow(top_mets_flux)
for (i in 1:nrow(top_mets)){
  top_mets_flux[i+num,1] <- met_m$h0[which(met_m$ID == top_mets[i,1])]
}

# 4h
num = nrow(top_mets_flux)
for (i in 1:nrow(top_mets)){
  top_mets_flux[i+num,1] <- met_e$h4[which(met_e$ID == top_mets[i,1])]
}

num = nrow(top_mets_flux)
for (i in 1:nrow(top_mets)){
  top_mets_flux[i+num,1] <- met_l$h4[which(met_l$ID == top_mets[i,1])]
}

num = nrow(top_mets_flux)
for (i in 1:nrow(top_mets)){
  top_mets_flux[i+num,1] <- met_m$h4[which(met_m$ID == top_mets[i,1])]
}

# 8h
num = nrow(top_mets_flux)
for (i in 1:nrow(top_mets)){
  top_mets_flux[i+num,1] <- met_e$h8[which(met_e$ID == top_mets[i,1])]
}

num = nrow(top_mets_flux)
for (i in 1:nrow(top_mets)){
  top_mets_flux[i+num,1] <- met_l$h8[which(met_l$ID == top_mets[i,1])]
}

num = nrow(top_mets_flux)
for (i in 1:nrow(top_mets)){
  top_mets_flux[i+num,1] <- met_m$h8[which(met_m$ID == top_mets[i,1])]
}

# 24h
num = nrow(top_mets_flux)
for (i in 1:nrow(top_mets)){
  top_mets_flux[i+num,1] <- met_e$h24[which(met_e$ID == top_mets[i,1])]
}

num = nrow(top_mets_flux)
for (i in 1:nrow(top_mets)){
  top_mets_flux[i+num,1] <- met_l$h24[which(met_l$ID == top_mets[i,1])]
}

num = nrow(top_mets_flux)
for (i in 1:nrow(top_mets)){
  top_mets_flux[i+num,1] <- met_m$h24[which(met_m$ID == top_mets[i,1])]
}

top_mets_flux$ID <- as.factor(rep(top_mets$Name,12))
colnames(top_mets_flux) <- c('Flux_sum','ID')
top_mets_flux$Ratio <- rep(c(rep('1:1',nrow(top_mets)),rep('1:1000',nrow(top_mets)),rep('1000:1',nrow(top_mets))),4)
top_mets_flux$Time <- as.factor(c(rep(rep('0h',nrow(top_mets)),3), rep(rep('4h',nrow(top_mets)),3), rep(rep('8h',nrow(top_mets)),3),
                                  rep(rep('24h',nrow(top_mets)),3)))

top_mets_flux$Flux_sum <- log10(top_mets_flux$Flux_sum)

library(viridis)
library(ggplot2)
library(ggpubr)
library(forcats)

ggplot(top_mets_flux, aes(Ratio, fct_inorder(ID), fill= Flux_sum)) + 
  geom_tile() +
  scale_fill_viridis(discrete=FALSE) +
  facet_grid(cols = vars(fct_inorder(top_mets_flux$Time)))+
  labs(fill = 'log10(flux sum)')+
  font("x.text", size = 15)+
  font("y.text", size = 15)+
  font('legend.text', size = 12)+
  font('legend.title', size = 15)+
  theme(axis.title = element_blank(),
        strip.text = element_text(size = 15))

ggsave(filename = "~/IMIC/Figure/Fig S10.svg",
       width = 30, height = 25,units = "cm")


