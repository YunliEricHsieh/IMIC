library(dplyr)
library(ggplot2)
library(forcats)
library(ggpubr)

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

met_number <- data.frame()
met_number <- data.frame(c(sum(met_e$h0 > 0),sum(met_e$h4 > 0),sum(met_e$h8 > 0),sum(met_e$h24 > 0),
                           sum(met_l$h0 > 0),sum(met_l$h4 > 0),sum(met_l$h8 > 0),sum(met_l$h24 > 0),
                           sum(met_m$h0 > 0),sum(met_m$h4 > 0),sum(met_m$h8 > 0),sum(met_m$h24 > 0)))

colnames(met_number) <- 'Number'

met_number$Ratio <- c(rep('1:1',4),rep('1:1000',4),rep('1000:1',4))
met_number$Timepoint <- as.factor(rep(c('0h','4h','8h','24h'),3))

ggplot(met_number, aes(x = fct_inorder(Timepoint), y = Number, fill = Ratio))+
  geom_bar(stat="identity",position='dodge')+
  geom_text(aes(label=Number), vjust=-0.3, size=5, position = position_dodge(0.9))+
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme_linedraw()+
  labs(y = 'Number of imported metabolite')+
  font("x.text", size = 15)+
  font("y.text", size = 15)+
  font("y.title", size = 17)+
  font("legend.text", size = 12)+
  theme(axis.title.x = element_blank(),
        legend.title = element_blank(),
        panel.grid = element_blank())

ggsave(filename = "~/IMIC/Figure/Fig S6.svg",
       width = 18, height = 13,units = "cm")
