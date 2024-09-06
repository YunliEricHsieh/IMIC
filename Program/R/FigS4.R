topDir <- "~/IMIC/table/key_reaction/"

data <- read.table(paste0(topDir, 'the impact of core reaction on community growth rate.csv'),
                   header = T, sep = ',')

library(tidyr)

# Reshape the data to long format
data_long <- pivot_longer(data, cols = everything(), names_to = "Reaction", values_to = "Value")

library(ggplot2)

ggplot(data_long, aes(x= Value, y = Reaction))+
  geom_boxplot()+
  geom_point(position='jitter', alpha = 1)+
  labs(y = '', x = 'Influence ratio')+
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.title.x = element_text(size = 15))

ggsave(filename = "~/IMIC/Figure/Fig S4.svg",
       width = 15, height = 15,units = "cm")
