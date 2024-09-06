topDir <- "~/IMIC/table/key_reaction/"

pathway <- read.table(paste0(topDir, 'The list of pathways from key rxns.csv'),
                      header = T, sep = ',')

library(dplyr)
pathway_counts <- pathway %>%
  count(Metabolic.Pathway) %>%
  arrange(desc(n)) 

library(ggplot2)
ggplot(pathway_counts, aes(x = n, y = reorder(`Metabolic.Pathway`,n)))+
  geom_col() +
  labs(x = "", y = "") +
  theme_linedraw()+
  theme(axis.text.y = element_text(size = 13),
        axis.text.x = element_text(size = 13))

ggsave(filename = "~/IMIC/Fig S3.svg",
       width = 15, height = 15,units = "cm")



