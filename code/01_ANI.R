library(ggplot2)
library(dplyr)



#=========================================================#
#                      Load data                          #
#=========================================================#

data = read.csv2('../data/ANI_plot.csv', sep = ';', row.names = NULL)


#==================================#
#                Plot              #
#==================================#
data$Group = factor(data$Group, levels=c('Rlc','R. anhuiense','leg-etli clade','Rhizobium', 'other genera'))

plot2 <- ggplot(data, aes(x=reorder(Strain, -ANI), y=ANI)) +
  geom_point(aes(color = Group, fill = Group, alpha = 0.8), size = 2) +
  theme_bw() +
  labs(x="Genome rank", y = "ANI (%)") +
  scale_colour_manual(values=c("#ed1174" ,"#762a83", "#548cc9", "#47ac95", "#bcbec0")) + 
  theme(legend.title = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.title=element_text(size = 22),
        legend.text = element_text(size = 20),
        panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())

 
plot2
ggsave('../output/ANI.pdf', plot = plot2, width = 25, height = 20, unit = 'cm', useDingbats = FALSE)

