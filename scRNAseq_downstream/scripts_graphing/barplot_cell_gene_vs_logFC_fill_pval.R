library(readxl)
library(tidyr)
library(ggplot2)
library(viridis)
library(viridisLite)
library(RColorBrewer)
library(hexbin)
library(reshape2)

# this script plots cell_gene vs logFC, filled with pvalue

df <- as.data.frame(read_excel('UC_compiled.xlsx', sheet = 'for_graph', col_names=TRUE))

df$cell_gene <- paste(df$Case.CellStates, df$gene, sep = " → ")

# remove rows with na
df <- df[complete.cases(df$average_logFC), ]

ggplot(data=df, aes(x=cell_gene, y=average_logFC, fill = average_FDR)) +
  geom_bar(stat="identity") + 
  coord_flip() +
  theme_classic() +
  #remove x and y axis labels
  labs(x="Cell-type DEG\n", y="logFC")+
  theme(
    # make the x axis labels rotate, hjust=0.5, vjust = 0.5 as necc
    axis.text.x = element_text(size = 11, color = 'black', face = 'bold'),
    axis.text.y = element_text(size = 11, color = 'black', face = 'bold'),
    #set thickness of axis ticks
    axis.ticks=element_line(size=1)
  ) +
  
  # legend stuff
  guides(fill = guide_colourbar(barwidth = 1, barheight = 3, 
                                title = 'p value'
                                #title.position = 'top', 
                                #title.hjust = 0.5, title.vjust=0.5)
  )) +
  theme(#legend.position='top', 
    #legend.direction = 'horizontal', 
    legend.title = element_text(color = "black", size = 9, face='bold'), 
    legend.text=element_text(size = 8, face = 'bold'),
    axis.title=element_text(size=11,face="bold"),
    axis.ticks.length = unit(.25, 'cm'),
    # change size of the graph
    aspect.ratio = 1
  ) +
  scale_fill_gradient2(high = "#F62D2D", mid = '#A2264B', low = "#412F88")

