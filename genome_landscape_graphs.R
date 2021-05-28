library(data.table)
library(tidyverse)
library(cowplot)

df <- fread('pcawg_long.txt')
#sampnum <- NROW(unique(df$Sample))
sampchoice <- sample(unique(df$Sample), 8, replace = FALSE)


genes <- fread('genes_hg19.tsv')
df <- left_join(df, genes, by = c('ensembl_gene_id'))
count <- 1
pl <- vector(mode = 'list', length = 8)
dfchoice <- df %>% filter(Sample %in% sampchoice)
for (i in sampchoice) {
  p <- df %>% 
    filter(Sample == i) %>%
    ggplot(aes(x = start_position, y = as.factor(Copy_Number), colour = as.factor(Copy_Number))) + geom_point() + theme(legend.position = 'none') + facet_wrap(~chromosome_name, scales = 'free_x', nrow = 1) + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) + ylab('') + xlab('') + theme(plot.margin = unit(c(0, 0, 0, 0), "cm"))
  ggsave(p, file = paste0('genome_landscape_', i, '.png'))
  pl[[count]] <- p
  count <- count + 1
}

#patch <- ( pl[[1]] | pl[[2]] | pl[[3]] | pl[[4]] ) / ( pl[[5]] | pl[[6]] | pl[[7]] | pl[[8]] ) / ( pl[[9]] | pl[[10]] | pl[[11]] | pl[[12]] ) / ( pl[[13]] | pl[[14]] | pl[[15]] | pl[[16]] )  
sav <- cowplot::plot_grid(plotlist = pl, nrow = 8)
ggsave(sav, file = 'sixteen_genome_landscapes.png', width = 25, height = 14)