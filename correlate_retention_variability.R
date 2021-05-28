library(data.table)
library(tidyverse)
library(patchwork)
library(ggpubr)

pl <- vector(mode = 'list', length = 3)
i = 1
for (source in c('ascat', 'absolute', 'pcawg')) {


  
wgd <- fread(paste0('genecentric_retention_', source, '.txt'))
dip <- fread(paste0('mapped_', source, '_all.bed'))
names(dip) <- c('chromosome_name', 'start_position', 'end_position', 'ensembl_gene_id', 'mean', 'median', 'stdev', 'min', 'max', 'count')

tog <- left_join(dip, wgd, by = c('ensembl_gene_id')) %>%
  mutate(Ohnolog = ifelse(ohno == T, 'Ohnolog', 'Non-Ohnolog'))

p <- tog %>%
  filter(stdev != -999) %>%
  ggplot(aes(x = stdev, y = perc4, colour = Ohnolog, group = 1)) + geom_point() + stat_cor() + ggtitle(toupper(source)) + ylab('Proportion of WGD Samples Where CN >= 4') + xlab('SD of Gene CN in Diploid Samples') + labs(colour = 'Human \nOhnolog')

ggsave(p, file = paste0('correlation_retention_variability_', source, '.png'))

pl[[i]] <- p
i = i + 1
}

p1 <- pl[[1]] | pl[[2]] | pl[[3]]
ggsave(p1, file = 'combine_correlation_retention_variability.png', width = 20, height = 12)


