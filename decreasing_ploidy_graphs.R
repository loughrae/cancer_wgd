library(data.table)
library(tidyverse)
library(ggpubr)

absect <- fread('overlaps_absolute_wgd_CN_cut.bed')
absamples_filtered <- fread('absolute_samples_WGD_filtered.txt')
ohnos <- read.table('~/Downloads/strict_ohnologs_full.csv', sep = '\t', header = T) #not actually comma-separated
oh <- unlist(ohnos)
uniq_oh <- unique(oh)

absect$ohno <- ifelse(absect$V1 %in% uniq_oh, T, F)

retain <- function(hh, number, source) {
  print(number)
  df <-  hh %>% group_by(ohno, V2) %>% summarize(total = n(), nums = length(V3[V3 >= number]), perc = length(V3[V3 >= number])/n() , exact = length(V3[V3 == number])/n()) %>% mutate(percent = nums/total) %>%
    mutate(Ohnolog = ifelse(ohno == T, 'Ohnolog', 'Non-Ohnolog')) %>%
    left_join(absamples_filtered, by = c('V2' = 'array')) 
  ggplot(df, aes(x = ploidy, y = percent, colour = Ohnolog)) + geom_point() + facet_wrap(~Ohnolog) + scale_x_reverse() + xlab('Sample mean copy number') + ylab(paste0('Proportion of genes in at least ', number, ' copies')) + theme(legend.position = 'none')
  ggsave(paste0('decreasing_ploidy_atleast_', number, '_', source, '.png'), width = 20, height = 15)
  ggplot(df, aes(x = ploidy, y = exact, colour = Ohnolog)) + geom_point() + facet_wrap(~Ohnolog) + scale_x_reverse() + xlab('Sample mean copy number') + ylab(paste0('Proportion of genes in exactly ', number, ' copies')) + theme(legend.position = 'none')
  ggsave(paste0('decreasing_ploidy_exactly_', number, '_', source, '.png'), width = 20, height = 15)
  ggplot(df, aes(x = as.factor(Ohnolog), y = percent)) + geom_boxplot() + stat_compare_means() + xlab('Gene Type') + ylab(paste0('Proportion of genes in at least ', number, ' copies'))
  ggsave(paste0('boxplot_atleast_', number, '_', source, '.png'))
  return(df)
}

retained_3 <- retain(absect, 3, 'absolute')
retained_4 <- retain(absect, 4, 'absolute')
