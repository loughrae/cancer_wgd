args = commandArgs(trailingOnly = TRUE)

library(data.table)
library(tidyverse)
library(ggpubr)
library(qs)

long <- fread(args[1])
sourcename <- args[2]
samp_means <- fread(args[3])



if (sourcename %in% c('ascat', 'absolute')) {
  names(long) <- c('ensembl_gene_id', 'Sample', 'Copy_Number')
}

if (sourcename == 'ascat') {
  samp_means$samplename <- as.character(samp_means$samplename)
}


ohnos <- read.table('~/Downloads/strict_ohnologs_full.csv', sep = '\t', header = T) #not actually comma-separated
oh <- unlist(ohnos)
uniq_oh <- unique(oh)

long$ohno <- ifelse(long$ensembl_gene_id %in% uniq_oh, T, F)

long %>% group_by(ohno, ensembl_gene_id) %>% 
  summarize(total = n(), fours = length(Copy_Number[Copy_Number >= 4]), perc4 = length(Copy_Number[Copy_Number >= 4])/n(), perc3 = length(Copy_Number[Copy_Number >= 3])/n(), exact4 = length(Copy_Number[Copy_Number == 4])/n()) %>%
  fwrite(paste0('genecentric_retention_', sourcename, '.txt'))

retain <- function(hh, number, source) {
  print(number)
  df <-  hh %>% group_by(ohno, Sample) %>% summarize(total = n(), nums = length(Copy_Number[Copy_Number >= number]), perc = length(Copy_Number[Copy_Number >= number])/n() , exact = length(Copy_Number[Copy_Number == number])/n()) %>% mutate(percent = nums/total) %>%
    mutate(Ohnolog = ifelse(ohno == T, 'Ohnolog', 'Non-Ohnolog')) %>%
    left_join(samp_means, by = c('Sample' = 'samplename')) 
  p1 <- ggplot(df, aes(x = ploidy, y = percent, colour = Ohnolog)) + geom_point() + facet_wrap(~Ohnolog) + scale_x_reverse() + xlab('Sample mean copy number') + ylab(paste0('Proportion of genes in at least ', number, ' copies')) + theme(legend.position = 'none')
  ggsave(p1, file = paste0('decreasing_ploidy_atleast_', number, '_', source, '.png'), width = 20, height = 15)
  ggplot(df, aes(x = ploidy, y = exact, colour = Ohnolog)) + geom_point() + facet_wrap(~Ohnolog) + scale_x_reverse() + xlab('Sample mean copy number') + ylab(paste0('Proportion of genes in exactly ', number, ' copies')) + theme(legend.position = 'none')
  ggsave(paste0('decreasing_ploidy_exactly_', number, '_', source, '.png'), width = 20, height = 15)
  p2 <- ggplot(df, aes(x = as.factor(Ohnolog), y = percent, fill = as.factor(Ohnolog))) + geom_boxplot() + stat_compare_means() + xlab('Gene Type') + ylab(paste0('Proportion of genes in at least ', number, ' copies')) + theme(legend.position = 'none')
  ggsave(p2, file = paste0('boxplot_atleast_', number, '_', source, '.png'))
 qsave(p1, file = paste0('decreasing_ploidy_atleast_', number, '_', source, '.qs'))
 qsave(p2, file = paste0('boxplot_atleast_', number, '_', source, '.qs'))
  d2 <- hh %>% group_by(ohno, ensembl_gene_id) %>% summarize(total = n(), nums = length(Copy_Number[Copy_Number >= number]), perc = length(Copy_Number[Copy_Number >= number])/n() , exact = length(Copy_Number[Copy_Number == number])/n()) %>% mutate(percent = nums/total) %>%
    mutate(Ohnolog = ifelse(ohno == T, 'Ohnolog', 'Non-Ohnolog'))
  p3 <- ggplot(d2, aes(x = Ohnolog, y = perc, fill = Ohnolog)) + geom_boxplot() + stat_compare_means() + theme(legend.position = 'none')
  ggsave(p3, file = paste0('boxplot_genecentred_', number, '_', source, '.png'))
  return(df)
}

retained_3 <- retain(long, 3, sourcename)
retained_4 <- retain(long, 4, sourcename)

