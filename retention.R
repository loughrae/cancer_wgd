library(data.table)
library(tidyverse)

hh = fread('overlaps_ascat_wgd_CN_cut.bed')

sam3 <- hh %>% filter(V3 >= 3) %>% count(V2)
sam4 <- hh %>% filter(V3 >= 4) %>% count(V2)

gene3 <- hh %>% filter(V3 >= 3) %>% count(V1)
gene4 <- hh %>% filter(V3 >= 4) %>% count(V1)


count_threshold <- function(x, kind, cutoff) {
if (kind == 'Sample') {
  total <- 57000
} else {
  total <- 3800
}
  i = 1
  lost <- vector(mode = 'list')
  for (thresh in seq(1, total, by = total/50)) {
    lost[[i]] <- data.frame(thresh, table(x$n >= thresh))
    i = i + 1
    
  }
  all = bind_rows(lost) %>% complete(thresh, Var1)
 ggplot(all, aes(x = thresh, y = Freq, fill = Var1)) + geom_col(position = 'dodge') + xlab('Threshold') + ylab('Frequency') + labs(fill = 'Over Threshold', subtitle = cutoff) + ggtitle(paste0('Number of ', kind, 's meeting each threshold')) +
 ggsave(paste0('threshold_', kind, '_', cutoff, '.png'), width = 20, height = 12)

 
}

count_threshold(sam3, kind = 'Sample', cutoff = 3)
count_threshold(sam4, kind = 'Sample', cutoff = 4)
count_threshold(gene3, kind = 'Gene', cutoff = 3)
count_threshold(gene4, kind = 'Gene', cutoff = 4)

#wazzup big cat
ohnos <- read.table('~/Downloads/strict_ohnologs_full.csv', sep = '\t', header = T) #not actually comma-separated
oh <- unlist(ohnos)
uniq_oh <- unique(oh)

hh$ohno <- ifelse(hh$V1 %in% uniq_oh, T, F)


retain <- function(hh, number) {
  print(number)
  df <-  hh %>% group_by(ohno, V2) %>% summarize(total = n(), nums = length(V3[V3 >= number]), perc = length(V3[V3 >= number])/n() , exact = length(V3[V3 == number])/n()) %>% mutate(percent = nums/total) %>%
    mutate(V2 = as.numeric(V2)) %>%
    mutate(Ohnolog = ifelse(ohno == T, 'Ohnolog', 'Non-Ohnolog')) %>%
    left_join(avg, by = c('V2' = 'index')) 
  ggplot(df, aes(x = mean_CN, y = percent, colour = Ohnolog)) + geom_point() + facet_wrap(~Ohnolog) + scale_x_reverse() + xlab('Sample mean copy number') + ylab(paste0('Proportion of genes in at least ', number, ' copies')) + theme(legend.position = 'none')
  ggsave(paste0('decreasing_ploidy_atleast_', number, '.png'), width = 20, height = 15)
  ggplot(df, aes(x = mean_CN, y = exact, colour = Ohnolog)) + geom_point() + facet_wrap(~Ohnolog) + scale_x_reverse() + xlab('Sample mean copy number') + ylab(paste0('Proportion of genes in exactly ', number, ' copies')) + theme(legend.position = 'none')
  ggsave(paste0('decreasing_ploidy_exactly_', number, '.png'), width = 20, height = 15)
  
}

retain(hh, 3)
retain(hh, 4)

