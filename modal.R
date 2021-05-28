library(data.table)
library(tidyverse)

all = fread('overlaps_ascat_wgd_CN_cut.bed')
avg <- fread('sample_average_ploidy_ASCAT.txt', header = T, sep = '\t')
genes <- fread('genes_hg38.bed')

modes <- all %>% filter(V3 != -1) %>% group_by(V2) %>% count(V3) %>% arrange(desc(n)) %>% filter(row_number() == 1)

table(modes$V3) %>% as.data.frame() %>% ggplot(aes(x = Var1, y = Freq, fill = Var1)) + geom_col() + xlab('Modal Copy Number') + ylab('Number of WGD Samples') + theme(legend.position = 'none')

ggplot(fif, aes(x = V2.y, y = as.factor(V3.x), color = as.factor(V3.x))) + geom_point() + facet_wrap(~V1.y, scales = 'free_x') + ylab('Copy Number') + ggtitle('Sample 537 (TCGA-DD-A1E9-11A-11D-A151-01)') + labs(color = '') + xlab('Gene start position') + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

twos <- modes %>% filter(V3 == 2) %>% pull(V2) 
all_twos <- all %>% filter(V2 %in% twos)

for (i in 1:10) {
  df <- 
}