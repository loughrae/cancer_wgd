library(data.table)
library(tidyverse)
library(patchwork)

pl <- vector(mode = 'list', length = 3)
p2 <- vector(mode = 'list', length = 3)
i = 1
for (sourcename in c('ascat', 'absolute', 'pcawg')) {
  
  if (sourcename == 'ascat') {
    GO <- fread('GOterms_hg38.txt')
  }
  else {
    GO <- fread('GOterms_hg19.txt')
  }
  GO <- GO %>% separate_rows(ANNOTATED_GENES, sep = ',') %>% mutate(ANNOTATED_GENES = trimws(ANNOTATED_GENES))
  wgd <- fread(paste0('genecentric_retention_', sourcename, '.txt'))
  wgd %>% arrange(desc(perc3)) %>% select(ensembl_gene_id) %>% write.table(file = paste0('ranked_genecentric_retention_', sourcename, '.txt'), sep = '\t', col.names = F, row.names = F, quote = F)
  wgd <- wgd %>% left_join(GO, by = c('ensembl_gene_id' = 'ANNOTATED_GENES')) 
  p <- wgd %>% ggplot(aes(x = reorder(TERM, perc3, FUN = median), y = perc3, fill = TERM)) + geom_boxplot() + ylab('Proportion of Samples Where Gene CN >= 3') + xlab('GO Biological Process') + coord_flip() + theme(legend.position = 'none') + ggtitle(sourcename) 
  ggsave(p, file = paste0('GOmap_', sourcename, '.png'))
  pl[[i]] <- p
  i = i + 1
}

p1 <- pl[[1]] | pl[[2]] | pl[[3]]
ggsave(p1, file = 'combined_GOmap.png', width = 20, height = 12)



