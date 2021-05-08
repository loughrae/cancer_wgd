map <- fread('mapped_hg38_ascat_withinfo_pasted.bed')
colnames(map) <- c('Chromosome', 'Start', 'End', 'ensembl_gene_id', 'Mean', 'Median', 'Min', 'Max', 'Count', 'SD', 'Info')
none <- map %>% filter(Count == 0) %>% pull(ensembl_gene_id)
test <- map %>% filter(Count > 0) %>% slice(1:10) 
table(map$Count == 0)
table(map$Mean == -999)


map %>% 
    filter(Count > 0) %>%
    slice(1:10000) %>% 
    group_by(ensembl_gene_id) %>% 
    separate_rows(col = Info, sep = ',') %>% 
    separate(col = Info, sep = '_', into = c('CN', 'Cancer', 'WGD')) %>% 
               filter(CN == 4) %>% 
               count() %>% 
               ggplot(aes(x = n)) + geom_histogram() 
  
  