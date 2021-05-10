library(tidyverse)
library(data.table)

mop <- fread('mapped_hg38_ascat_withinfo_pasted.bed')
colnames(mop) <- c('Chromosome', 'Start', 'End', 'ensembl_gene_id', 'Mean', 'Median', 'Min', 'Max', 'Count', 'SD', 'Info')
none <- mop %>% filter(Count == 0) %>% pull(ensembl_gene_id)
table(mop$Count == 0)
table(mop$Mean == -999)

# 
# mop %>% 
#     filter(Count > 0) %>%
#     slice(1:1000) %>% 
#     group_by(ensembl_gene_id) %>% 
#     separate_rows(Info, sep = ',') %>% 
#   data.table %>%
#   .[,c('CN', 'Cancer', 'WGD') := tstrsplit(Info, "_", fixed = TRUE)] %>%
#                filter(CN == 4) %>% 
#                count() %>% 
#                ggplot(aes(x = n)) + geom_histogram() 
#   
#   st = Sys.time()
#   tro = setDT(trii)[,paste0("Info", 1:3) := tstrsplit(Info, "_", fixed = TRUE)] 
#   en = Sys.time()
#   
#   st2 = Sys.time()
#   tro = separate(trii, col = Info, into = c('Info1', 'Info2', 'Info3'), sep = '_')
#   en2 = Sys.time()
#   
#   
#   st3 = Sys.time()
#   tri1 = mop %>% filter(Count > 0) %>% slice(1:1000) %>% separate_rows(Info, sep = ',') 
#   tri2 = setDT(tri1)[,c('CN', 'Cancer', 'WGD') := tstrsplit(Info, "_", fixed = TRUE)]
#   tri3 = tri2 %>% count(ensembl_gene_id, WGD, CN)
#   en3 = Sys.time()
#   
  mop <- mop %>% filter(Count > 0)
  splitmop <- split(mop , cut(1:nrow(mop), 60))
  process <- function(df) {
    df %>% 
      separate_rows(Info, sep = ',') %>%
      data.table %>%
      .[, c('CN', 'Cancer', 'WGD') := tstrsplit(Info, '_', fixed = TRUE)] %>%
      mutate(copy = case_when(CN > 4 ~ '>4', CN == 4 ~ '4', CN == 3 ~ '3', CN == 2 ~ '2', CN < 2 ~ '<2')) %>%
      count(ensembl_gene_id, WGD, copy)
  }
  
  processed <- lapply(splitmop, process)
  counts <- bind_rows(processed)