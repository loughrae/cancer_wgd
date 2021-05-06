library(TCGAbiolinks)
library(tidyverse)

#### Download allele-specific copy number data #### 

allproj <- getGDCprojects()
projs <- allproj[startsWith(allproj$id, 'TCGA'),]$id
dall <- GDCquery(project = projs, data.category = 'Copy Number Variation', workflow.type = 'ASCAT2', data.type = 'Allele-specific Copy Number Segment')
GDCdownload(dall)


#### Case metadata ####

meta <- dall$results[[1]][, c('cases', 'project', 'file_name', 'id')] %>%
  separate(col = 1, into = c('TCGA', 'TSS', 'Individual', 'Sample', 'Portion', 'Plate', 'Center'), sep = '-', remove = FALSE) %>%
  mutate(Patient = paste(TCGA, TSS, Individual, sep = '-')) %>%
  separate(project, into = c('drop', 'proj'), sep = '-', remove = FALSE)
  write.table(meta, file = 'GDCquery_results_allcancers.txt', col.names = T, row.names = F, quote = F, sep = '\t')

#### FFPE Clinical Data ####

ffpe <- lapply(projs, FUN = GDCquery_clinic, 'biospecimen') %>% 
  bind_rows() %>%
  select(sample_type_id, sample_id, sample_type, submitter_id, is_ffpe)
  write.table(ffpe, 'ffpe.txt', col.names = T, row.names = F, quote = F, sep= '\t')




