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


#### Download gene-level copy number data ####
  
glcn <- GDCquery(project = projs, data.category = 'Copy Number Variation', workflow.type = 'ASCAT2')
extract_glcn <- glcn$results[[1]] %>% filter(data_type == 'Gene Level Copy Number')
glcn$results[[1]] <- extract_glcn #guess i already did know how to use subassignment
GDCdownload(glcn)

#### Use TCGAutils to connect filename --> case --> aliquot

all_files <- list.files(path = "/home/elle/PhD/cancer_wgd/GDCdata/", recursive = TRUE) #22208
glcn_files <- all_files[grep('gene_level_copy_number.tsv', all_files)] %>% #11104 
  as.data.frame() %>%
  separate(col = 1, sep = '/', into = c('Cancer', 'harmonized', 'CNV', 'GLCN', 'folder', 'file'), remove = FALSE) 

cases <- filenameToBarcode(glcn_files$file, legacy = FALSE) #22208
aliquots <- barcodeToUUID(cases$aliquots.submitter_id)

glcn_files <- glcn_files %>% left_join(aliquots, by = c('file' = 'file_name')) 
