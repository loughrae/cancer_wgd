library(tidyverse)
library(TCGAutils)

avg <- fread('sample_average_ploidy_ASCAT.txt', header = T, sep = '\t')



all_files <- list.files(path = "/home/elle/PhD/cancer_wgd/GDCdata/", recursive = TRUE) #22208
glcn_files <- all_files[grep('gene_level_copy_number.tsv', all_files)] %>% #11104 
   as.data.frame() %>%
   separate(col = 1, sep = '/', into = c('Cancer', 'harmonized', 'CNV', 'GLCN', 'folder', 'file'), remove = FALSE) 

cases <- filenameToBarcode(glcn_files$file, legacy = FALSE) #22208
aliquots <- barcodeToUUID(cases$aliquots.submitter_id)

glcn_files <- glcn_files %>% left_join(aliquots, by = c('file' = 'file_name')) 


   mutate(filter = ifelse(file %in% avg$file_name, 'Yes', 'No'))

st = Sys.time()
for (i in 1:1000) {
  x = fread(paste0("/home/elle/PhD/cancer_wgd/GDCdata/", al2[i]))
  print(mean(x$copy_number, na.rm = T))
  x <- x %>% filter(chromosome %in% paste0('chr', 1:22))
}
en = Sys.time()