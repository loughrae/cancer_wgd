library(TCGAbiolinks)
library(tidyverse)
allproj <- getGDCprojects()
projs <- allproj[startsWith(allproj$id, 'TCGA'),]$id
dall <- GDCquery(project = projs, data.category = 'Copy Number Variation', workflow.type = 'ASCAT2')
glcn <- dall
ascs <- dall
glc <- glcn$results[[1]] %>% filter(data_type == 'Gene Level Copy Number')
glcn$results[[1]] <- glc
asc <- ascs$results[[1]] %>% filter(data_type == 'Allele-specific Copy Number Segment')
ascs$results[[1]] <- asc
GDCdownload(ascs)
#GDCdownload(glcn)