library(biomaRt)
library(tidyverse)

#### download genes ####
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl") 

geneset <- function(hg) {
  if (hg == 'hg19') {
    ensembl = useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl', host = 'grch37.ensembl.org')
  }
  genes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'chromosome_name', 'start_position', 'end_position', 'ensembl_peptide_id', 'gene_biotype'), mart = ensembl) %>%
    filter(chromosome_name %in% c(1:22)) %>%
    mutate(protein = case_when(ensembl_peptide_id != '' ~ 'Coding', TRUE ~ 'Non-Coding')) %>%
    filter(external_gene_name != '') %>%
    distinct(external_gene_name, .keep_all = TRUE)
  write.table(genes, file = paste0('genes_', hg, '.tsv'), row.names = F, col.names = T, quote = F, sep = '\t')
  
  genes %>% 
    mutate(chromosome_name = paste0('chr', chromosome_name)) %>%
    mutate(start_position = start_position - 1) %>%
    dplyr::select(chromosome_name, start_position, end_position, ensembl_gene_id) %>% 
    write.table(file = paste0('genes_', hg, '.bed'), row.names = F, col.names = F, quote = F, sep = '\t')
    
  return(genes)  
}

hg38 <- geneset('hg38')
hg19 <- geneset('hg19')


