library(biomaRt)
library(tidyverse)
library(janitor)
library(data.table)

### annotations ####
classB <- read.table('~/Downloads/ClassBGenes.txt')

ohnos <- read.table('~/Downloads/strict_ohnologs_full.csv', sep = '\t', header = T) #not actually comma-separated
oh <- unlist(ohnos)
uniq_oh <- unique(oh)

DBOs <- read.table('~/Downloads/DBO_list.txt')

species <- read.table('~/Downloads/speciesCopyNumberAnalysis.1-Y.tsv', header = T, sep = '\t')
conserved <- species[species$Unchanged == 13,]$Ensembl.Gene.ID
#### download genes ####

geneset <- function(hg) {
  if (hg == 'hg19') {
    genes <- fread('~/Downloads/mart_export(2).txt') %>% clean_names() 
    names(genes) <- c('ensembl_gene_id', 'chromosome_name', 'external_gene_name', 'start_position', 'end_position', 'gene_biotype')
    }
  if (hg == 'hg38') {
    ensembl = useEnsembl(biomart = 'genes', dataset = 'hsapiens_gene_ensembl')
    genes <- getBM(attributes = c('ensembl_gene_id', 'external_gene_name', 'chromosome_name', 'start_position', 'end_position', 'gene_biotype'), mart = ensembl) 
  }
    genes <- genes %>%
    filter(chromosome_name %in% c(1:22)) %>%
    mutate(type = case_when(gene_biotype == 'protein_coding' ~ 'Coding', endsWith(gene_biotype, 'pseudogene') ~ 'Pseudogene', endsWith(gene_biotype, 'RNA') ~ 'RNA', TRUE ~ 'Other')) %>%
    distinct(ensembl_gene_id, .keep_all = TRUE) %>%
    mutate(Ohnolog = case_when(ensembl_gene_id %in% uniq_oh ~ 'Ohnolog', TRUE ~ 'Non-Ohnolog')) %>%
    mutate(BCNV = case_when(ensembl_gene_id %in% classB$V1 ~ 'Class B', TRUE ~ 'Not Class B')) %>%
    mutate(DBO = case_when(ensembl_gene_id %in% DBOs$V1 ~ 'DBO', TRUE ~ 'Non-DBO')) %>%
    mutate(copyconserved = case_when(ensembl_gene_id %in% conserved ~ 'Conserved', TRUE ~ 'Not Conserved'))
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


