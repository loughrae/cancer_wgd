library(data.table)
library(ggplot2)
library(tidyverse)

genes <- fread('genes_hg19.tsv')

gn <- fread('~/Downloads/all_samples.consensus_CN.by_gene.170214.txt.gz', header = T, sep = '\t') 
gn <- gn[, -c('Locus ID','Cytoband')] %>% 
  as.data.frame() %>%
  separate(col = 1, into = c('Gene', 'Chr'), sep = '\\|', fill = 'right', remove = TRUE) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  left_join(genes[, c('chromosome_name', 'external_gene_name')], by = c('Gene' = 'external_gene_name')) %>%
  filter(chromosome_name %in% 1:22) %>% #no sex chrs
  dplyr::select(-c(chromosome_name, Chr))

#samples <- gn %>% dplyr::select(contains('-')) %>% colnames()

sam <- read.table('~/Downloads/pcawg_sample_sheet.tsv', header = T, sep = '\t')
sam <- separate(sam, col = dcc_project_code, into = c('cancer', 'place'), sep = '-')

ploi <- read.table('~/Downloads/consensus.20170217.purity.ploidy.txt.gz', header = T)

samples <- sapply(gn[, 2:ncol(gn)], function(x) sum(is.na(x)))
samples <- as.data.frame(samples)
samples$Sample <- rownames(samples)
names (samples)[1] <- 'nas' #not using the letters NA lol
#should i put a limit on number of NAs? like a hard one?
include <- samples %>%
  left_join(sam, by = c('Sample' = 'aliquot_id')) %>% 
  filter(donor_wgs_exclusion_white_gray == 'Whitelist') %>%
  filter(library_strategy == 'WGS') %>% #not RNA-Seq
  filter(grepl('Primary tumour', dcc_specimen_type, fixed = TRUE)) %>%
  dplyr::select(Sample, nas, donor_unique_id, cancer, dcc_specimen_type) %>%
  left_join(ploi, by = c('Sample' = 'samplename')) %>%
  filter(wgd_uncertain == FALSE) %>%
  dplyr::select(-purity_conf_mad, -wgd_uncertain) %>%
  arrange(nas) %>%
  distinct(donor_unique_id, .keep_all = TRUE)

filt_samples <- include$Sample
all.equal(include$Sample, unique(include$Sample))
gn <- gn %>% select(Gene, all_of(filt_samples))

long <- pivot_longer(gn, cols = !Gene, names_to = 'Sample', values_to = 'CN') %>%
  filter(!is.na(CN)) %>%
  left_join(sam, by = c('Sample' = 'aliquot_id')) %>% 
  dplyr::select(Gene, Sample, CN, donor_unique_id, cancer, dcc_specimen_type) %>%
  left_join(ploi, by = c('Sample' = 'samplename')) %>%
  dplyr::select(-purity_conf_mad, -wgd_uncertain) 

fwrite(long, file = 'pcawg_long.txt', sep = '\t', col.names = T, row.names = F, quote = F)

NROW(unique(long$Sample)) #2349
NROW(unique(long$donor_unique_id)) #2349

#re,pve x amd y!! or else female cancers could get higher CNs



sample_means <- long %>% group_by(cancer, Sample) %>% summarize(mean = mean(CN, na.rm = T), med = median(CN, na.rm = T), std = sd(CN, na.rm = T))
samples_plot <- ggplot(sample_means, aes(x = mean)) + geom_density() + facet_wrap(~cancer, scales = 'free_y') + xlab('Mean Gene Copy Number') + ylab('Density (independent y scales)') + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) + ggtitle('Mean Copy Number of PCAWG Samples')
ggsave(samples_plot, file = 'sample_means_plot.png')
write.table(sample_means, file = 'pcawg_sample_means.txt', sep = '\t', quote = F, col.names = T, row.names = F)
