#### Setup ####
library(TCGAutils)
library(data.table)
library(tidyverse)
library(janitor)
options(scipen = 999) #without this bedtools complains about the occasional scientific notation

#### Import GDC metadata ####

meta <- read.table('GDCquery_results_allcancers.txt', header = TRUE)
ffpe <- read.table('ffpe.txt', header = TRUE, sep = '\t')

#### Recommended exclusions by ASCAT team ####

ascat_contam_swap <- "TCGA-04-1351, TCGA-04-1371, TCGA-06-0178, TCGA-57-1586, TCGA-V5-AASX-b" 
ascat_germline <- "TCGA-05-4403, TCGA-06-2559, TCGA-12-3651, TCGA-25-1877, TCGA-26-5135, TCGA-2G-AAEW, TCGA-77-6842-a, TCGA-86-8054, TCGA-91-A4BD, TCGA-97-8179, TCGA-A5-A2K4, TCGA-A8-A07L, TCGA-AB-2812, TCGA-AB-2843, TCGA-AB-2847, TCGA-AB-2855, TCGA-AB-2856, TCGA-AB-2887, TCGA-AB-2891, TCGA-AB-2909, TCGA-AB-2913, TCGA-AB-2924, TCGA-AB-2944, TCGA-AB-2954, TCGA-AB-2997, TCGA-AB-3008, TCGA-AC-A2QJ, TCGA-AH-6547, TCGA-AP-A0LQ, TCGA-AY-4071, TCGA-BH-A1EW, TCGA-BJ-A28W-a, TCGA-BK-A139-e, TCGA-BM-6198, TCGA-BT-A20U, TCGA-C5-A1MK, TCGA-CG-5717, TCGA-CH-5753, TCGA-CH-5769, TCGA-CV-6441, TCGA-D8-A146, TCGA-D8-A1XC, TCGA-DD-AAEH, TCGA-DJ-A2Q2, TCGA-DM-A1D6, TCGA-E2-A15G, TCGA-EE-A29D, TCGA-EY-A1G8, TCGA-FD-A62N, TCGA-FS-A4F0, TCGA-GS-A9U4, TCGA-P5-A72U, TCGA-UZ-A9PJ" 
ascat_exclus <- "TCGA-05-4403, TCGA-06-2559, TCGA-12-3651, TCGA-25-1877, TCGA-26-5135, TCGA-2G-AAEW, TCGA-77-6842-a, TCGA-86-8054, TCGA-91-A4BD, TCGA-97-8179, TCGA-A5-A2K4, TCGA-A8-A07L, TCGA-AB-2812, TCGA-AB-2843, TCGA-AB-2847, TCGA-AB-2855, TCGA-AB-2856, TCGA-AB-2887, TCGA-AB-2891, TCGA-AB-2909, TCGA-AB-2913, TCGA-AB-2924, TCGA-AB-2944, TCGA-AB-2954, TCGA-AB-2997, TCGA-AB-3008, TCGA-AC-A2QJ, TCGA-AH-6547, TCGA-AP-A0LQ, TCGA-AY-4071, TCGA-BH-A1EW, TCGA-BJ-A28W-a, TCGA-BK-A139-e, TCGA-BM-6198, TCGA-BT-A20U, TCGA-C5-A1MK, TCGA-CG-5717, TCGA-CH-5753, TCGA-CH-5769, TCGA-CV-6441, TCGA-D8-A146, TCGA-D8-A1XC, TCGA-DD-AAEH, TCGA-DJ-A2Q2, TCGA-DM-A1D6, TCGA-E2-A15G, TCGA-EE-A29D, TCGA-EY-A1G8, TCGA-FD-A62N, TCGA-FS-A4F0, TCGA-GS-A9U4, TCGA-P5-A72U, TCGA-UZ-A9PJ, TCGA-04-1351, TCGA-04-1371, TCGA-06-0178, TCGA-57-1586, TCGA-V5-AASX-b" 
ascat_exclus <- gsub(',', '', trimws(scan(text = ascat_exclus, what = ','))) 

absolute <- fread('TCGA_mastercalls.abs_segtabs.fixed.txt')
absamples <- fread('TCGA_mastercalls.abs_tables_JSedit.fixed.txt')

absamples_filtered <- absamples %>% 
  clean_names() %>%
  filter(call_status == 'called') %>%
  #filter(genome_doublings == 1) %>% 
  separate(sample, into = c('TCGA', 'TSS', 'Individual', 'Sample', 'Portion', 'Plate', 'Center'), sep = '-', remove = FALSE) %>%
  separate(Sample, into = c('Sample.Type', 'Vial'), remove = FALSE, sep = 2) %>%
  mutate(Patient = paste(TCGA, TSS, Individual, sep = '-')) %>%
  mutate(Specimen = paste(Patient, Sample, sep = '-')) %>%
  filter(Sample.Type %in% c('01', '03', '09')) %>% #had 01, 02, 03, 05, 06
  filter(!Patient %in% ascat_exclus) %>%
  filter(Specimen %in% ffpe[ffpe$is_ffpe == FALSE,]$submitter_id) %>%  #remove FFPE samples
  arrange(Vial, desc(Plate), Portion) %>% 
  distinct(Patient, .keep_all = TRUE) 

write.table(absamples_filtered, 'absolute_samples_WGD_filtered.txt', sep = '\t', col.names = T, row.names = F, quote = F)
absamples_filtered$samplename <- absamples_filtered$array



filtered_absolute <- absolute %>%
  filter(Chromosome %in% 1:22) %>%
  filter(Sample %in% absamples_filtered$array) %>%
  mutate(Chromosome = paste0('chr', Chromosome)) %>%
  left_join(absamples_filtered, by = c('Sample' = 'array')) #635283

write.table(filtered_absolute, 'filtered_absolute_WGD.txt', sep = '\t', col.names = T, row.names = F, quote = F) #2883 samples

filtered_absolute %>%
  filter(genome_doublings == 1)%>%
  mutate(Start = Start - 1) %>%
  dplyr::select(Chromosome, Start, End, Sample, Modal_Total_CN) %>%
  write.table(file = 'absolute_filtered_WGD_CN.bed', quote = F, col.names = F, row.names = F, sep = '\t')


filtered_absolute %>%
  filter(genome_doublings == 0) %>%
  mutate(Start = Start - 1) %>%
  dplyr::select(Chromosome, Start, End, Sample, Modal_Total_CN) %>%
  write.table(file = 'absolute_filtered_diploid_CN.bed', quote = F, col.names = F, row.names = F, sep = '\t')


#time bedtools map -a genes_hg19_sorted.bed -b absolute_filtered_diploid_CN_sorted.bed -c 5,5,5,5,5,5 -o mean,median,min,max,count,stdev -f 0.9 -null -999 > mapped_hg19_absolute_diploid.bed







  




