#### Setup ####
library(TCGAutils)
library(data.table)
library(tidyverse)
options(scipen = 999) #without this bedtools complains about the occasional scientific notation

#### Import GDC metadata ####

meta <- read.table('GDCquery_results_allcancers.txt', header = TRUE)
ffpe <- read.table('ffpe.txt', header = TRUE, sep = '\t')

#### Recommended exclusions by ASCAT team ####

ascat_contam_swap <- "TCGA-04-1351, TCGA-04-1371, TCGA-06-0178, TCGA-57-1586, TCGA-V5-AASX-b" 
ascat_germline <- "TCGA-05-4403, TCGA-06-2559, TCGA-12-3651, TCGA-25-1877, TCGA-26-5135, TCGA-2G-AAEW, TCGA-77-6842-a, TCGA-86-8054, TCGA-91-A4BD, TCGA-97-8179, TCGA-A5-A2K4, TCGA-A8-A07L, TCGA-AB-2812, TCGA-AB-2843, TCGA-AB-2847, TCGA-AB-2855, TCGA-AB-2856, TCGA-AB-2887, TCGA-AB-2891, TCGA-AB-2909, TCGA-AB-2913, TCGA-AB-2924, TCGA-AB-2944, TCGA-AB-2954, TCGA-AB-2997, TCGA-AB-3008, TCGA-AC-A2QJ, TCGA-AH-6547, TCGA-AP-A0LQ, TCGA-AY-4071, TCGA-BH-A1EW, TCGA-BJ-A28W-a, TCGA-BK-A139-e, TCGA-BM-6198, TCGA-BT-A20U, TCGA-C5-A1MK, TCGA-CG-5717, TCGA-CH-5753, TCGA-CH-5769, TCGA-CV-6441, TCGA-D8-A146, TCGA-D8-A1XC, TCGA-DD-AAEH, TCGA-DJ-A2Q2, TCGA-DM-A1D6, TCGA-E2-A15G, TCGA-EE-A29D, TCGA-EY-A1G8, TCGA-FD-A62N, TCGA-FS-A4F0, TCGA-GS-A9U4, TCGA-P5-A72U, TCGA-UZ-A9PJ" 
ascat_exclus <- "TCGA-05-4403, TCGA-06-2559, TCGA-12-3651, TCGA-25-1877, TCGA-26-5135, TCGA-2G-AAEW, TCGA-77-6842-a, TCGA-86-8054, TCGA-91-A4BD, TCGA-97-8179, TCGA-A5-A2K4, TCGA-A8-A07L, TCGA-AB-2812, TCGA-AB-2843, TCGA-AB-2847, TCGA-AB-2855, TCGA-AB-2856, TCGA-AB-2887, TCGA-AB-2891, TCGA-AB-2909, TCGA-AB-2913, TCGA-AB-2924, TCGA-AB-2944, TCGA-AB-2954, TCGA-AB-2997, TCGA-AB-3008, TCGA-AC-A2QJ, TCGA-AH-6547, TCGA-AP-A0LQ, TCGA-AY-4071, TCGA-BH-A1EW, TCGA-BJ-A28W-a, TCGA-BK-A139-e, TCGA-BM-6198, TCGA-BT-A20U, TCGA-C5-A1MK, TCGA-CG-5717, TCGA-CH-5753, TCGA-CH-5769, TCGA-CV-6441, TCGA-D8-A146, TCGA-D8-A1XC, TCGA-DD-AAEH, TCGA-DJ-A2Q2, TCGA-DM-A1D6, TCGA-E2-A15G, TCGA-EE-A29D, TCGA-EY-A1G8, TCGA-FD-A62N, TCGA-FS-A4F0, TCGA-GS-A9U4, TCGA-P5-A72U, TCGA-UZ-A9PJ, TCGA-04-1351, TCGA-04-1371, TCGA-06-0178, TCGA-57-1586, TCGA-V5-AASX-b" 
ascat_exclus <- gsub(',', '', trimws(scan(text = ascat_exclus, what = ','))) 

#### Import and format combined ASCAT segment data ####

ascat <- fread('~/cat_all_AS.txt', header = TRUE) %>%
  filter(GDC_Aliquot != 'GDC_Aliquot') #remove interspersed headers from cat

#convert UUIDs to barcodes using TCGAutils to query GDC (?)
aliquots <- unique(ascat$GDC_Aliquot) #11104 rows
codes <- UUIDtoBarcode(aliquots, from_type = 'aliquot_ids') #11104 rows

codes <- codes %>%
  separate(col = 'portions.analytes.aliquots.submitter_id', into = c('TCGA', 'TSS', 'Individual', 'Sample', 'Portion', 'Plate', 'Center'), sep = '-', remove = FALSE) %>%
  mutate(Patient = paste(TCGA, TSS, Individual, sep = '-')) %>%
  mutate(Specimen = paste(Patient, Sample, sep = '-')) %>%
  separate(col = Sample, into = c('Sample.Type', 'Vial'), sep = 2) #11104 rows

#### Filter ASCAT samples ####

filtered_codes <- codes %>%
    filter(Sample.Type %in% c('01', '03', '09')) %>% #remove metastatic and recurrent; keep solid primary tumour and primary blood cancers
    filter(!Patient %in% ascat_exclus) %>% 
    filter(Specimen %in% ffpe[ffpe$is_ffpe == FALSE,]$submitter_id) %>%  #remove FFPE samples
    arrange(Vial, desc(Plate), Portion) %>%  #Plate is sorted in descending order according to standard
    left_join(meta, by = c('Patient')) %>%
    distinct(Patient, .keep_all = TRUE)  

filtered_ascat <- ascat %>%
  filter(!Chromosome %in% c('chrX', 'chrY')) %>% #remove sex chromosomes
  filter(GDC_Aliquot %in% filtered_codes$portions.analytes.aliquots.aliquot_id) %>% #keep only preferred samples
  mutate_at(c('Start', 'End', 'Copy_Number', 'Major_Copy_Number', 'Minor_Copy_Number'), as.numeric)  %>% #convert these columns to numeric
  left_join(filtered_codes, by = c('GDC_Aliquot' = 'portions.analytes.aliquots.aliquot_id')) #get barcodes

write.table(filtered_ascat, 'filtered_ascat.txt', sep = '\t', col.names = T, row.names = F, quote = F) #10346 samples 


#### Sample Ploidy ####

sample_avg = filtered_ascat %>%
  mutate(len = (End - Start) + 1) %>%
  mutate(prod = len * Copy_Number) %>%
  group_by(GDC_Aliquot, cases, proj) %>%
  summarize(mean_CN = sum(prod)/sum(len)) 
sample_avg$index <- 1:nrow(sample_avg)

ggplot(sample_avg, aes(x = mean_CN, fill = proj)) + geom_density() + facet_wrap(~proj, scales = 'free_y') + theme(legend.position = 'none') + ylab('') + ggtitle('Sample Ploidy') + labs(subtitle = 'ASCAT Copy Number data from GDC')
ggsave('sample_average_ASCAT.png') 
ggplot(sample_avg, aes(x = mean_CN)) + geom_density() + ylab('') + ggtitle('Sample Ploidy (All Cancers)') + labs(subtitle = 'ASCAT Copy Number data from GDC')
ggsave('sample_average_ASCAT_pancancer.png') 
write.table(sample_avg, file = 'sample_average_ploidy_ASCAT.txt', sep = '\t', quote = F, col.names = T, row.names = F)


#### BED format ####

filtered_ascat %>%
  mutate(Start = Start - 1) %>%
  left_join(sample_avg, by = c('GDC_Aliquot')) %>%
  mutate(sample_info = paste(Copy_Number, proj.x, round(mean_CN,3), sep = '_')) %>%
  dplyr::select(Chromosome, Start, End, sample_info, Copy_Number) %>%
  write.table(file = 'ascat_filtered.bed', quote = F, col.names = F, row.names = F, sep = '\t')
  

#### Make bed file with binary-encoded CN, recoded sample ID, and only WGD samples ####

filtered_ascat %>%
  mutate(Start = Start - 1) %>%
  left_join(sample_avg, by = c('GDC_Aliquot')) %>%
  mutate(sample_info = paste(Copy_Number, proj.x, round(mean_CN,3), sep = '_')) %>%
  filter(mean_CN >= 2.7) %>%
  dplyr::select(Chromosome, Start, End, index, Copy_Number) %>%
  mutate(Copy_Number = ifelse(Copy_Number >= 3, 1, 0)) %>%
  write.table(file = 'ascat_filtered_binary.bed', quote = F, col.names = F, row.names = F, sep = '\t')


  