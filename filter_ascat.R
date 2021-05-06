#### Setup ####
library(TCGAutils)
library(data.table)
library(tidyverse)

#### Import GDC metadata ####

meta <- read.table('GDCquery_results_allcancers.txt', header = TRUE)
ffpe <- read.table('ffpe.txt', header = TRUE)

#### Recommended exclusions by ASCAT team ####

ascat_contam_swap <- "TCGA-04-1351, TCGA-04-1371, TCGA-06-0178, TCGA-57-1586, TCGA-V5-AASX-b" 
ascat_germline <- "TCGA-05-4403, TCGA-06-2559, TCGA-12-3651, TCGA-25-1877, TCGA-26-5135, TCGA-2G-AAEW, TCGA-77-6842-a, TCGA-86-8054, TCGA-91-A4BD, TCGA-97-8179, TCGA-A5-A2K4, TCGA-A8-A07L, TCGA-AB-2812, TCGA-AB-2843, TCGA-AB-2847, TCGA-AB-2855, TCGA-AB-2856, TCGA-AB-2887, TCGA-AB-2891, TCGA-AB-2909, TCGA-AB-2913, TCGA-AB-2924, TCGA-AB-2944, TCGA-AB-2954, TCGA-AB-2997, TCGA-AB-3008, TCGA-AC-A2QJ, TCGA-AH-6547, TCGA-AP-A0LQ, TCGA-AY-4071, TCGA-BH-A1EW, TCGA-BJ-A28W-a, TCGA-BK-A139-e, TCGA-BM-6198, TCGA-BT-A20U, TCGA-C5-A1MK, TCGA-CG-5717, TCGA-CH-5753, TCGA-CH-5769, TCGA-CV-6441, TCGA-D8-A146, TCGA-D8-A1XC, TCGA-DD-AAEH, TCGA-DJ-A2Q2, TCGA-DM-A1D6, TCGA-E2-A15G, TCGA-EE-A29D, TCGA-EY-A1G8, TCGA-FD-A62N, TCGA-FS-A4F0, TCGA-GS-A9U4, TCGA-P5-A72U, TCGA-UZ-A9PJ" 
ascat_exclus <- "TCGA-05-4403, TCGA-06-2559, TCGA-12-3651, TCGA-25-1877, TCGA-26-5135, TCGA-2G-AAEW, TCGA-77-6842-a, TCGA-86-8054, TCGA-91-A4BD, TCGA-97-8179, TCGA-A5-A2K4, TCGA-A8-A07L, TCGA-AB-2812, TCGA-AB-2843, TCGA-AB-2847, TCGA-AB-2855, TCGA-AB-2856, TCGA-AB-2887, TCGA-AB-2891, TCGA-AB-2909, TCGA-AB-2913, TCGA-AB-2924, TCGA-AB-2944, TCGA-AB-2954, TCGA-AB-2997, TCGA-AB-3008, TCGA-AC-A2QJ, TCGA-AH-6547, TCGA-AP-A0LQ, TCGA-AY-4071, TCGA-BH-A1EW, TCGA-BJ-A28W-a, TCGA-BK-A139-e, TCGA-BM-6198, TCGA-BT-A20U, TCGA-C5-A1MK, TCGA-CG-5717, TCGA-CH-5753, TCGA-CH-5769, TCGA-CV-6441, TCGA-D8-A146, TCGA-D8-A1XC, TCGA-DD-AAEH, TCGA-DJ-A2Q2, TCGA-DM-A1D6, TCGA-E2-A15G, TCGA-EE-A29D, TCGA-EY-A1G8, TCGA-FD-A62N, TCGA-FS-A4F0, TCGA-GS-A9U4, TCGA-P5-A72U, TCGA-UZ-A9PJ, TCGA-04-1351, TCGA-04-1371, TCGA-06-0178, TCGA-57-1586, TCGA-V5-AASX-b" 
ascat_exclus <- gsub(',', '', trimws(scan(text = ascat_exclus, what = ','))) 

#### Import and format combined ASCAT segment data ####

ascat <- fread('~/cat_all_AS.txt', header = TRUE) %>%
  filter(GDC_Aliquot != 'GDC_Aliquot')

#convert UUIDs to barcodes using TCGAutils to query GDC (?)
aliquots <- unique(ascat$GDC_Aliquot) #11104 rows
codes <- UUIDtoBarcode(aliquots, from_type = 'aliquot_ids') #11104 rows

aliquots <- codes %>%
  separate(col = 'portions.analytes.aliquots.submitter_id', into = c('TCGA', 'TSS', 'Individual', 'Sample', 'Portion', 'Plate', 'Center'), sep = '-', remove = FALSE) %>%
  mutate(Patient = paste(TCGA, TSS, Individual)) %>%
  mutate(Specimen = paste(Patient, Sample)) %>%
  separate(col = Sample, into = c('Sample.Type', 'Vial'), sep = 2) #11104 rows

#### Filter ASCAT samples ####

filtered_codes <- codes %>%
    filter(Sample.Type %in% c('01', '03', '09')) %>% #remove metastatic and recurrent; keep solid primary tumour and primary blood cancers
    filter(!Patient %in% ascat_exclus) %>% 
    filter(Specimen %in% ffpe[ffpe$is_ffpe == FALSE,]$submitter_id) %>%  
    arrange(Vial, desc(Plate), Portion, by.group = TRUE) %>%  
    distinct(Patient, .keep_all = TRUE)  
  return(asc) 
} 

filters <- filter_barcodes()

  