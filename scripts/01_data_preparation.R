library(readxl)
library(edgeR)
# library(hipathia)
library(survminer)
library(survival)
library(tidyverse)

#################
### read data ###
#################

setwd("/media/mpaya/DATA/immunosarc/ARCHIVE/dataIMMUNOSARC")

## clinical results
clinical <- read_excel("Clinical Results GEIS-52.xlsx")
colnames(clinical) <- unlist(lapply(str_split(colnames(clinical), " "), "[", 1))

## cytokines
cyt1 <- read_excel("65Cyt_ALL_IMMUNOSARC_01072021 with days.xlsx", sheet = 1, col_names = F)
cyt2 <- read_excel("65Cyt_ALL_IMMUNOSARC_01072021 with days.xlsx", sheet = 2, col_names = F)
cyt3 <- read_excel("65Cyt_ALL_IMMUNOSARC_01072021 with days.xlsx", sheet = 3, col_names = F)
cyt4 <- read_excel("65Cyt_ALL_IMMUNOSARC_01072021 with days.xlsx", sheet = 4, col_names = F)

cyt_data <- rbind(filter(cyt1, !is.na(...1)), filter(cyt2, !is.na(...1)), 
                  filter(cyt3, !is.na(...1)), filter(cyt4, !is.na(...1)))
cyt_data <- cyt_data[,is.na(as.list(cyt1[1,]))]
cyt_data <- cyt_data[,2:ncol(cyt_data)]
colnames(cyt_data) <- c(as.character(cyt1[2,2:6]), cyt1[1,][!is.na(cyt1[1,])])
cyt_data <- cyt_data %>% rename(ID = `Sample ID`, days = `days from basal`) %>% 
  mutate(ID = as.numeric(ID)) %>% filter(!is.na(ID)) %>% 
  mutate(days = ifelse(is.na(days), 0, days))

## density
density <- read_excel("Sarcomas_IBis_Results.xlsx")
colnames(density) <- c("SUBTYPE", "ID", "time", 
                       paste0("density_", unlist(lapply(str_split(colnames(density)[4:10], " "), '[', 2))), 
                       paste0("perc_", unlist(lapply(str_split(colnames(density)[11:17], " "), '[', 1))))
density <- mutate(density, time = ifelse(time == "Basal", "BASAL", time))

## HTG gene expression
files <- list.files("HTG", full.names = TRUE)
htg1 <- read_excel(files[1], sheet = "Raw", col_names = F,  skip = 7) %>% t
htg2 <- read_excel(files[2], sheet = "Raw", col_names = F, skip = 7) %>% t
htg3 <- read_excel(files[3], sheet = "Raw", col_names = F, skip = 7) %>% t
htg4 <- read_excel(files[4], sheet = "Raw", col_names = F, skip = 7) %>% t
colnames(htg1) <- htg1[1,]
htg1 <- htg1[-1,]
htg1 <- htg1 %>% as.data.frame %>% rename(`HTG ID` = HTG) %>% 
  add_column(batch = c(rep(1,10), rep(2, 9), rep(3, 10), rep(4, 8)), .before = 7)
colnames(htg2) <- htg2[1,]
htg2 <- add_column(as.data.frame(htg2), batch = 5, .before = 7)
htg2 <- htg2[-1,colnames(htg1)]
colnames(htg3) <- htg3[1,]
htg3 <- add_column(as.data.frame(htg3), batch = 6, .before = 7)
htg3 <- htg3[-1,colnames(htg1)]
colnames(htg4) <- htg4[1,]
htg4 <- add_column(as.data.frame(htg4), batch = 7, .before = 7)
htg4 <- htg4[-1,colnames(htg1)]

all_htg <- rbind(htg1, htg2, htg3, htg4) %>% .[,-c(1,2,3,8)] %>% 
  rename(ID = `GEIS-52 ID`, time = Time) %>%  # select(-contains("CTRL")) %>% 
  mutate(ID = as.numeric(ID)) %>% filter(!is.na(ID)) %>% 
  mutate(time = ifelse(time == "Basal", "BASAL", "W13"))

## RETURN TO PROJECT FOLDER

setwd("/media/mpaya/DATA/immunosarc/01_revisited/results")


## join pheno data
cyt_pheno <- cyt_data[1:5]
dens_pheno <- density[1:3]
htg_pheno <- all_htg[1:4]
all_pheno <- full_join(htg_pheno, cyt_pheno, by = c("ID", "time")) %>% 
  full_join(dens_pheno, by = c("ID", "time")) %>% arrange(ID, time) %>% 
  add_column(custom_id = seq(1, nrow(.)), .before = 1)

### fill missing values
type_pairs <- select(all_pheno, diagnosis, SUBTYPE) %>% unique %>% 
  filter(!is.na(diagnosis) & !is.na(SUBTYPE))
type_names <- type_pairs$diagnosis
names(type_names) <- type_pairs$SUBTYPE
all_pheno <- all_pheno %>% 
  mutate(diagnosis = ifelse(SUBTYPE %in% type_pairs$SUBTYPE, type_names[SUBTYPE], diagnosis))

patient_type <- select(all_pheno, ID, diagnosis) %>% unique %>% filter(!is.na(diagnosis))
type_names <- patient_type$diagnosis
names(type_names) <- patient_type$ID
all_pheno <- all_pheno %>% 
  mutate(diagnosis = ifelse(ID %in% patient_type$ID, type_names[as.character(ID)], diagnosis))

type_names <- type_pairs$SUBTYPE
names(type_names) <- type_pairs$diagnosis
all_pheno <- all_pheno %>% 
  mutate(SUBTYPE = ifelse(diagnosis %in% type_pairs$diagnosis, type_names[diagnosis], SUBTYPE))

all_pheno <- left_join(all_pheno, clinical, by = "ID")

write_tsv(all_pheno, "data_pheno.tsv")  # further manual corrections on subtypes
write_tsv(clinical, "data_clinical.tsv")


## add custom ids to data tables

### 65 cyt
cyt_fixed <- inner_join(select(all_pheno, ID, days, custom_id), cyt_data, 
                        by = c("ID", "days")) %>% 
  select(-ID, -days, -time, -R_RECIST, -diagnosis) %>% 
  apply(., 2, function(x) str_replace(x, "<=0", "0")) %>% 
  as.data.frame
write_tsv(cyt_fixed, "data_cyt.tsv")

### density
density_fixed <- inner_join(select(all_pheno, ID, time, custom_id), density, 
                            by = c("ID", "time")) %>% 
  select(-ID, -time, -SUBTYPE)
write_tsv(density_fixed, "data_density.tsv")

## HTG
htg_fixed <- inner_join(select(all_pheno, `Sample Name`, time, custom_id), 
                        all_htg, by = c("Sample Name", "time")) %>% 
  select(-c(1,2,4,5)) %>% mutate_all(as.numeric)

write_tsv(htg_fixed, "data_htg.tsv")

