
library(hipathia)
library(org.Hs.eg.db)
library(limma)
library(clusterProfiler)
library(tidyverse)
library(ggVennDiagram)
library(RColorBrewer)

setwd("/media/mpaya/DATA/immunosarc/01_revisited/results")

## read data
pheno <- read_tsv("data_pheno.tsv", col_types = cols())
expr <- read_tsv("data_htg_tmm.tsv", col_types = cols())
expr <- as.data.frame(expr) %>% column_to_rownames("custom_id") %>% t

path_list <- read_tsv( "/media/mpaya/DATA/data/physiological_paths.tsv", 
                       col_types = cols(), col_names = FALSE)
pathways <- load_pathways(species = "hsa", path_list[[2]])

if (!dir.exists("04_drug")) {dir.create("04_drug")}
setwd("04_drug/")


# prepare expression table

sunit <- c("PDGFRB","FLT1","KIT","KDR","FLT4","FLT3","CSF1R","PDGFRA")

ko <- 0.6

ko_cnt <- cnt[,as.character(filter(pheno, time == "BASAL", !is.na(`Sample Name`))$custom_id)]
ko_06 <- ko_cnt
colnames(ko_06) <- paste0(colnames(ko_cnt), "_ko06")
ko_06[sunit,] <- ko * ko_06[sunit,]
ko_01 <- ko_cnt
colnames(ko_01) <- paste0(colnames(ko_cnt), "_ko01")
ko_01[sunit,] <- 0.1 * ko_01[sunit,]

colnames(cnt) <- paste0(colnames(cnt), "_noko")

cnt_all <- cbind(cnt, ko_01, ko_06)

# prepare and run hipathia
cnt_all <- translate_data(cnt_all, species = "hsa")
cnt_all <- normalize_data(cnt_all)
results <- hipathia(cnt_all, pathways)
path_vals <- get_paths_data(results, matrix = T)
path_vals <- normalize_paths(path_vals, pathways)

write.table(path_vals, "tmm_ko01_ko06-path_vals.tsv", sep = "\t", col.names = NA)
path_vals <- read.table("tmm_ko01_ko06-path_vals.tsv", header = T, sep = '\t', row.names = 1) %>% 
  rename_all(.funs = function(x) str_extract(x, "[0-9].*"))

## prepare pheno data with ko information
ko_samples <- str_split_fixed(colnames(path_vals), "_", 2) %>% as.data.frame
colnames(ko_samples) <- c("custom_id", "ko")
ko_samples$ko_id <- colnames(path_vals)
ko_pheno <- left_join(ko_samples, mutate(pheno, custom_id = as.character(custom_id)), by = "custom_id")

## sigmoidal curve of drug perturbation
ko_patients <- filter(ko_pheno, ko == "ko06") %>% pull(custom_id)
sig_data <- path_vals %>% as.data.frame %>% rownames_to_column("circuit") %>% 
  pivot_longer(-circuit, names_to = "ko_id", values_to = "activity") %>% 
  inner_join(ko_pheno, by = "ko_id") %>% filter(custom_id %in% ko_patients) %>% 
  select(-ko_id) %>% pivot_wider(names_from = "ko", values_from = "activity") %>% 
  mutate(logFC = log2(ko06) - log2(noko)) %>% filter(logFC != 0) %>% 
  group_by(ID) %>% mutate(mean_fc = mean(abs(logFC))) %>% ungroup %>% 
  mutate(SUBTYPE = ifelse(is.na(SUBTYPE), "ND", SUBTYPE)) %>% 
  mutate(ID = factor(ID, levels = arrange(unique(select(., ID, mean_fc)), mean_fc)$ID))

ggplot(sig_data) + 
  geom_point(aes(x = ID, y = mean_fc, color = Progression, shape = SUBTYPE)) + 
  scale_color_brewer(palette = "Paired") + 
  scale_shape_manual(values = seq(1,length(unique(sig_data$SUBTYPE)))) + 
  ggtitle("Mean response of patients to in silico treatment with sunitinib")

ggsave("sigmoidal_curve_basalVSko06.png", width = 7.5, height = 4.5)


ggplot(sig_data %>% mutate(ID = factor(ID, levels = rev(levels(ID))))) + 
  geom_point(aes(x = ID, y = mean_fc, color = Progression, shape = SUBTYPE)) + 
  scale_color_brewer(palette = "Paired") + 
  scale_shape_manual(values = seq(1,length(unique(sig_data$SUBTYPE)))) + 
  ggtitle("Mean response of patients to in silico treatment with sunitinib")

ggsave("sigmoidal_curve_basalVSko06_desc.png", width = 7.5, height = 4.5)


## plot distribution of annotations

data_sub <- select(sig_data, ID, mean_fc) %>% unique %>% arrange(desc(mean_fc))
rank_fc <- data_sub$mean_fc
names (rank_fc) <- data_sub$ID 

term2gene <- unique(select(sig_data, ID, Progression))
term2gene <- data.frame(term = term2gene[,2], gene = term2gene[,1])

gse <- GSEA(
  rank_fc,
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = term2gene, 
  scoreType = "pos"
)

gtitle <- paste("Progression", gse@result$ID[1])
p <- gseaplot(gse, by = "runningScore", title = gtitle, geneSetID = 1) + 
  xlab("Position in the Ranked List of Patients") + scale_x_continuous(expand = c(0,0))
ggsave("gsea_sigmoidal_ko06VSnoko_prog1.png", p, width = 8, height = 2)
ggsave("gsea_sigmoidal_ko06VSnoko_prog1_yaxis.png", p, width = 8.2, height = 3)

gtitle <- paste("Progression", gse@result$ID[2])
p <- gseaplot(gse, by = "runningScore", title = gtitle, geneSetID = 2) + 
  xlab("Position in the Ranked List of Patients") + scale_x_continuous(expand = c(0,0))
ggsave("gsea_sigmoidal_ko06VSnoko_prog0.png", p, width = 8.05, height = 2)



## sigmoidal curve of actual response

sig_data2 <- path_vals %>% as.data.frame %>% rownames_to_column("circuit") %>% 
  pivot_longer(-circuit, names_to = "ko_id", values_to = "activity") %>% 
  inner_join(select(ko_pheno, -custom_id, -`Sample Name`, -batch, -days, -R_RECIST), by = "ko_id") %>% 
  filter(ko == "noko") %>% select(-ko_id, -ko) %>% 
  pivot_wider(names_from = "time", values_from = "activity", values_fn = function(x) mean(x)) %>% 
  filter(!is.na(BASAL) & !is.na(W13)) %>% mutate(OS = as.numeric(OS), PFS = as.numeric(PFS)) %>% 
  mutate(cat_pfs = ifelse(PFS < median(PFS, na.rm = T), "low PFS", "high PFS")) %>% 
  mutate(cat_pfs = ifelse(is.na(cat_pfs), "ND", cat_pfs)) %>% 
  mutate(cat_OS = ifelse(OS < median(OS, na.rm = T), "low OS", "high OS")) %>% 
  mutate(cat_OS = ifelse(is.na(cat_OS), "ND", cat_OS)) %>% 
  mutate(logFC = log2(W13 + 0.5) - log2(BASAL + 0.5)) %>% filter(logFC != 0) %>% 
  group_by(ID) %>% mutate(mean_fc = mean(abs(logFC))) %>% ungroup %>% 
  mutate(SUBTYPE = ifelse(is.na(SUBTYPE), "ND", SUBTYPE)) %>% 
  mutate(ID = factor(ID, levels = arrange(unique(select(., ID, mean_fc)), mean_fc)$ID))

ggplot(sig_data2) + 
  geom_point(aes(x = ID, y = mean_fc, color = Progression, shape = SUBTYPE)) + 
  scale_color_brewer(palette = "Paired") + 
  scale_shape_manual(values = seq(1,length(unique(sig_data$SUBTYPE)))) + 
  ggtitle("Mean response of patients to treatment with sunitinib")

ggsave("sigmoidal_curve_basalVSw13.png", width = 6.5, height = 4.5)

ggplot(sig_data2) + 
  geom_point(aes(x = ID, y = mean_fc, color = cat_pfs, shape = cat_OS)) + 
  scale_color_brewer(palette = "Paired") + 
  scale_shape_manual(values = c(3, 4, 5)) + 
  ggtitle("Mean response of patients to treatment with sunitinib")

ggsave("sigmoidal_curve_basalVSw13_surv.png", width = 6.5, height = 4.5)

## plot distribution of annotations

data_sub <- select(sig_data2, ID, mean_fc) %>% unique %>% arrange(desc(mean_fc))
rank_fc <- data_sub$mean_fc
names (rank_fc) <- data_sub$ID 

term2gene <- unique(select(sig_data2, ID, Progression))
term2gene <- data.frame(term = term2gene[,2], gene = term2gene[,1])

gse <- GSEA(
  rank_fc,
  minGSSize = 5,
  maxGSSize = 500,
  pvalueCutoff = 1,
  pAdjustMethod = "BH",
  TERM2GENE = term2gene, 
  scoreType = "pos"
)

gtitle <- paste("Progression", gse@result$ID[1])
p <- gseaplot(gse, by = "runningScore", title = gtitle, geneSetID = 1) + 
  xlab("Position in the Ranked List of Patients") + scale_x_continuous(expand = c(0,0))
ggsave("gsea_sigmoidal_basalVSw13_prog1.png", p, width = 8, height = 2)
ggsave("gsea_sigmoidal_basalVSw13_prog1_yaxis.png", p, width = 8.2, height = 3)

gtitle <- paste("Progression", gse@result$ID[2])
p <- gseaplot(gse, by = "runningScore", title = gtitle, geneSetID = 2) + 
  xlab("Position in the Ranked List of Patients") + scale_x_continuous(expand = c(0,0))
ggsave("gsea_sigmoidal_basalVSw13_prog0.png", p, width = 8.05, height = 2)


# compare simulated and real responses

bad_resp <- filter(sig_data, mean_fc < 0.15) %>% pull(ID) %>% unique
good_resp <- filter(sig_data, mean_fc > 0.15) %>% pull(ID) %>% unique

bad_resp2 <- filter(sig_data2, mean_fc < 0.06) %>% pull(ID) %>% unique
good_resp2 <- filter(sig_data2, mean_fc > 0.06) %>% pull(ID) %>% unique

length(intersect(good_resp, good_resp2))
length(intersect(bad_resp, bad_resp2))
length(intersect(good_resp, bad_resp2))
length(intersect(bad_resp, good_resp2))

x <- list(real_high = good_resp2, real_low = bad_resp2, pred_high = good_resp, pred_low = bad_resp)
x<- lapply(x, as.character)
p <- ggVennDiagram(x, label_alpha = 0, label = "count") + 
  ggplot2::scale_fill_gradient(low = "blue", high = "yellow")
p <- p + ggtitle("Comparison of real and predicted response level")
p
ggsave("venn-real_pred_resp_lvl.png", p, width = 4, height = 4)



## see differences by circuit

circ_diff <- sig_data %>% group_by(circuit) %>% mutate(mean_circ = mean(logFC)) %>% 
  ungroup %>% arrange(mean_circ) %>% 
  mutate(circuit = factor(circuit, levels = unique(circuit)))

ggplot(circ_diff) + 
  geom_boxplot(aes(x = circuit, y = logFC), outlier.size = 0.5) + 
  theme(axis.text.x=element_blank()) +
  ggtitle("Mean response by circuit to in silico treatment with sunitinib 0.6")
ggsave("circuitFC_drug06.png", width = 6, height = 4)


circ_diff01 <- path_vals %>% as.data.frame %>% rownames_to_column("circuit") %>% 
  pivot_longer(-circuit, names_to = "ko_id", values_to = "activity") %>% 
  inner_join(ko_pheno, by = "ko_id") %>% filter(custom_id %in% ko_patients) %>% 
  select(-ko_id) %>% pivot_wider(names_from = "ko", values_from = "activity") %>% 
  mutate(logFC = log2(ko01) - log2(noko)) %>% filter(logFC != 0) %>% 
  group_by(circuit) %>% mutate(mean_circ = mean(logFC)) %>% ungroup %>% 
  arrange(mean_circ) %>% mutate(circuit = factor(circuit, levels = unique(circuit)))

ggplot(circ_diff01) + 
  geom_boxplot(aes(x = circuit, y = logFC), outlier.size = 0.5) + 
  theme(axis.text.x=element_blank()) +
  ggtitle("Mean response by circuit to in silico treatment with sunitinib 0.1")
ggsave("circuitFC_drug01.png", width = 6, height = 4)


## compare circuit changes by response

sig_data4 <- path_vals %>% as.data.frame %>% rownames_to_column("circuit") %>% 
  pivot_longer(-circuit, names_to = "ko_id", values_to = "activity") %>% 
  inner_join(select(ko_pheno, -custom_id, -`Sample Name`, -batch, -days, -R_RECIST), by = "ko_id") %>% 
  filter(ko == "noko") %>% select(-ko_id, -ko) %>% 
  pivot_wider(names_from = "time", values_from = "activity", values_fn = function(x) mean(x)) %>% 
  full_join(
    select(sig_data, circuit, ID, noko, ko01, ko06) %>% mutate(ID = as.numeric(ID)), 
    by = c("circuit", "ID")
    )

dat4.1 <- filter(sig_data4, !is.na(BASAL) & !is.na(W13) & !is.na(noko)) %>% 
  mutate(real_fc = log2(W13 + 0.5) - log2(BASAL + 0.5), ko_fc = log2(ko06) - log2(noko)) %>% 
  group_by(ID) %>% mutate(mean_fc = mean(abs(real_fc)), mean_fcko = mean(abs(ko_fc))) %>% ungroup %>% 
  mutate(SUBTYPE = ifelse(is.na(SUBTYPE), "ND", SUBTYPE)) %>% 
  mutate(ID = factor(ID, levels = arrange(unique(select(., ID, mean_fc)), mean_fc)$ID))

ggplot(dat4.1) + 
  geom_point(aes(x = real_fc, y = ko_fc), pch = 1, alpha = 0.5) + 
  geom_abline(intercept = 0, slope = 1, color = "blue", size = 1) + 
  ggtitle("Relation of observed FC and predicted FC (ko  * 0.6) with sunitinib", 
          paste("Plotting", length(unique(dat4.1$circuit)), "circuits in", length(unique(dat4.1$ID)), "patients."))

ggsave("circuitFC_realVSko.png")

ggplot(dat4.1) + 
  geom_point(aes(x = ID, y = mean_fc, color = Progression, shape = SUBTYPE)) + 
  geom_point(aes(x = ID, y = mean_fcko, color = Progression, shape = SUBTYPE), 
             size = 0.8, show.legend = F) + 
  scale_color_brewer(palette = "Paired") + 
  scale_shape_manual(values = seq(1,length(unique(dat4.1$SUBTYPE)))) + 
  ggtitle("Mean response of patients to treatment with sunitinib", 
          "Ordered by observed patient FC and added FC of predicted kos")

ggsave("sigmoidal_curve_realVSko.png", width = 6.5, height = 4.5)

dev.off()
