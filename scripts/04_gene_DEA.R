
library(org.Hs.eg.db)
library(limma)
library(clusterProfiler)
library(enrichR)
library(tidyverse)

setwd("/media/mpaya/DATA/immunosarc/01_revisited/results")


## read data
pheno <- read_tsv("data_pheno.tsv", col_types = cols())
expr <- read_tsv("data_htg_tmm.tsv", col_types = cols())

if (!dir.exists("03_genes")) {dir.create("03_genes")}
setwd("03_genes/")
if (!dir.exists("func_annot")) {dir.create("func_annot")}


########################
### tmm and model ID ###
########################
# paired DEA by time

load("expr_filt_tmm.Rdata")
htg_pheno <- filter(pheno, custom_id %in% colnames(y))
all.equal(htg_pheno$custom_id, colnames(y))

Subject <- factor(htg_pheno$ID)
Time <- factor(htg_pheno$time, levels=c("BASAL","W13"))
design <- model.matrix(~ Subject + Time)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit)
tts_time <- topTags(qlf, n = Inf) %>% as.data.frame

deg_time <- rownames(filter(tts_time, FDR < 0.05))   # 13

write_tsv(rownames_to_column(tts_time, "symbol"), "tts_time.tsv")
write_tsv(tibble(symbol = deg_time), "geneList-deg_time.tsv")

# functional enrichment

enriched <- enrichr(deg_time, dbs)
write_tsv(enriched$GO_Biological_Process_2021, "func_annot/enrichR-table-time.tsv")
maxchar <- max(nchar(enriched$GO_Biological_Process_2021$Term[1:25]))
w <- maxchar * 8 + 380
p <- plotEnrich(enriched[[1]], showTerms = 25, numChar = maxchar, y = "Count", 
                title = paste("GO enrichment of W13 vs BASAL")) + 
  theme(axis.text = element_text(size = 15), plot.title = element_text(size = 20))
ggsave("func_annot/enrichR-DEG-paired-time.png", p, scale = 4, 
       width = w, height = 700, units = "px")



# DEA by time and progression status

time <- factor(htg_pheno$time)
prog <- factor(htg_pheno$Progression)
batch <- factor(htg_pheno$batch)
timeprog <- factor(paste0(time, prog))
ID <- factor(htg_pheno$ID)

design <- model.matrix(~ 0 + timeprog + ID)
colnames(design)[1:length(levels(timeprog))] <- levels(timeprog)

## do contrasts
contrast.matrix <- makeContrasts(W130 - BASAL0, W131 - BASAL1, BASAL1 - BASAL0, 
                                 W131 - W130, levels = design) 

## fit

cnt <- expr %>% column_to_rownames("custom_id") %>% t

fit <- lmFit(cnt, design)
fit <- contrasts.fit(fit, contrast.matrix) 
fit <- eBayes(fit)
tts_1_resp <- topTable(fit, coef = 1, number = Inf) 
tts_2_prog <- topTable(fit, coef = 2, number = Inf) 
tts_3_bas <- topTable(fit, coef = 3, number = Inf) 
tts_4_w13 <- topTable(fit, coef = 4, number = Inf) 

deg_1_resp <- rownames(filter(tts_1_resp, adj.P.Val < 0.05))  # 0
deg_2_prog <- rownames(filter(tts_2_prog, adj.P.Val < 0.05))  # 1
deg_3_bas <- rownames(filter(tts_3_bas, adj.P.Val < 0.05))    # 517
deg_4_w13 <- rownames(filter(tts_4_w13, adj.P.Val < 0.05))    # 824

write_tsv(tibble(symbol = deg_1_resp), "geneList-deg_1_resp.tsv")
write_tsv(tibble(symbol = deg_2_prog), "geneList-deg_2_prog.tsv")
write_tsv(tibble(symbol = deg_3_bas), "geneList-deg_3_bas.tsv")
write_tsv(tibble(symbol = deg_4_w13), "geneList-deg_4_w13.tsv")

## functional enrichment

enriched <- enrichr(deg_3_bas, dbs)
write_tsv(enriched$GO_Biological_Process_2021, "func_annot/enrichR-table-basal.tsv")
maxchar <- max(nchar(enriched$GO_Biological_Process_2021$Term[1:25]))
w <- maxchar * 8 + 380
p <- plotEnrich(enriched[[1]], showTerms = 25, numChar = maxchar, y = "Count", 
                title = paste("GO enrichment of DEGs at basal time")) + 
  theme(axis.text = element_text(size = 15), plot.title = element_text(size = 20))
ggsave("func_annot/enrichR-DEG-basal.png", p, scale = 4, 
       width = w, height = 700, units = "px")

enriched <- enrichr(deg_4_w13, dbs)
write_tsv(enriched$GO_Biological_Process_2021, "func_annot/enrichR-table-w13.tsv")
maxchar <- max(nchar(enriched$GO_Biological_Process_2021$Term[1:25]))
w <- maxchar * 8 + 380
p <- plotEnrich(enriched[[1]], showTerms = 25, numChar = maxchar, y = "Count", 
                title = paste("GO enrichment of DEGs at w13")) + 
  theme(axis.text = element_text(size = 15), plot.title = element_text(size = 20))
ggsave("func_annot/enrichR-DEG-w13.png", p, scale = 4, 
       width = w, height = 700, units = "px")


#############################
### RESULTS ON 2 CLUSTERS ###
#############################
# On the heatmap, 2 clear patient and gene clusters were formed
# here, with DEA and functional enrichment we'll see their functions

setwd("/media/mpaya/DATA/immunosarc/01_revisited/results/03_genes/")
load(file = "../01_qc/heatmap_clusters.Rdata")
colclus <- as.factor(colclus)

design <- model.matrix(~ 0 + colclus)
colnames(design) <- c("clus1", "clus2")

contrast.matrix <- makeContrasts(clus2 - clus1, levels = design) 

fit <- lmFit(cnt, design)
fit <- contrasts.fit(fit, contrast.matrix) 
fit <- eBayes(fit)
tts_clus <- topTable(fit, coef = 1, number = Inf) 
deg_clus <- rownames(filter(tts_clus, adj.P.Val < 0.05))  # 2050/2559

write_tsv(tibble(symbol = deg_clus), "geneList-deg_clus.tsv")

design <- model.matrix(~ colclus)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit)
tts_clus <- topTags(qlf, n = Inf) %>% as.data.frame
deg_clus <- rownames(filter(tts_clus, FDR < 0.05))  # 1986/2559

write_tsv(rownames_to_column(tts_clus, "symbol"), "tts_clus.tsv")
write_tsv(tibble(symbol = deg_clus), "geneList-deg_clus.tsv")

# functional enrichment of DEGs by gene clusters

degclus <- rowclus[deg_clus]
degclus_1 <- degclus[degclus == 1]
degclus_2 <- degclus[degclus == 2]

# dbs <- listEnrichrDbs()  # 189 x 6 table
dbs <- "GO_Biological_Process_2021"
moredbs <- c(
  "KEGG_2021_Human", "GO_Molecular_Function_2021", "GO_Biological_Process_2021", 
  "GWAS_Catalog_2019", "ClinVar_2019", "DisGeNET", "Human_Phenotype_Ontology"
)

enriched <- enrichr(names(degclus_1), dbs)
write_tsv(enriched$GO_Biological_Process_2021, "func_annot/enrichR-table-cluster1.tsv")
maxchar <- max(nchar(enriched$GO_Biological_Process_2021$Term[1:25]))
w <- maxchar * 8 + 380
p <- plotEnrich(enriched[[1]], showTerms = 25, numChar = maxchar, y = "Count", 
                title = paste("GO enrichment of cluster 1 DEGs")) + 
  theme(axis.text = element_text(size = 15), plot.title = element_text(size = 20))
ggsave("func_annot/enrichR-DEG-cluster1.png", p, scale = 4, 
       width = w, height = 700, units = "px")

enriched <- enrichr(names(degclus_2), dbs)
write_tsv(enriched$GO_Biological_Process_2021, "func_annot/enrichR-table-cluster2.tsv")
maxchar <- max(nchar(enriched$GO_Biological_Process_2021$Term[1:25]))
w <- maxchar * 8 + 380
p <- plotEnrich(enriched[[1]], showTerms = 25, numChar = maxchar, y = "Count", 
                title = paste("GO enrichment of cluster 2 DEGs")) + 
  theme(axis.text = element_text(size = 15), plot.title = element_text(size = 20))
ggsave("func_annot/enrichR-DEG-cluster2.png", p, scale = 4, 
       width = w, height = 700, units = "px")


# SAVE RESULTS

save(list = ls(pattern = "^tts"), file = "DEA_tts_res.Rdata")

load("DEA_tts_res.Rdata")
