
library(hipathia)
library(org.Hs.eg.db)
library(limma)
library(tidyverse)

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

##############################
### join tables + hipathia ###
##############################

# separate groups by time and progression
control_basal <- filter(pheno, time == "BASAL", Progression == 0) %>% pull(custom_id)
control_w13 <- filter(pheno, time == "W13", Progression == 0) %>% pull(custom_id)
nonresp_basal <- filter(pheno, time == "BASAL", Progression == 1) %>% pull(custom_id)
nonresp_w13 <- filter(pheno, time == "W13", Progression == 1) %>% pull(custom_id)

# prepare expression tables for each group
sunit <- c("PDGFRB","FLT1","KIT","KDR","FLT4","FLT3","CSF1R","PDGFRA")
ko <- 0.6

## expression of respondents
ctrl13_expr <- expr[,intersect(control_w13, colnames(expr))]
colnames(ctrl13_expr) <- paste0("ctrl13_", colnames(ctrl13_expr))
ctrlb_expr <- expr[,intersect(control_basal, colnames(expr))]
drugctrl_expr <- ctrlb_expr
drugctrl_expr[sunit,] <- ko * drugctrl_expr[sunit,]
colnames(ctrlb_expr) <- paste0("ctrlB_", colnames(ctrlb_expr))
colnames(drugctrl_expr) <- paste0("drugresp_", colnames(drugctrl_expr))

## add a smaller ko to get DACs
exp_ko01 <- expr[,intersect(control_basal, colnames(expr))]
exp_ko01[sunit,] <- 0.1 * expr[sunit,intersect(control_basal, colnames(expr))]
colnames(exp_ko01) <- paste0("drug01_", intersect(control_basal, colnames(expr)))

## expression of non-respondents
nr13_expr <- expr[,intersect(nonresp_w13, colnames(expr))]
colnames(nr13_expr) <- paste0("nr13_", colnames(nr13_expr))
nrb_expr <- expr[,intersect(nonresp_basal, colnames(expr))]
drugnr_expr <- nrb_expr
drugnr_expr[sunit,] <- ko * drugnr_expr[sunit,]
colnames(nrb_expr) <- paste0("nrB_", colnames(nrb_expr))
colnames(drugnr_expr) <- paste0("drugnr_", colnames(drugnr_expr))

## join all tables
sarc_drug <- cbind(ctrl13_expr, ctrlb_expr, drugctrl_expr, exp_ko01, nr13_expr, 
                   nrb_expr, drugnr_expr)


# prepare and run hipathia
sarc_drug <- translate_data(sarc_drug, species = "hsa")
sarc_drug <- normalize_data(sarc_drug)
results <- hipathia(sarc_drug, pathways)
path_vals <- get_paths_data(results, matrix = T)
path_vals <- normalize_paths(path_vals, pathways)

write.table(path_vals, "hipathia-path_vals.tsv", sep = "\t", col.names = NA)

length(intersect(rownames(sarc_drug), pathways$all.genes))  # 1144/2516

# plot activity on perturbed circuits

types <- unlist(lapply(str_split(colnames(path_vals), pattern = "_"), "[", 1))

# compare non-treated respondents with drug simulation
# this way we intend to see the effect of the drug

levels <- c("ctrlB", "drug01")
sel_vals <- path_vals[,colnames(path_vals)[types %in% levels]]
sel_group <- factor(types[types %in% levels], levels = levels)
des <- model.matrix(~sel_group)
fit <- lmFit(sel_vals, des)
fit <- eBayes(fit)
tts_drug <- topTable(fit, coef = ncol(des), number = nrow(sel_vals), sort.by = "none")


paths <- rownames(filter(tts_drug, adj.P.Val < 0.05))
path_names <- get_path_names(pathways, paths) %>% str_replace(":", ":\n")
names(path_names) <- paths
levels <- c("ctrlB", "ctrl13", "drugresp", "nrB", "nr13", "drugnr")
vals <- t(path_vals)[types %in% levels, paths] %>% as.data.frame %>% 
  rownames_to_column("sample") %>% 
  add_column(type = factor(types[types %in% levels], levels = levels)) %>% 
  pivot_longer(starts_with("P-"), names_to = "pathway", values_to = "activity") %>% 
  add_column(pathway_names = path_names[.$pathway]) %>% group_by(pathway) %>% 
  mutate(minmax_activ = (activity - min(activity)) / (max(activity) - min(activity)))

p_title <- paste("Circuit activity after perturbation with Sunitinib *", ko, "on all patients")
p <- ggplot(vals) + 
  geom_boxplot(aes(x = pathway_names, y = minmax_activ, fill = type), outlier.size = 0.1) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank(), 
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle(p_title) + 
  scale_fill_discrete(name = "Case",
                      breaks = levels,
                      labels = c("Respondents - BASAL", "Respondents - W13", 
                                 "Respondents - Drug", "Non resp - BASAL", 
                                 "Non resp - W13", "Non resp - Drug"))
ggsave("plot_effects_sunitinib_all_patients.png", p, width = 12, height = 6)




#####################
### drug gradient ###
#####################

setwd("/media/mpaya/DATA/immunosarc/01_revisited/results")

# read data

pheno <- read_tsv("data_pheno.tsv", col_types = cols())
expr <- read.table("data_htg_tmm.tsv", header = T, row.names = 1, sep = "\t") %>% t

path_list <- read_tsv( "/media/mpaya/DATA/data/physiological_paths.tsv", 
                       col_types = cols(), col_names = FALSE)
pathways <- load_pathways(species = "hsa", path_list[[2]])

setwd("04_drug/")

# separate groups by time and progression
control_ids <- filter(pheno, time == "W13", Progression == 0) %>% pull(custom_id)
control_basal <- filter(pheno, time == "BASAL", Progression == 0) %>% pull(custom_id)

# prepare expression tables for each group
sunit <- c("PDGFRB","FLT1","KIT","KDR","FLT4","FLT3","CSF1R","PDGFRA")

## expression of respondents
ctrl_expr <- expr[,intersect(control_ids, colnames(expr))]
colnames(ctrl_expr) <- paste0("ctrl13_", colnames(ctrl_expr))
ctrlb_expr <- expr[,intersect(control_basal, colnames(expr))]
colnames(ctrlb_expr) <- paste0("ctrlB_", colnames(ctrlb_expr))

exp_ko001 <- expr[,intersect(control_basal, colnames(expr))]
exp_ko001[sunit,] <- 0.01 * expr[sunit,intersect(control_basal, colnames(expr))]
colnames(exp_ko001) <- paste0("drug001_", intersect(control_basal, colnames(expr)))
exp_ko01 <- expr[,intersect(control_basal, colnames(expr))]
exp_ko01[sunit,] <- 0.1 * expr[sunit,intersect(control_basal, colnames(expr))]
colnames(exp_ko01) <- paste0("drug01_", intersect(control_basal, colnames(expr)))
exp_ko05 <- expr[,intersect(control_basal, colnames(expr))]
exp_ko05[sunit,] <- 0.5 * expr[sunit,intersect(control_basal, colnames(expr))]
colnames(exp_ko05) <- paste0("drug05_", intersect(control_basal, colnames(expr)))
exp_ko08 <- expr[,intersect(control_basal, colnames(expr))]
exp_ko08[sunit,] <- 0.8 * expr[sunit,intersect(control_basal, colnames(expr))]
colnames(exp_ko08) <- paste0("drug08_", intersect(control_basal, colnames(expr)))

sarc_drug <- cbind(ctrl_expr, ctrlb_expr, exp_ko001, exp_ko01, exp_ko05, exp_ko08)

# prepare and run hipathia
sarc_drug <- translate_data(sarc_drug, species = "hsa")
sarc_drug <- normalize_data(sarc_drug)
results <- hipathia(sarc_drug, pathways)
path_vals <- get_paths_data(results, matrix = T)
path_vals <- normalize_paths(path_vals, pathways)

types <- unlist(lapply(str_split(colnames(path_vals), pattern = "_"), "[", 1))

# compare non-treated respondents with drug simulation
# this way we intend to see the effect of the drug
levels <- c("ctrlB", "drug001")
sel_vals <- path_vals[,colnames(path_vals)[types %in% levels]]
sel_group <- factor(types[types %in% levels], levels = levels)
des <- model.matrix(~sel_group)
fit <- lmFit(sel_vals, des)
fit <- eBayes(fit)
tts_drug <- topTable(fit, coef = ncol(des), number = nrow(sel_vals), sort.by = "none")


paths <- rownames(filter(tts_drug, adj.P.Val < 0.05))
path_names <- get_path_names(pathways, paths) %>% str_replace(":", ":\n")
names(path_names) <- paths
levels <- c("ctrlB", "ctrl13", rev(unique(types)[3:6]))
vals <- t(path_vals)[types %in% levels, paths] %>% as.data.frame %>% 
  rownames_to_column("sample") %>% 
  add_column(type = factor(types[types %in% levels], levels = levels)) %>% 
  pivot_longer(starts_with("P-"), names_to = "pathway", values_to = "activity") %>% 
  add_column(pathway_names = path_names[.$pathway]) %>% group_by(pathway) %>% 
  mutate(minmax_activ = (activity - min(activity)) / (max(activity) - min(activity)))

p <- ggplot(vals) + 
  geom_boxplot(aes(x = pathway_names, y = minmax_activ, fill = type), outlier.size = 0.1) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.title = element_blank(), 
        plot.title = element_text(hjust = 0.5)) + 
  ggtitle("Circuit activity after gradient perturbation with Sunitinib on respondents") + 
  scale_fill_discrete(name = "Case",
                      breaks = levels,
                      labels = c("Respondents - BASAL", "Respondents - W13", 
                                 "Drug 0.8", "Drug 0.5", "Drug 0.1", "Drug 0.01"))
ggsave("plot_effects_sunitinib_controls_gradient.png", p, width = 10, height = 6)



