
library(edgeR)
library(hipathia)
library(Rtsne)
library(pheatmap)
library(RColorBrewer)
library(tidyverse)


######################
### pre-processing ###
######################

setwd("/media/mpaya/DATA/immunosarc/01_revisited/results")
if (!dir.exists("01_qc")) {dir.create("01_qc")}

# load data
pheno <- read_tsv("data_pheno.tsv", show_col_types = F) %>% as.data.frame
rownames(pheno) <- pheno$custom_id
htg <- read.table("data_htg.tsv", header = T, row.names = 1, sep = "\t")


## quality control with negative controls
# convert htg counts to log2(cpm) | y = log2(x + 3 / mean(colSums(cnt)) * 10^6)
cpm <- DGEList(counts = t(htg))
cpm <- calcNormFactors(cpm, method = "none")
cpm <- cpm(cpm, log = TRUE, prior.count = 3)
write_tsv(as.data.frame(t(cpm)) %>% rownames_to_column("custom_id"), "data_htg_cpm.tsv")

boxplot(as.vector(unlist(cpm)))
png("01_qc/boxplot_expr_cpm.png")
boxplot(cpm[-c(1:8),1:60])
dev.off()

## process negative control values
ant_vals <- t(cpm[1:4,]) %>% as.data.frame
ant_vals$avg <- rowMeans(ant_vals)
ant_vals$grand_avg <- mean(ant_vals$avg)
ant_vals$delta_mean <- ant_vals$avg - ant_vals$grand_avg
ant_vals$sd_dmean <- sd(ant_vals$delta_mean)
ant_vals$good <- abs(ant_vals$delta_mean) < 2*abs(ant_vals$sd_dmean)
table(ant_vals$good)  ## all 102 TRUE

## add metadata and plot points
ant_vals <- ant_vals %>% rownames_to_column("sample")
ant_vals <- cbind(ant_vals, pheno[ant_vals$sample,])
p <- ggplot(data = ant_vals) + 
  geom_point(aes(x = ID, y = delta_mean, shape = time)) + 
  geom_hline(yintercept = 2 * ant_vals$sd_dmean[1], color = "blue", linetype = "dashed") + 
  geom_hline(yintercept = -2 * ant_vals$sd_dmean[1], color = "blue", linetype = "dashed") + 
  scale_y_continuous(sec.axis = dup_axis(breaks = c(-2 * ant_vals$sd_dmean[1], 2 * ant_vals$sd_dmean[1]), 
                                         labels = c("-2SD", "2SD"), name = NULL)) + 
  ggtitle("Quality control chart for ANT negative controls of HTG")
ggsave("01_qc/plot_qc_ant.tiff", p, compression = "lzw")


###########
### TMM ###
###########
## normalization of HTG data with TMM
# control samples are removed before normalization

x <- t(htg[,-c(1:8)])
group <- factor(pheno[colnames(x),]$SUBTYPE)
y <- DGEList(counts = x, group = group)
keep <- filterByExpr(y)   # all 2567 probes are TRUE
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method = "TMM")

save(y, file = "expr_filt_tmm.Rdata")

y <- cpm(y, log=TRUE, prior.count=1)

png("01_qc/boxplot_expr_tmm.png")
boxplot(y[,1:60])
dev.off()
write_tsv(as.data.frame(t(y)) %>% rownames_to_column("custom_id"), "data_htg_tmm.tsv")



#####################
### VISUALIZATION ###
#####################

setwd("01_qc/")

### Apply z-scaling to data and draw heatmap and MDS

z_data <- as.data.frame(t(scale(t(y))))


# heatmap

annotation_col = data.frame(
  Batch = factor(pheno$batch),
  Time = factor(pheno$time),
  Progression = factor(pheno$Progression),
  Subtype = factor(pheno$SUBTYPE)
)
rownames(annotation_col) = rownames(pheno)

getColors <- colorRampPalette(brewer.pal(12, "Paired"))
colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colors1 <- colorBlindBlack8[1:length(levels(annotation_col$Batch))]
names(colors1) <- levels(annotation_col$Batch)
colors2 <- getColors(length(levels(annotation_col$Time)))
names(colors2) <- levels(annotation_col$Time)
prog_colors <- colorBlindBlack8[1:length(levels(annotation_col$Progression))]
names(prog_colors) <- levels(annotation_col$Progression)
subtype_colors <- getColors(length(unique(pheno$SUBTYPE)))
names(subtype_colors) <- unique(pheno$SUBTYPE)
annot_colors <- list(Batch = colors1, Time = colors2, Progression = prog_colors, 
                     Subtype = subtype_colors)

h <- pheatmap(z_data, annotation_col = annotation_col, annotation_colors =annot_colors,
              show_rownames = F, show_colnames = F, clustering_method = "complete", 
              clustering_distance_rows = "correlation", clustering_distance_cols = "correlation", 
              breaks = seq(-1.5, 1.5, length.out = 101), width = 14, height = 10,
              filename = "heatmap_TMM.png")

gc()

rowclus <- cutree(h$tree_row, 2)
colclus <- cutree(h$tree_col, 2)

table(pheno[names(colclus[colclus == 1]),]$SUBTYPE)
table(pheno[names(colclus[colclus == 2]),]$SUBTYPE)

save(rowclus, colclus, file = "heatmap_clusters.Rdata")


## check clusters

colclus_pheno <- pheno %>% 
  filter(custom_id %in% colnames(z_data)) %>% 
  add_column(colclus = colclus[as.character(.$custom_id)])
colclus_pheno %>% filter(Progression == 0) %>% group_by(colclus) %>% 
  summarize(n_resp = n())  # 1 9 (right); 2 13 (left)

clus_data <- z_data %>% rownames_to_column("symbol") %>% 
  pivot_longer(-symbol, names_to = "custom_id", values_to = "value") %>% 
  add_column(colclus = as.factor(colclus[.$custom_id])) %>% 
  add_column(rowclus = as.factor(rowclus[.$symbol]))

ggplot(clus_data) + 
  geom_violin(aes(colclus, value, fill = rowclus))

# t-SNE new classification

set.seed(142)
tSNE_fit <- z_data %>% t %>% Rtsne(check_duplicates = F)
tSNE_df <- tSNE_fit$Y %>% 
  as.data.frame() %>%
  rename(tSNE1 = "V1",
         tSNE2 = "V2") %>%
  cbind(pheno[colnames(z_data),])

tSNE_df %>%
  ggplot(aes(x = tSNE1, 
             y = tSNE2,
             color = SUBTYPE,
             shape = SUBTYPE, 
             size = time)) +
  geom_point() + 
  scale_color_manual(values = subtype_colors) + 
  scale_shape_manual(values = seq(1, 18)) + 
  scale_size_manual(values = c(4,2)) + 
  theme_bw() + 
  ggtitle("tSNE plot of sarcoma expression data")
ggsave("tSNE_subtype.png", width = 11, height = 9)

tSNE_df %>%
  ggplot(aes(x = tSNE1, 
             y = tSNE2,
             color = as.character(batch),
             shape = as.character(batch))) +
  geom_point() + 
  scale_shape_manual(values = seq(15, 27)) + 
  scale_color_manual(values = colors1) + 
  theme_bw() + 
  ggtitle("tSNE plot of sarcoma batch effect")
ggsave("tSNE_tmm_batch.png", width = 11, height = 9)


########################################
### VISUALIZATION with modeled batch ###
########################################

batch <- factor(pheno[colnames(y),]$batch)
y2 <- removeBatchEffect(y, batch)
z_data2 <- as.data.frame(t(scale(t(y2))))

set.seed(142)
tSNE_fit <- z_data2 %>% t %>% Rtsne(check_duplicates = F)
tSNE_df <- tSNE_fit$Y %>% 
  as.data.frame() %>%
  rename(tSNE1 = "V1",
         tSNE2 = "V2") %>%
  cbind(pheno[colnames(z_data),])

tSNE_df %>%
  ggplot(aes(x = tSNE1, 
             y = tSNE2,
             color = SUBTYPE,
             shape = SUBTYPE, 
             size = time)) +
  geom_point() + 
  scale_color_manual(values = subtype_colors) + 
  scale_shape_manual(values = seq(1, 18)) + 
  scale_size_manual(values = c(4,2)) + 
  theme_bw() + 
  ggtitle("tSNE plot of sarcoma expression data")
ggsave("tSNE_subtype_rmbatch.png", width = 11, height = 9)

dev.off()


#########################
### k-means / jaccard ###
#########################
# an option to measure differences in methods is to create clusters
# they are looked for similarity with subtypes
# find best k which maximizes jaccard scores

setwd("/media/mpaya/DATA/immunosarc/01_revisited/results")

pheno <- read_tsv("data_pheno.tsv", col_types = cols())
y <- read_tsv("data_htg_tmm.tsv", col_types = cols()) %>% 
  column_to_rownames("custom_id") %>% t
htg_pheno <- pheno[match(pheno$custom_id, colnames(y))]
all.equal(as.character(htg_pheno$custom_id), colnames(y))

z_data <- as.data.frame(t(scale(t(y))))

batch <- factor(pheno[colnames(y),]$batch)
y2 <- removeBatchEffect(y, batch)
z_data2 <- as.data.frame(t(scale(t(y2))))

setwd("01_qc/")

n_subtypes <- length(unique(htg_pheno$SUBTYPE))

jac_res <- mat.or.vec(nr = 0, nc = 11)
for (i in seq(n_subtypes-5, n_subtypes)) {
  set.seed(123)
  res.km <- kmeans(t(z_data), i, nstart = 25)
  jac_df <- tibble(k = i, Subtype = htg_pheno$SUBTYPE, cluster = res.km$cluster) %>% 
    add_count(cluster, name = "len_clu") %>% 
    add_count(Subtype, cluster, name = "intersect") %>% 
    add_count(Subtype, name = "len_sub") %>% unique %>% 
    mutate(perc = intersect/len_sub) %>% 
    mutate(union = len_clu + len_sub - intersect, jaccard = intersect/union) %>% 
    group_by(cluster) %>% mutate(max_clu = max(jaccard)) %>% ungroup %>% 
    group_by(Subtype) %>% mutate(max_sarc = max(jaccard)) %>% ungroup
  jac_res <- rbind(jac_res, jac_df)
}

jac_res2 <- mat.or.vec(nr = 0, nc = 11)
for (i in seq(n_subtypes-5, n_subtypes)) {
  set.seed(123)
  res.km <- kmeans(t(z_data2), i, nstart = 25)
  jac_df <- tibble(k = i, Subtype = htg_pheno$SUBTYPE, cluster = res.km$cluster) %>% 
    add_count(cluster, name = "len_clu") %>% 
    add_count(Subtype, cluster, name = "intersect") %>% 
    add_count(Subtype, name = "len_sub") %>% unique %>% 
    mutate(perc = intersect/len_sub) %>% 
    mutate(union = len_clu + len_sub - intersect, jaccard = intersect/union) %>% 
    group_by(cluster) %>% mutate(max_clu = max(jaccard)) %>% ungroup %>% 
    group_by(Subtype) %>% mutate(max_sarc = max(jaccard)) %>% ungroup
  jac_res2 <- rbind(jac_res2, jac_df)
}

jac_res <- add_column(jac_res, rmbatch = "no", .before = 1)
jac_res2 <- add_column(jac_res2, rmbatch = "yes", .before = 1)
jac_all <- rbind(jac_res, jac_res2)

write_tsv(jac_all, "jaccard_clusters_kmeans.tsv")

summ_res <- filter(jac_all, jaccard == max_sarc) %>% group_by(k, rmbatch) %>% 
  summarise(mean_perc = mean(perc), mean_jac = mean(jaccard)) %>% ungroup
summ_res%>% 
  pivot_longer(starts_with("mean"), names_to = "stat", values_to = "score") %>% 
  ggplot(.) +
  geom_line(aes(x = k, y = score, color = stat, linetype = rmbatch)) + 
  ggtitle("Jaccard score and % occupancy of 9 subtypes on k clusters")
ggsave("jaccard_plot_by_k.png")

filter(jac_res, jaccard == max_sarc, k == filter(summ_res, mean_jac == max(mean_jac))$k)

# max jaccard is at k = 8. I'll look at the grouping of samples
library(factoextra)
k <- filter(summ_res, mean_jac == max(mean_jac))$k
set.seed(123)
res.km <- kmeans(t(z_data), k, nstart = 25)

# Dimension reduction using PCA
res.pca <- prcomp(t(z_data),  scale = TRUE)
# Coordinates of individuals
ind.coord <- as.data.frame(get_pca_ind(res.pca)$coord)
# Add Species groups from the original data sett
ind.coord$Subtype <- htg_pheno$SUBTYPE
# Add clusters obtained using the K-means algorithm
ind.coord$cluster <- factor(res.km$cluster)

# Percentage of variance explained by dimensions
eigenvalue <- round(get_eigenvalue(res.pca), 1)
variance.percent <- eigenvalue$variance.percent
head(eigenvalue)

ggscatter(
  ind.coord, x = "Dim.1", y = "Dim.2", 
  color = "cluster", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "Subtype", size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", variance.percent[1], "% )" ),
  ylab = paste0("Dim 2 (", variance.percent[2], "% )" )
) +
  stat_mean(aes(color = cluster), size = 4) + 
  scale_shape_manual(values = seq(1:10)) + 
  ggtitle(paste0("PCA and k-means clustering of STS (",k," clusters)"))
ggsave("PCA_kmeans_STS.png", width = 8, height = 6)

ggscatter(
  ind.coord, x = "Dim.1", y = "Dim.2", 
  color = "cluster", palette = "npg", ellipse = TRUE, ellipse.type = "convex",
  shape = "cluster", size = 1.5,  legend = "right", ggtheme = theme_bw(),
  xlab = paste0("Dim 1 (", variance.percent[1], "% )" ),
  ylab = paste0("Dim 2 (", variance.percent[2], "% )" )
) +
  stat_mean(aes(color = cluster), size = 4) + 
  scale_shape_manual(values = seq(1:10)) + 
  ggtitle(paste0("PCA and k-means clustering of STS (",k," clusters)"))
ggsave("PCA_kmeans_cluster.png", width = 8, height = 6)

dev.off()
