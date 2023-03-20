library(tidyverse)
library(edgeR)
# library(hipathia)
library(survminer)
library(survival)

setwd("/media/mpaya/DATA/immunosarc/01_revisited/results")
if (!dir.exists("02_survival/survival_plots")) 
  {dir.create("02_survival/survival_plots", recursive = T)}

## read data
pheno <- read_tsv("data_pheno.tsv", col_types = cols())
expr <- read_tsv("data_htg_tmm.tsv", col_types = cols())
dens <- read_tsv("data_density.tsv", col_types = cols())
cyt <- read_tsv("data_cyt.tsv", col_types = cols())

#########################
### survival analysis ###
#########################

setwd("02_survival/")

## join all patient data at basal point

pheno_basal <- pheno %>% filter(!is.na(days)) %>% group_by(ID) %>% 
  mutate(max_days = max(days)) %>% ungroup %>% filter(time == "BASAL") %>% 
  add_column(colclus = colclus[as.character(.$custom_id)])
data_basal <- pheno_basal %>% 
  mutate(Death = as.numeric(Death)) %>% filter(!is.na(Death)) %>% 
  # inner_join(cyt, by = "custom_id") %>%
  # inner_join(dens, by = "custom_id") %>%
  inner_join(expr, by = "custom_id") %>%
  anti_join(select(., ID, time) %>% .[duplicated(.),], by = c("ID", "time"))
colnames(data_basal) <- make.names(colnames(data_basal))

ncol(expr) - length(intersect(colnames(expr), colnames(data_basal))) # identical names

## survival analysis

covariates <- colnames(data_basal)[15:ncol(data_basal)]
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(max_days, Death) ~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = data_basal)})

# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value <- signif(x$wald["pvalue"], digits = 3)
                         wald.test <- signif(x$wald["test"], digits = 3)
                         beta <- signif(x$coef[1], digits = 3);  # coefficient beta
                         HR <- signif(x$coef[2], digits = 3);  # exp(beta)
                         HR.ci.lower <- signif(x$conf.int[,"lower .95"], 3)
                         HR.ci.upper <- signif(x$conf.int[,"upper .95"], 3)
                         HR2 <- paste0(HR, " (", HR.ci.lower, "-", HR.ci.upper, ")")
                         conc <- signif(x$concordance["C"], digits = 3)
                         res <- c(beta, HR2, HR, HR.ci.lower, HR.ci.upper, wald.test, conc, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "HR", "HR_CI_lower",  
                                       "HR_CI_upper","wald.test", "concordance", "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
res$p.adj <- signif(p.adjust(res$p.value, method = "fdr"), 3)
res$p.zph <- unlist(lapply(univ_models, function(x) tryCatch(cox.zph(x)$table[1,3], error=function(err) NA)))


head(arrange(res, p.value))[,-2]
dim(filter(res, p.adj < 0.1))

res <- arrange(res, p.value) %>% rownames_to_column("symbol")
write_tsv(res, "surv_htg_basal.tsv")


### Survival plots of significant genes
expr_basal <- column_to_rownames(expr, "custom_id") %>% 
  .[as.character(data_basal$custom_id),]
z_data <- as.data.frame(t(scale(expr_basal)))
event_data <- apply(z_data, 1, function(x) 
  ifelse(x > 0.3, "High", ifelse(x < -0.3, "Low", NA))) %>% as.data.frame
event_data <- cbind(data_basal[,1:14], event_data)

table(unlist(event_data[,15:ncol(event_data)]))

signif <- filter(res, p.adj < 0.1)$symbol
for (s in signif) {
  gene_sym <- s
  out <- paste0("survival_plots/surv_coxsignif-", s, ".png")
  fit <- survfit(Surv(max_days, Death) ~ eval(parse(text = s)), data = event_data)
  p <- ggsurvplot(fit, data = event_data, conf.int = TRUE, pval = TRUE, 
                  surv.median.line = "hv", legend.labs = c("High expr", "Low expr"), 
                  palette = c("#E7B800", "#2E9FDF")) + ggtitle(s)
  ggsave(out, p$plot, width = 5.4, height = 4)
}


out <- "survival_plots/surv_sample_clusters_basal.png"
fit <- survfit(Surv(max_days, Death) ~ colclus, data = data_basal)
p <- ggsurvplot(fit, data = data_basal, conf.int = TRUE, pval = TRUE, 
                surv.median.line = "hv", palette = c("#E7B800", "#2E9FDF")) + 
  ggtitle("Survival by sample clusters: BASAL")
ggsave(out, p$plot, width = 5.4, height = 4)


#######################
### SURVIVAL AT W13 ###
#######################

pheno_w13 <- pheno %>% filter(!is.na(days)) %>% group_by(ID) %>% 
  mutate(max_days = max(days)) %>% ungroup %>% filter(time == "W13") %>% 
  add_column(colclus = colclus[as.character(.$custom_id)])
data_w13 <- pheno_w13 %>% 
  mutate(Death = as.numeric(Death)) %>% filter(!is.na(Death)) %>% 
  inner_join(expr, by = "custom_id") %>%
  anti_join(select(., ID, time) %>% .[duplicated(.),], by = c("ID", "time"))
colnames(data_w13) <- make.names(colnames(data_w13))

## survival analysis

covariates <- colnames(data_w13)[15:ncol(data_w13)]
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(max_days, Death) ~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = data_w13)})

# Extract data 
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value <- signif(x$wald["pvalue"], digits = 3)
                         wald.test <- signif(x$wald["test"], digits = 3)
                         beta <- signif(x$coef[1], digits = 3);  # coefficient beta
                         HR <- signif(x$coef[2], digits = 3);  # exp(beta)
                         HR.ci.lower <- signif(x$conf.int[,"lower .95"], 3)
                         HR.ci.upper <- signif(x$conf.int[,"upper .95"], 3)
                         HR2 <- paste0(HR, " (", HR.ci.lower, "-", HR.ci.upper, ")")
                         conc <- signif(x$concordance["C"], digits = 3)
                         res <- c(beta, HR2, HR, HR.ci.lower, HR.ci.upper, wald.test, conc, p.value)
                         names(res)<-c("beta", "HR (95% CI for HR)", "HR", "HR_CI_lower",  
                                       "HR_CI_upper","wald.test", "concordance", "p.value")
                         return(res)
                         #return(exp(cbind(coef(x),confint(x))))
                       })
res_w13 <- as.data.frame(t(as.data.frame(univ_results, check.names = FALSE)))
res_w13$p.adj <- signif(p.adjust(res_w13$p.value, method = "fdr"), 3)
res_w13$p.zph <- unlist(lapply(univ_models, function(x) tryCatch(cox.zph(x)$table[1,3], error=function(err) NA)))

head(arrange(res_w13, p.value))[,-2]
dim(filter(res_w13, p.value < 0.005))

res_w13 <- arrange(res_w13, p.value) %>% rownames_to_column("symbol")
write_tsv(res_w13, "surv_htg_w13.tsv")


### Survival plots of significant genes
expr_w13 <- column_to_rownames(expr, "custom_id")[as.character(data_w13$custom_id),]
z_data <- as.data.frame(t(scale(expr_w13)))
event_data <- apply(z_data, 1, function(x) 
  ifelse(x > 0.3, "High", ifelse(x < -0.3, "Low", NA))) %>% as.data.frame
event_data <- cbind(data_w13[,1:14], event_data)

table(unlist(event_data[,15:ncol(event_data)]))

signif <- filter(res_w13, p.value < 0.005)$symbol
for (s in signif) {
  gene_sym <- s
  out <- paste0("survival_plots/surv_psignif-w13-", s, ".png")
  fit <- survfit(Surv(max_days, Death) ~ eval(parse(text = s)), data = event_data)
  p <- ggsurvplot(fit, data = event_data, conf.int = TRUE, pval = TRUE, 
                  surv.median.line = "hv", legend.labs = c("High expr", "Low expr"), 
                  palette = c("#E7B800", "#2E9FDF")) + ggtitle(s)
  ggsave(out, p$plot, width = 5.4, height = 4)
}

out <- "survival_plots/surv_sample_clusters_w13.png"
fit <- survfit(Surv(max_days, Death) ~ colclus, data = data_w13)
p <- ggsurvplot(fit, data = data_w13, conf.int = TRUE, pval = TRUE, 
                surv.median.line = "hv", palette = c("#E7B800", "#2E9FDF")) + 
  ggtitle("Survival by sample clusters: W13")
ggsave(out, p$plot, width = 5.4, height = 4)

dev.off()


###################
### FOREST PLOT ###
###################
## join all data and create forest plot

univ_res <- rbind(res %>% add_column(week = "BASAL", .before = 1), 
                  res_w13 %>% add_column(week = "W13", .before = 1)) %>% 
  mutate(HR = as.numeric(HR), HR_CI_lower = as.numeric(HR_CI_lower), 
         HR_CI_upper = as.numeric(HR_CI_upper), week = factor(week, levels = c("W13", "BASAL")),
         predictor = factor(symbol, levels = unique(filter(., p.value < 0.005)$symbol)))
univ_res$HR_CI_upper[univ_res$HR_CI_upper > 10] = 10
predictors <- levels(univ_res$predictor)
p <- ggplot(data = univ_res %>% filter(predictor %in% predictors),
            aes(x = week, y = HR, ymin = HR_CI_lower, ymax = HR_CI_upper ))+
  geom_pointrange(aes(col = week)) +
  geom_hline(yintercept = 1, linetype = 2) +
  xlab('Time')+ ylab("Risk Ratio (95% Confidence Interval)") +
  geom_errorbar(aes(ymin=HR_CI_lower, ymax=HR_CI_upper, col = week), linewidth = 0.5, cex = 1) + 
  facet_wrap(~predictor, strip.position = "left", nrow = 9, scales = "free_y") + 
  scale_y_log10() + 
  coord_flip() + 
  scale_color_manual(values = c("firebrick", "lightblue3")) + 
  theme(plot.title=element_text(size = 16,face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(angle = 90, face="bold"),
        axis.title = element_text(size=12,face="bold"),
        strip.text.y = element_text(hjust=0,vjust = 1,angle=180,face="bold"))
p
ggsave("plots/univ_conf_int.png", p, width = 12, height = 7.5)

