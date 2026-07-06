##############################################################################
# quick_setup.R
#
# One-shot script: train model + embed into package source + install.
# Run this from an R session where all original data files are available.
#
# After running, install with:
#   devtools::install("/Users/HXE/Desktop/MMimmuneRisk")
##############################################################################

library(survival)
library(survminer)
library(randomForestSRC)
library(plsRcox)
library(rms)
library(dplyr)

setwd("/Users/HXE/Desktop/机器学习鉴别功能性高危患者")

# ---- 1. Load data ----
load("原始数据2/model_data_fpkm.Rdata")
load("原始数据2/deg_gene.Rdata")
load("原始数据/JSPH-MM.Rdata")
load("原始数据2/exp_test2_cel.Rdata")
ptm_gene <- openxlsx::read.xlsx("原始数据2/2499_IRG.xlsx", sheet = 2)
load("原始数据2/train_sample.Rdata")
rm(test_Data2, exp_test2)

# ---- 2. Prepare train/test split ----
test_Data2 <- train %>% filter(!sample %in% train_sample)
train      <- train %>% filter(sample %in% train_sample)
exp_test2  <- exp_train[, test_Data2$sample]
exp_train_full <- exp_train[, train$sample]

exp_train_full <- exp_train_full[rownames(exp_train_full) %in% ptm_gene$gene, ]
exp_test1 <- exp_test1[rownames(exp_test1) %in% ptm_gene$gene, ]
exp_test2 <- exp_test2[rownames(exp_test2) %in% ptm_gene$gene, ]
exp_test3 <- exp_test3[rownames(exp_test3) %in% ptm_gene$gene, ]

intersects_fn <- function(...) Reduce(base::intersect, list(...))
intersects <- intersects_fn(
  rownames(exp_train_full), rownames(exp_test1),
  rownames(exp_test2), rownames(exp_test3), deg_gene
)

expr <- as.data.frame(t(exp_train_full[intersects, ])) %>%
  tibble::rownames_to_column("sample")
df <- train %>% inner_join(expr) %>% tibble::column_to_rownames("sample")

# ---- 3. Univariate Cox filter ----
pfilter <- 0.05
uniresult <- data.frame()
for (i in colnames(df[, 5:ncol(df)])) {
  unicox <- coxph(Surv(time = os_time, event = os_status) ~ df[, i], data = df)
  unisum <- summary(unicox)
  pvalue <- round(unisum$coefficients[, 5], 3)
  if (pvalue < pfilter) {
    uniresult <- rbind(uniresult, cbind(
      gene = i, HR = as.numeric(unisum$coefficients[, 2]),
      L95CI = as.numeric(unisum$conf.int[, 3]),
      H95CI = as.numeric(unisum$conf.int[, 4]),
      pvalue = as.numeric(unisum$coefficients[, 5])
    ))
  }
}

gene <- intersect(uniresult$gene, ptm_gene$gene)
pre_var <- gene  # <<< This is model_genes: all genes fed to plsRcox

cat("Number of model genes (pre_var):", length(pre_var), "\n")

# ---- 4. Prepare modeling data ----
unigene <- subset(df, select = c("os_status", "os_time", gene))
unigene <- unigene %>% tibble::rownames_to_column("ID")
colnames(unigene)[1:3] <- c("ID", "OS", "OS.time")

mm <- list(Training_Dataset = unigene)
mm <- lapply(mm, function(x) {
  x[, -c(1:3)] <- scale(x[, -c(1:3)])
  return(x)
})

est_data <- mm[["Training_Dataset"]]
est_dd   <- est_data[, c("OS.time", "OS", pre_var)]

# ---- 5. RSF feature selection (for reference only) ----
seed <- 1234
rf_nodesize <- 10
set.seed(seed)
rsf_fit <- rfsrc(Surv(OS.time, OS) ~ ., data = est_dd,
                 ntree = 1000, nodesize = rf_nodesize,
                 splitrule = "logrank", importance = TRUE,
                 proximity = TRUE, forest = TRUE, seed = seed)
rid <- var.select(rsf_fit)$topvars
cat("RSF selected features:", paste(rid, collapse = ", "), "\n")

# ---- 6. plsRcox model (trained on ALL pre_var) ----
ddist <- datadist(est_data)
options(datadist = "ddist")
cv_res <- cv.plsRcox(list(x = est_dd[, pre_var],
                           time = est_dd$OS.time,
                           status = est_dd$OS),
                     nt = 10, verbose = FALSE)
plsrcox_fit <- plsRcox(est_dd[, pre_var],
                        time = est_dd$OS.time,
                        event = est_dd$OS,
                        nt = as.numeric(cv_res[5]))

# ---- 7. Risk score & optimal cutoff ----
rs_vec <- as.numeric(predict(plsrcox_fit, type = "lp",
                             newdata = est_dd[, pre_var]))
train$rs <- rs_vec
res_cut <- surv_cutpoint(train, time = "os_time",
                         event = "os_status", variables = "rs")
sur_cutoff <- summary(res_cut)[1, 1]
cat("Optimal risk cutoff:", sur_cutoff, "\n")

# ---- 8. Cox nomogram model (rs + Age) ----
mmrf_sur <- openxlsx::read.xlsx("TCGA-MMRF/MMRF_survival_info.xlsx", sheet = 1)
mmrf_sur <- mmrf_sur %>% dplyr::select(2:6)
colnames(mmrf_sur) <- c("sample", "os_time", "os_status", "pfs_time", "pfs_status")
mmrf_sur$os_time  <- as.numeric(mmrf_sur$os_time) / 30
mmrf_sur$pfs_time <- as.numeric(mmrf_sur$pfs_time) / 30

mmrf_clinical <- openxlsx::read.xlsx("TCGA-MMRF/1_Individual patient features.xlsx",
                                      sheet = 5)
colnames(mmrf_clinical)[1] <- "sample"

train_clin <- mmrf_sur %>%
  filter(sample %in% train$sample) %>%
  inner_join(train[, c("sample", "rs")]) %>%
  inner_join(mmrf_clinical[, c(1, grep("Age", colnames(mmrf_clinical)))])

ddist2 <- datadist(train_clin)
options(datadist = "ddist2")
cox_mod <- cph(Surv(os_time, os_status) ~ rs + Age,
               surv = TRUE, x = TRUE, y = TRUE, data = train_clin)

# ---- 9. Save into package ----
highlight_genes <- c("FABP5", "BIRC5", "RASGRP3", "FGF13", "RAC3",
                      "BLNK", "MBL2", "CSF2", "SOD1")

# Scale params from the training expression matrix (ALL model genes)
train_sub <- as.data.frame(t(exp_train_full[pre_var, ]))
scale_params <- list(
  mean = colMeans(train_sub),
  sd   = apply(train_sub, 2, sd)
)

plsrcox_model <- plsrcox_fit
cox_model     <- cox_mod
risk_cutoff   <- sur_cutoff
model_genes   <- pre_var

pkg_model_dir <- "/Users/HXE/Desktop/MMimmuneRisk/inst/model"
dir.create(pkg_model_dir, recursive = TRUE, showWarnings = FALSE)
save(plsrcox_model, cox_model, scale_params,
     risk_cutoff, model_genes, highlight_genes,
     file = file.path(pkg_model_dir, "model_objects.rda"))

cat("\n========================================\n")
cat("Model objects saved successfully!\n")
cat("  Model genes:", length(model_genes), "\n")
cat("  Risk cutoff:", risk_cutoff, "\n")
cat("  Location:", file.path(pkg_model_dir, "model_objects.rda"), "\n")
cat("\nNext step — install the package:\n")
cat('  devtools::install("/Users/HXE/Desktop/MMimmuneRisk")\n')
cat("\nThen use:\n")
cat("  library(MMimmuneRisk)\n")
cat("  run_dynnom()           # Launch dynamic nomogram\n")
cat("  classify_mmrisk(expr)  # Batch classification\n")
cat("========================================\n")
