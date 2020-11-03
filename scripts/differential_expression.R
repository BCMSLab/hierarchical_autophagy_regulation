# loading required libraries
library(tidyverse)
library(SummarizedExperiment)
library(DESeq2)
library(limma)

# timecourse, no perturbations
# loading data
adipo_counts <- read_rds('data/adipo_counts.rds')
adipo_counts$time <- relevel(factor(adipo_counts$time), ref = '0')

# may be filter low reads
# may be use combat for batch effects removal

dds <- DESeqDataSet(adipo_counts, ~time)
dds <- DESeq(dds)
None <- map(resultsNames(dds)[c(5, 6, 9, 11, 14:17)],
              function(x) {
                  results(dds,
                          name = x,
                          pAdjustMethod = 'fdr',
                          tidy = TRUE) %>%
                      mutate(time = x)
              }) %>%
    bind_rows() %>%
    select(time,
           gene_id = row,
           fold_change = log2FoldChange,
           pvalue,
           padj) %>%
    as_tibble() %>%
    mutate(time = as.numeric(str_split(time, '_', simplify = TRUE)[, 2]))

# timecourse, with perturbations
se <- read_rds('data/counts_cebpb_kd.rds')
se <- se[rowSums(assay(se)) > 0,]

CEBPB <- map(c('0' = 0, '4' = 4), function(x) {
    e <- se[,se$time == x]
    dds <- DESeqDataSet(e, design = ~group-1)
    dds <- DESeq(dds)
    
    results(dds,
            tidy = TRUE,
            pAdjustMethod = 'fdr',
            cooksCutoff=FALSE)
}) %>%
    bind_rows(.id = 'time') %>%
    select(time,
           gene_id = row,
           fold_change = log2FoldChange,
           pvalue,
           padj) %>%
    mutate(time = as.numeric(time)) %>%
    as_tibble()

# pparg 'GSE12929'
eset <- read_rds('data/arrays_pparg_kd.rds')
eset <- eset[, !is.na(eset$time)]
exprs(eset) <- normalizeBetweenArrays(exprs(eset))

exprs(eset) <- log2(exprs(eset) + 1)

PPARG <- map(c('0' = 0, '2' = 2, '4' = 4, '5' = 5, '6' = 6), function(x) {
    e <- eset[, eset$time == x]
    mat <- exprs(e)
    mod <- model.matrix(~group, data = pData(e))
    fit <- lmFit(mat, design = mod, method = 'robust')
    fit <- eBayes(fit)
    
    se <- sqrt(fit$s2.post) * fit$stdev.unscaled
    
    topTable(fit,
             sort.by = 'none',
             number = Inf,
             adjust.method = 'fdr') %>%
        rownames_to_column('gene_id') %>%
        mutate(se = se[, 2]) %>%
        as_tibble()
}) %>%
    bind_rows(.id = 'time') %>%
    select(time,
           gene_id,
           fold_change = logFC,
           pvalue = P.Value,
           padj = adj.P.Val) %>%
    mutate(time = as.numeric(time) * 24) %>%
    as_tibble()

list(None = None,
     CEBPB = CEBPB,
     PPARG = PPARG) %>%
    bind_rows(.id = 'course') %>%
    mutate(course = factor(course, levels = c('None', 'PPARG', 'CEBPB'))) %>%
    write_rds('data/differential_expression.rds')
