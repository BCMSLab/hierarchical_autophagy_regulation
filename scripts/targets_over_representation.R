# loading required libraries
library(tidyverse)
library(clusterProfiler)

# loading data
gene_annotations <- read_rds('data/gene_annotations.rds') 
gene_annotations <- split(gene_annotations$ENTREZID,
                          gene_annotations$category)

binding_peaks <- read_rds('data/binding_peaks.rds')
   
peak_groups <- as_tibble(binding_peaks) %>%
    dplyr::select(factor, time, gene = geneId) %>%
    unite(term, c('factor', 'time')) %>%
    unique()

# enrichment all
peak_list <- split(peak_groups$gene, peak_groups$term)

# all <- map_df(peak_list,function(x) {
#     enricher(x,
#              TERM2GENE = peak_groups,
#              pAdjustMethod = 'fdr',
#              minGSSize = 5,
#              pvalueCutoff = 1)@result
# }, .id = 'group')

# autophagy
autophagy_peak_groups <- peak_groups %>%
    filter(gene %in% gene_annotations$Autophagy)

autophagy_peak_list <- split(autophagy_peak_groups$gene,
                             autophagy_peak_groups$term)

autophagy <- map_df(autophagy_peak_list, function(x) {
    enricher(x,
             TERM2GENE = autophagy_peak_groups,
             pAdjustMethod = 'fdr',
             pvalueCutoff = 1,
             minGSSize = 1)@result
}, .id = 'group')

# lipogenesis
lipogenesis_peak_groups <- peak_groups %>%
    filter(gene %in% gene_annotations$Lipogenesis)
lipogenesis_peak_list <- split(lipogenesis_peak_groups$gene, lipogenesis_peak_groups$term)

lipogenesis <- map_df(lipogenesis_peak_list, function(x) {
    enricher(x,
             TERM2GENE = lipogenesis_peak_groups,
             pAdjustMethod = 'fdr',
             pvalueCutoff = 1,
             minGSSize = 1)@result
}, .id = 'group')

# write to disk
list(
    # 'All' = all,
    'Autophagy' = autophagy,
    'Lipogenesis' = lipogenesis) %>%
    bind_rows(.id = 'category') %>%
    write_rds('data/targets_over_representation.rds')
