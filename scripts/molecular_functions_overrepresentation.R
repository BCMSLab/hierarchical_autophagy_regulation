# loading required libraries
library(tidyverse)
library(clusterProfiler)

# loading data
gene_annotations <- read_rds('data/gene_annotations.rds') 
gene_annotations <- split(gene_annotations$ENTREZID,
                          gene_annotations$category)
mf_annotations <- read_rds('data/mf_annotations.rds')
mf_annotations <- dplyr::select(mf_annotations, term, gene = gene_id)

binding_peaks <- read_rds('data/binding_peaks.rds')

peak_groups <- as_tibble(binding_peaks) %>%
    with(split(geneId, paste(factor, time, sep = '_'))) %>%
    map(unique)

# autophagy
autophagy_peak_groups <- map(peak_groups,
                             function(x) intersect(x, 
                                                   gene_annotations$Autophagy))

autophagy <- map_df(autophagy_peak_groups, function(x) {
    enricher(x,
             TERM2GENE = mf_annotations,
             pAdjustMethod = 'fdr',
             pvalueCutoff = 1,
             minGSSize = 1)@result
}, .id = 'group')

lipogenesis_peak_groups <- map(peak_groups,
                             function(x) intersect(x, 
                                                   gene_annotations$Lipogenesis))

lipogenesis <- map_df(lipogenesis_peak_groups, function(x) {
    enricher(x,
             TERM2GENE = mf_annotations,
             pAdjustMethod = 'fdr',
             pvalueCutoff = 1,
             minGSSize = 1)@result
}, .id = 'group')

list(
    'Autophagy' = autophagy,
    'Lipogenesis' = lipogenesis) %>%
    bind_rows(.id = 'category') %>%
    write_rds('data/molecular_functions_overrepresentation.rds')
