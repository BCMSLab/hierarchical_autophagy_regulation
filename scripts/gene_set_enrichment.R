# loading required libraries
library(tidyverse)
library(fgsea)

# loading data
diff_exprs <- read_rds('data/differential_expression.rds')

# get fold change for all genes
stats <- diff_exprs %>%
    na.omit() %>%
    group_by(course, time) %>%
    group_split() %>%
    map(function(x) {
        vec <- x$fold_change
        names(vec) <- x$gene_id
        vec
    })

names(stats) <- unite(diff_exprs, 'group', course, time) %>%
    pull(group) %>%
    unique()

# calculate enrichment
## annotations
gene_annotations <- read_rds('data/gene_annotations.rds')
gene_annotations %>% head()

processes_sets <- split(gene_annotations$SYMBOL, gene_annotations$GO)
processes_sets <- map(processes_sets, unique)

set.seed(12345)
processes_enrichment <- map_df(stats,
       function(x) {
           fgsea::fgsea(processes_sets,
                        x,
                        nperm = 1000,
                        nproc = 4)
       }, .id = 'group') %>%
    as_tibble()
## binding peaks
binding_peaks <- read_rds('data/binding_peaks.rds')

peaks_sets <- binding_peaks %>%
    as_tibble() %>%
    na.omit() %>%
    dplyr::select(factor, time, geneId) %>%
    unique() %>%
    inner_join(tibble(geneId = gene_annotations$ENTREZID, 
                      gene_id = gene_annotations$SYMBOL)) %>%
    unite('group', c('factor', 'time')) %>%
    with(split(gene_id, group)) %>%
    map(unique)

set.seed(12345)
peaks_set_enrichment <- map_df(stats,
                               function(x) {
                                   fgsea::fgsea(peaks_sets,
                                                x,
                                                nperm = 1000,
                                                nproc = 4)
                               }, .id = 'group') %>%
    as_tibble()

# autophagy sets
autophagy_peaks <- map(peaks_sets, function(x) {
    intersect(x, gene_annotations$SYMBOL[gene_annotations$category == 'Autophagy'])
})

set.seed(12345)
autophagy_peaks_enrichment <- map_df(stats,
                               function(x) {
                                   fgsea::fgsea(autophagy_peaks,
                                                x,
                                                nperm = 1000,
                                                nproc = 4)
                               }, .id = 'group') %>%
    as_tibble()

# lipogenesis sets
lipogenesis_peaks <- map(peaks_sets, function(x) {
    intersect(x, gene_annotations$SYMBOL[gene_annotations$category == 'Lipogenesis'])
})

set.seed(12345)
lipogenesis_peaks_enrichment <- map_df(stats,
                                     function(x) {
                                         fgsea::fgsea(lipogenesis_peaks,
                                                      x,
                                                      nperm = 1000,
                                                      nproc = 4)
                                     }, .id = 'group') %>%
    as_tibble()

# write to disk
list('GO Biological Processes' = processes_enrichment,
     'Peak Sets' = peaks_set_enrichment,
     'Autophagy Peaks' = autophagy_peaks_enrichment,
     'Lipogenesis Peaks' = lipogenesis_peaks_enrichment) %>%
    bind_rows(.id = 'categories') %>%
    write_rds('data/gene_set_enrichment.rds')
