# loading required libraries
library(org.Mm.eg.db)
library(tidyverse)
library(GO.db)

# make a list of terms of interest
gene_annotations <- read_rds('data/gene_annotations.rds')
gene_annotations <- with(gene_annotations, split(ENTREZID, category))
gene_annotations <- gene_annotations[c('Autophagy', 'Lipogenesis')]

# get autophagy subsets
mf <- map_df(
    gene_annotations,
    function(x) {
        AnnotationDbi::select(org.Mm.eg.db,
                              x,
                              'GOALL',
                              'ENTREZID'
        )
    }
    , .id = 'category') %>%
    filter(ONTOLOGYALL == 'MF') %>%
    group_by(category, GOALL) %>%
    mutate(n = length(unique(ENTREZID))) %>%
    filter(n > 3) %>%
    mutate(term = unlist(map(GOALL, ~GOTERM[[.x]]@Term))) %>%
    as_tibble() %>%
    dplyr::select(category, term, gene_id = ENTREZID)

write_rds(mf, 'data/mf_annotations.rds')
