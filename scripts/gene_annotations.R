# loading required libraries
library(org.Mm.eg.db)
library(tidyverse)

# make a list of terms of interest
terms <- list(
    'Autophagy' = 'GO:0006914',
    'Autophagy Regulation' = c('GO:0010507', 'GO:0010508'),
    'Autophagy Subtypes' = c('GO:0000422', 'GO:0030242', 'GO:0061684',
                             'GO:0061738', 'GO:0016236'),
    'Selective Autophagy' = c('GO:0035973', 'GO:0000423',
                              'GO:0061709', 'GO:0098792'),
    'Lipogenesis' = 'GO:0006629',
    'Lipogenesis Regulation' = c('GO:0045833', 'GO:0045834')
)

# get autophagy subsets
gene_annotations <- map_df(
    terms,
    function(x) {
        AnnotationDbi::select(org.Mm.eg.db,
            x,
            c('SYMBOL', 'ENTREZID'),
            'GO'
        )
    }
, .id = 'category')

write_rds(gene_annotations, 'data/gene_annotations.rds')

