# loading libraries
library(tidyverse)
library(reshape2)
library(target)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

# load data
## prepare regions
### gene annotations
gene_annotations <- AnnotationDbi::select(org.Mm.eg.db,
                                          keys = keys(org.Mm.eg.db, 'SYMBOL'),
                                          'ENTREZID',
                                          'SYMBOL') %>%
    setNames(c('symbol', 'gene_id')) %>%
    na.omit() %>%
    as_tibble()

### differential expression
differential_expression <- read_rds('data/differential_expression.rds') %>%
    na.omit() %>%
    filter(course != 'None') %>%
    dcast(gene_id ~ course + time, value.var = 'fold_change') %>%
    dplyr::select(symbol = gene_id, everything()) %>%
    left_join(gene_annotations) %>%
    as_tibble()

### gene coordinates
gene_coordinates <- genes(TxDb.Mmusculus.UCSC.mm10.knownGene,
                          filter = list(gene_id = gene_annotations$gene_id))
ind <- match(gene_coordinates$gene_id, differential_expression$gene_id)

### make regions granges
regions <- granges(gene_coordinates)
mcols(regions) <- differential_expression[ind,]

## prepare binding peaks
bindin_peaks <- read_rds('data/binding_peaks.rds')
group <- paste(bindin_peaks$factor, bindin_peaks$time, sep = '_')
bindin_peaks <- split(bindin_peaks, group)

# apply target
pp <-names(bindin_peaks)[grepl('CEBPB|PPARG', names(bindin_peaks))]
rr <- names(mcols(regions))

targets <- intersect(pp, rr)
names(targets) <- targets

dt <- map_df(targets, 
       function(x) {
           direct_targets(bindin_peaks[[x]],
                          regions,
                          'symbol',
                          x) %>%
               as_tibble() %>% 
               dplyr::select(symbol,
                             starts_with('score'),
                             starts_with('stat'),
                             starts_with('rank'))
       }, .id = 'group') %>%
    separate(group, c('factor', 'time'))

write_rds(dt, 'data/direct_targets.rds')
