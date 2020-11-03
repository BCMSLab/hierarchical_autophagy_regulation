# loading required libraries
library(tidyverse)
library(GenomicRanges)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ChIPseeker)

# loading metadata
md <- map_df(
    c('data/segmentation/input_files_factors_long.tsv',
      'data/segmentation/input_files_factors_short.tsv'),
    read_tsv,
    col_names = c('time', 'factor', 'gsm')
) %>%
    mutate(time = str_split(time, '_', simplify = TRUE)[, 2],
           gsm = str_split(gsm, '\\.', simplify = TRUE)[, 1],
           file = file.path('data/peaks', paste0(gsm, '_peaks.xls'))) %>%
    filter(file.exists(file))

fls <- md$file
names(fls) <- md$gsm

# reading peaks
peaks <- map(fls, read_tsv, comment = '#')
peaks <- bind_rows(peaks, .id = 'gsm')
peaks <- left_join(peaks, md)

gr <- makeGRangesFromDataFrame(peaks,
                               keep.extra.columns = TRUE)

# annotate peaks
gr <- annotatePeak(gr,
                   TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene,
                   level = 'gene')
gr <- as.GRanges(gr)
gr <- gr[abs(gr$distanceToTSS) < 50000]

write_rds(gr, 'data/binding_peaks.rds')
