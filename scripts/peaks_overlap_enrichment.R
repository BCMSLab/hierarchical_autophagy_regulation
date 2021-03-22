# loading libraries
library(tidyverse)

# loading binding sites
binding_peaks <- read_rds('data/binding_peaks.rds')

dir.create('data/binding_peaks')
dir.create('data/binding_peaks_segmentation')

binding_peaks %>%
    as_tibble() %>%
    unite(group, c('factor', 'time'), remove = FALSE) %>%
    group_split(group) %>%
    map(function(x) {
        nm <- unique(x$group)
        dplyr::select(x, seqnames, start, end) %>%
            write_tsv(paste0('data/binding_peaks/', nm, '.bed'))
    })

segments_files <- list.files('data/segmentation/hm', 
           pattern = '*_segments.bed',
           full.names = TRUE)
names(segments_files) <- str_split(segments_files, '/|_9', simplify = TRUE)[, 4]

cmd <- paste(segments_files,
             'data/binding_peaks',
             paste0('data/binding_peaks_segmentation/', names(segments_files)))

map(cmd, function(x) {
    system2(
        command = "java",
        args = paste(
            "-mx12000M -jar",
            'ChromHMM/ChromHMM.jar',
            'OverlapEnrichment',
            '-noimage',
            x))
})
nms <- str_split(list.files('data/binding_peaks/'), '\\.', simplify = TRUE)[ ,1]

peaks_overlap_enrichment <- list.files('data/binding_peaks_segmentation',
           full.names = TRUE) %>%
    set_names(c(0, 168, 4, 48)) %>%
    map_df(read_tsv,
           col_names = c('state', 'genome', nms),
           skip = 1,
           n_max = 9,
           .id = 'hour')

write_rds(peaks_overlap_enrichment,
          'data/peaks_overlap_enrichment.rds')
