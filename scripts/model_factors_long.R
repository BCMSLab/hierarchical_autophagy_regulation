# load required libraries
library(segmentr)

# load model
factors_long <- learn_model(
    inputdir = 'data/segmentation/factors_long/bins_100/',
    output = 'data/segmentation/factors_long/',
    coordsdir = 'data/mm10/COORDS',
    anchorsdir = 'data/mm10/ANCHORFILES',
    chromsizefile = 'data/mm10/mm10.txt',
    assembly = 'mm10',
    numstates = 8,
    binsize = 100,
    cells = paste('hour', c(0, 24, 48, 72, 96, 144, 168), sep = '_'),
    annotation = 'autophagy',
    read_only = TRUE)

# get txdb object
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene

# annotate segments
slot(factors_long, 'segment') <- as.list(
    annotate_segments(
        segment(factors_long),
        TxDb = txdb,
        level = 'gene',
        assignGenomicAnnotation = FALSE,
        verbose = FALSE
    )
)

# write to object
readr::write_rds(factors_long, 'data/model_factors_long.rds')
