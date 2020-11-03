# load required libraries
library(segmentr)

# load model
factors_short <- learn_model(
    inputdir = 'data/segmentation/factors_short/bins_100/',
    output = 'data/segmentation/factors_short/',
    coordsdir = 'data/mm10/COORDS',
    anchorsdir = 'data/mm10/ANCHORFILES',
    chromsizefile = 'data/mm10/mm10.txt',
    assembly = 'mm10',
    numstates = 10,
    binsize = 100,
    cells = paste('hour', c(0, 4), sep = '_'),
    annotation = 'autophagy',
    read_only = TRUE)

# get txdb object
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene

# annotate segments
slot(factors_short, 'segment') <- as.list(
    annotate_segments(
        segment(factors_short),
        TxDb = txdb,
        level = 'gene',
        assignGenomicAnnotation = FALSE,
        verbose = FALSE
    )
)

# write to object
readr::write_rds(factors_short, 'data/model_factors_short.rds')
