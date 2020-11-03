# load required libraries
library(segmentr)

# load model
model_hm <- learn_model(
    inputdir = 'data/segmentation/hm/bins_200',
    output = 'data/segmentation/hm',
    coordsdir = 'data/mm10/COORDS',
    anchorsdir = 'data/mm10/ANCHORFILES',
    chromsizefile = 'data/mm10/mm10.txt',
    assembly = 'mm10',
    numstates = 9,
    binsize = 200,
    cells = paste('hour', c(0, 48, 168), sep = '_'),
    annotation = 'autophagy',
    read_only = TRUE)

# get txdb object
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene::TxDb.Mmusculus.UCSC.mm10.knownGene

# annotate segments
slot(model_hm, 'segment') <- as.list(
    annotate_segments(
        segment(model_hm),
        TxDb = txdb,
        level = 'gene',
        assignGenomicAnnotation = FALSE,
        verbose = FALSE
    )
)

# write to object
readr::write_rds(model_hm, 'data/model_hm.rds')
