FROM rocker/verse

RUN Rscript -e "install.packages('reshape2')"
RUN Rscript -e "install.packages('scales')"
RUN Rscript -e "install.packages('xtable')"
RUN Rscript -e "install.packages('circlize')"
RUN Rscript -e "install.packages('ggupset')"
RUN Rscript -e "install.packages('cowplot')"
RUN Rscript -e "install.packages('MASS')"
RUN Rscript -e "install.packages('boot')"
RUN Rscript -e "install.packages('survival')"
RUN Rscript -e "install.packages('cluster')"
RUN Rscript -e "install.packages('BiocManager')"

RUN Rscript -e "BiocManager::install('ChIPseeker')"
RUN Rscript -e "BiocManager::install('ComplexHeatmap')"
RUN Rscript -e "BiocManager::install('DESeq2')"
RUN Rscript -e "BiocManager::install('limma')"
RUN Rscript -e "BiocManager::install('fgsea')"
RUN Rscript -e "BiocManager::install('clusterProfiler')"
RUN Rscript -e "BiocManager::install('GO.db')"
RUN Rscript -e "BiocManager::install('TxDb.Mmusculus.UCSC.mm10.knownGene')"
RUN Rscript -e "BiocManager::install('org.Mm.eg.db')"
RUN Rscript -e "BiocManager::install('SummarizedExperiment')"
RUN Rscript -e "BiocManager::install('GenomicRanges')"

ADD ./segmentr /home/segmentr
RUN ls home
RUN Rscript -e "devtools::install('/home/segmentr')"
