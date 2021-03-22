#!/bin/bash

# Define commands
RDATA=R CMD BATCH --vanilla $< logs/$(<F).Rout

all: data/gene_annotations.rds \
	data/mf_annotations.rds \
	data/model_hm.rds \
	data/model_factors_short.rds \
	data/model_factors_long.rds \
	data/differential_expression.rds \
	data/binding_peaks.rds \
	data/gene_set_enrichment.rds \
	data/peaks_overlap_enrichment.rds \
	data/targets_over_representation.rds \
	data/molecular_functions_overrepresentation.rds \
	manuscript.pdf
	
data/gene_annotations.rds: scripts/gene_annotations.R
	$(RDATA)
	
data/mf_annotations.rds: scripts/mf_annotations.R \
	data/gene_annotations.rds
	$(RDATA)
	
data/model_hm.rds: scripts/model_hm.R \
	$(bash find data/mm10/*) \
	$(bash find data/segmentation/hm/bins_200/*)
	$(RDATA)

data/model_factors_short.rds: scripts/model_factors_short.R \
	$(bash find data/mm10/*) \
	$(bash find data/segmentation/factors_short/bins_200/*)
	$(RDATA)

data/model_factors_long.rds: scripts/model_factors_long.R \
	$(bash find data/mm10/*) \
	$(bash find data/segmentation/factors_long/bins_200/*)
	$(RDATA)

data/differential_expression.rds: scripts/differential_expression.R \
	data/adipo_counts.rds \
	data/counts_cebpb_kd.rds \
	data/arrays_pparg_kd.rds
	$(RDATA)

data/binding_peaks.rds: scripts/binding_peaks.R \
	data/segmentation/input_files_factors_long.tsv \
    data/segmentation/input_files_factors_short.tsv \
    $(bash find data/peaks/*)
	$(RDATA)

data/gene_set_enrichment.rds: scripts/gene_set_enrichment.R \
	data/differential_expression.rds \
    data/gene_annotations.rds \
    data/binding_peaks.rds
	$(RDATA)
data/peaks_overlap_enrichment.rds: scripts/peaks_overlap_enrichment.R \
	data/binding_peaks.rds
	$(RDATA)

data/targets_over_representation.rds: scripts/targets_over_representation.R \
    data/gene_annotations.rds \
	data/binding_peaks.rds
	$(RDATA)

data/molecular_functions_overrepresentation.rds: scripts/molecular_functions_overrepresentation.R \
	data/gene_annotations.rds \
	data/mf_annotations.rds \
	data/binding_peaks.rds
	$(RDATA)

manuscript.pdf: manuscript.Rmd \
	$(bash find data/*)
	Rscript -e 'rmarkdown::render("manuscript.Rmd")' > logs/
	