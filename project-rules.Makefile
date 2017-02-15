PYTHON=python


########################################################################
# DATA                                                                 #
########################################################################

data: data/processed/gene-expression.nc data/processed/mri-features.nc
.PHONY: data


#---------
# Raw Data

# Gene Expression (RNA-seq)
data/raw/gene-expression.nc:
	mkdir -p data/raw
	. ./.env; echo "user = $$BEEHUB_USERNAME:$$BEEHUB_PASSWORD" | curl -K - https://beehub.nl/home/tychobismeijer/Imagene/gene_expression/2017-02-01-gene-expression-imagene.nc -o $@

# Sample annotation
data/raw/sample-tracking.tsv:
	mkdir -p data/raw
	. ./.env; echo "user = $$BEEHUB_USERNAME:$$BEEHUB_PASSWORD" | curl -K - https://beehub.nl/home/tychobismeijer/Imagene/gene_expression/2016-06-14-sample-tracking.tsv -o $@

# MRI Features
data/raw/mri-features.xlsx:
	mkdir -p data/raw
	. ./.env; echo "user = $$BEEHUB_USERNAME:$$BEEHUB_PASSWORD" | curl -K - https://beehub.nl/home/tychobismeijer/Imagene/mri/2016-03-31-Tumor_Parenchym_Features_variablenamesupdated.xlsx -o $@



#--------------
# External Data

# Sparse-factor analysis results on TCGA Breast
data/external/tcga-breast-gexp+rppa+cn-sfa-solution.h5:
	. ./.env; echo "user = $$BEEHUB_USERNAME:$$BEEHUB_PASSWORD" | curl -K - https://beehub.nl/home/tychobismeijer/Imagene/external/tcga-breast-gexp+rppa+cn-sfa-solution.h5 -o $@

# Ensembl Annotation
data/external/ensembl_annotation.tsv: src/data/query_ensembl_reference.py
	$(PYTHON) $< $@

#---------------
# Processed Data

# MRI Features
data/processed/mri-features.nc: src/data/process_mri.py data/raw/mri-features.xlsx
	mkdir -p data/processed
	$(PYTHON) $^ $@

# Gene Expression, add annotation and log2 CPM
data/processed/gene-expression.nc: src/data/process_gene_expression.py data/raw/gene-expression.nc data/raw/sample-tracking.tsv data/external/ensembl_annotation.tsv
	mkdir -p data/interim
	$(PYTHON) $^ $@


########################################################################
# MODELS                                                               #
########################################################################

models: models/sfa.nc
.PHONY: models

models/sfa.nc : src/models/apply_sfa.py data/processed/gene-expression.nc data/external/tcga-breast-gexp+rppa+cn-sfa-solution.h5
	$(PYTHON) $^ $@
