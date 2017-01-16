.PHONY: clean data lint requirements

#################################################################################
# GLOBALS                                                                       #
#################################################################################

PROJECT_NAME = mri-rnaseq-integration
PYTHON = python3

#################################################################################
# COMMANDS                                                                      #
#################################################################################

all: requirements data models
.PHONY: all

data: requirements data/processed/gene-expression.nc data/processed/mri-features.nc
.PHONY: data

models: models/sfa.nc
.PHONY: models

# Environment

requirements:
	pip install -r requirements.txt
.PHONY: requirements

update-requirements:
	pip install -Ur requirements.txt
.PHONY: requirements

# Tools

clean:
	find . -name "*.pyc" -exec rm {} \;

lint:
	flake8 --exclude=docs/conf.py,venv .

# Help

help:
	@echo "Available make commands:"
	@echo "help                Show this help"
	@echo "clean               Delete all compiled Python files"
	@echo "requirements        Install Python Dependencies"
	@echo "update-requirements Update Python Dependencies"
	@echo "data                Make Dataset"
	@echo "lint                Lint using flake8"


########################################################################
# DATA                                                                 #
########################################################################

# Gene Expression (RNA-seq)
data/raw/gene-expression.nc:
	mkdir -p data/raw
	. ./.env; echo "user = $$BEEHUB_USERNAME:$$BEEHUB_PASSWORD" | curl -K - https://beehub.nl/home/tychobismeijer/Imagene/gene_expression/2016-03-09-gene-expression-imagene.nc -o $@

# Sample annotation
data/raw/sample-tracking.tsv:
	mkdir -p data/raw
	. ./.env; echo "user = $$BEEHUB_USERNAME:$$BEEHUB_PASSWORD" | curl -K - https://beehub.nl/home/tychobismeijer/Imagene/gene_expression/2016-06-14-sample-tracking.tsv -o $@

# MRI Features
data/raw/mri-features.xlsx:
	mkdir -p data/raw
	. ./.env; echo "user = $$BEEHUB_USERNAME:$$BEEHUB_PASSWORD" | curl -K - https://beehub.nl/home/tychobismeijer/Imagene/mri/2016-03-31-Tumor_Parenchym_Features_variablenamesupdated.xlsx -o $@

# Processed Data #

# MRI Features
data/processed/mri-features.nc: src/data/process_mri.py data/raw/mri-features.xlsx
	mkdir -p data/processed
	$(PYTHON) $^ $@

# Gene Expression
data/processed/gene-expression.nc: src/data/process_gene_expression.py data/raw/gene-expression.nc data/raw/sample-tracking.tsv data/annotation/ensembl_annotation.tsv
	mkdir -p data/processed
	$(PYTHON) $^ $@

# Annotation Data #

# Ensembl Annotation
data/annotation/ensembl_annotation.tsv: src/data/query_ensembl_reference.py
	$(PYTHON) $< $@

# External Data #

# Sparse-factor analysis results on TCGA Breast
data/external/tcga-breast-gexp+rppa+cn-sfa-solution.h5:
	. ./.env; echo "user = $$BEEHUB_USERNAME:$$BEEHUB_PASSWORD" | curl -K - https://beehub.nl/home/tychobismeijer/Imagene/external/tcga-breast-gexp+rppa+cn-sfa-solution.h5 -o $@

########################################################################
# MODELS                                                               #
########################################################################
models/sfa.nc : src/models/apply_sfa.py data/processed/gene-expression.nc data/external/tcga-breast-gexp+rppa+cn-sfa-solution.h5
	$(PYTHON) $^ $@
