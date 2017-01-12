.PHONY: clean data lint requirements

#################################################################################
# GLOBALS                                                                       #
#################################################################################

PROJECT_NAME = mri-rnaseq-integration
PYTHON = python3

#################################################################################
# COMMANDS                                                                      #
#################################################################################

all: requirements data
.PHONY: all

data: requirements
.PHONY: data

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


#################################################################################
# PROJECT RULES                                                                 #
#################################################################################

