#################################################################################
# GLOBALS                                                                       #
#################################################################################

PYTHON = python3

#################################################################################
# COMMANDS                                                                      #
#################################################################################

all: venv requirements data models
.PHONY: all

# Environment

venv:
	$(PYTHON) -m venv venv

ACTIVATE_ENV=. venv/bin/activate

requirements: venv
	$(ACTIVATE_ENV); pip install -q -r requirements.txt
.PHONY: requirements

update-requirements: venv
	$(ACTIVATE_ENV); pip install -Ur requirements.txt
.PHONY: update-requirements

# Tools

clean:
	find . -name "*.pyc" -exec rm {} \;
.PHONY: clean

lint:
	$(ACTIVATE_ENV); flake8 --exclude=docs/conf.py,venv .
.PHONY: lint

# Help

help:
	@echo "Available make commands:"
	@echo "help                Show this help"
	@echo "clean               Delete all compiled Python files"
	@echo "venv                Make virtual environment"
	@echo "requirements        Install Python dependencies"
	@echo "update-requirements Update Python dependencies to their latest version"
	@echo "data                Make all datasets"
	@echo "models              Make all models"
	@echo "lint                Lint using flake8"

# Load environment and make files with project rules.

Makefile: ;
%: requirements
	$(ACTIVATE_ENV); $(MAKE) -f project-rules.Makefile $@
