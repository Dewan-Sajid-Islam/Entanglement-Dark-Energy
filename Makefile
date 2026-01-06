.PHONY: all install test clean figures paper help

# Variables
PYTHON = python
PIP = pip
NOTEBOOK = jupyter notebook

# Default target
all: install test figures

# Help
help:
	@echo "Available commands:"
	@echo "  make install     Install dependencies"
	@echo "  make test        Run tests"
	@echo "  make figures     Generate all figures"
	@echo "  make paper       Compile LaTeX paper"
	@echo "  make clean       Clean generated files"
	@echo "  make all         Install, test, and generate figures"

# Installation
install:
	$(PIP) install -r requirements.txt

# Testing
test:
	$(PYTHON) -m pytest tests/ -v

# Generate figures
figures:
	$(PYTHON) scripts/generate_figures.py

# Compile paper
paper:
	cd latex && \
	pdflatex main.tex && \
	bibtex main && \
	pdflatex main.tex && \
	pdflatex main.tex

# Clean up
clean:
	rm -rf __pycache__/
	rm -rf scripts/__pycache__/
	rm -rf tests/__pycache__/
	rm -f latex/*.aux latex/*.log latex/*.out latex/*.toc
	rm -f latex/*.synctex.gz latex/*.fls latex/*.fdb_latexmk
	rm -f figures/*.pdf figures/*.png
	rm -f data/chains/*.h5 data/chains/*.npy
	find . -name "*.pyc" -delete
	find . -name ".DS_Store" -delete

# Quick tests
quick-test:
	$(PYTHON) -c "from scripts.background_evolution import EDECosmology; ede = EDECosmology(); print('✓ EDE model loaded')"

# Open notebooks
notebooks:
	$(NOTEBOOK) notebooks/

# Check environment
check:
	$(PYTHON) --version
	$(PIP) list | grep -E "numpy|scipy|matplotlib"

# Update requirements
update-reqs:
	$(PIP) freeze > requirements.txt
