#!/bin/bash -x

# this script will reproduce the paper: 
#
# "Pervasive coinfection, reassortment, and recombination inspire a new *Reptarenavirus* taxonomy"
#
# it does this by running several nextflow workflows, several R scripts.  
# These workflows create data that is in turn turned into figures and tables by the R scripts.
# The figures, tables, and some R-calculated values are read in to the quarto markdown paper (paper/paper.qmd)
# The paper is finally rendered to HTML and .docx using quarto
#
#  Dependencies (version used for paper): 
#   - nextflow (v25.10.0)
#   - singularity (ce version 4.2.1)
# 
# Mark Stenglein 10/2026

# The nextflow pipeline in main.nf reproduces analyses and renders the paper

nextflow run -resume main.nf -output-dir workflow_output

