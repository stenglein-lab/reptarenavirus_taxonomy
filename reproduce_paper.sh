#!/bin/bash -x

# this script will reproduce the paper: 
#
# "A proposal for a simplified reptarenavirus taxonomy based on reassortment compatibility"
#
# it does this by running several nextflow workflows, several R scripts.  
# These workflows create data that is in turn turned into figures and tables by the R scripts.
# The figures, tables, and some R-calculated values are read in to the quarto markdown paper (paper/paper.qmd)
# The paper is finally rendered to HTML and .docx using quarto
#
#
# This is implemented in this shell script and I'm not wild about this strategy.  
# For one thing, this shell script will just re-run everything every time, which is 
# clearly not ideal.  It also does not handle dependencies using best practice tools
# like singularity or conda. 
#
# Alternative approaches could include:
# - putting all of this into a nextflow workflow
# - using makefiles like Pat Schloss does: 
#    e.g.: https://github.com/SchlossLab/Lesniak_restoreCR_mBio_2022/blob/main/Makefile
# - using snakefiles like Pat Schloss does: 
#    e.g.: https://github.com/SchlossLab/Schloss_Rarefaction_mSphere_2024/blob/main/Snakefile
#
#
#  Dependencies (version used for paper): 
#   - nextflow (v25.10.0)
#   - singularity (ce version 4.2.1)
#   - quarto (v1.8.26)
#   - R (4.5.2)and a number of R packages, loaded in the R scripts in the bin/ subdirectory
# 
#  The requirement of all these dependencies makes this not as reproducible as it could be. 
#  However I am putting this out there following the Schloss Suck until you don't principle.
#
# Mark Stenglein 11/2025

# run nextflow workflow to extract metadata re: reptarenavirus sequences from genbank files
cd metadata
./run_extract
cd ..

# run nextflow workflow to create alignments and trees
cd make_alignments_and_trees
./run_align
cd ..

# run R script to use Entrez E-Utils to fetch NCBI accessions linked to particular reptarenavirus taxids
cd bin
Rscript fetch_reptarenavirus_taxa_accessions.R
cd ..

# run R script to generate figures, data values, and tables 
cd bin
Rscript reptarenavirus_taxonomy.R
cd ..

# render paper using quarto
cd paper
quarto render paper.qmd --to docx
cd ..
