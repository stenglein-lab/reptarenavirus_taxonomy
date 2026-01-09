# A proposal for a simplified reptarenavirus taxonomy based on reassortment compatibility

This repository contains the code and data used to create the paper titled "A proposal for a simplified reptarenavirus taxonomy based on reassortment compatibility". 

### Reproducible implementation

This [paper](./paper/paper.docx) ([paper in qmd format](./paper/paper.qmd)]is implemented as a reproducible workflow that uses nextflow, singularity, and quarto markdown.

This approach was inspired by papers and workflows used by the Schloss lab at the University of Michigan, such as [this paper](https://github.com/SchlossLab/Lesniak_restoreCR_mBio_2022/tree/main/submission).

### How to reproduce this paper

```
# clone the repository
git clone https://github.com/stenglein-lab/reptarenavirus_taxonomy.git

# cd into the repository directory 
cd reptarenavirus_taxonomy

# run main entry point bash script to reproduce paper
./reproduce_paper.sh
```



