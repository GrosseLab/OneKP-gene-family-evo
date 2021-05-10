# OneKP-gene-family-evo
In this repo you find the scripts that produced the results and figures of the gene family evolution analysis from the [OneKP capstone paper](https://www.nature.com/articles/s41586-019-1693-2), as well as the validation of the approach.

## Bofore cloning
This repo contains large files that are stored with [git-lfs](https://git-lfs.github.com). You need to install git-lfs before cloning so the large files will automatically be downloaded. Please also see [this Stack Overflow question](https://stackoverflow.com/questions/56283415/how-do-i-clone-a-repository-that-includes-git-lfs-files).

## Getting started
First, you can install all programms and packages needed via conda ([installation instructions](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) for conda).
```
conda env create -f conda_env.yml
source activate onekp-gf-evo
```
Then you can run the main analysis from the paper by starting R and running
```
source("analysis.r")
```
To run the two validations uncomment the respective studies in the `analysis.r` script and rerun th above command. Finally, run the two other r-scripts
```
source("angiosperm_bias.r")
source("validation_on_model_organisms.R")
```
