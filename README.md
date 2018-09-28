# OneKP-gene-family-evo
In this repo you find the scripts that produced the results and figures of the gene family evolution analysis from the OneKP capstone paper (Ref will follow), as well as the validation of the approach.

## Getting started
First, you can install all programms and packages needed via conda ([installation instructions](https://conda.io/docs/user-guide/install/index.html) for conda).
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
