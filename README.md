# Kappa-casein evolutionary diversity for the self-assembly and stabilization of casein micelles


This repository contains the code and data for the paper: "Kappa-casein evolutionary diversity for the self-assembly and stabilization of casein micelles".

## To run the analysis

We use ProjectTemplate to manage the analysis. The code was tested only on Linux.

### Install dependencies

1. IUPred2A
2. R packages (see `packages.txt`

### Generates the figures only

1. go to root of this repository

2. run:

```
library(ProjectTemplate)
load.project()
```
### Run the rmarkdown file

#### With Rstudio

1. open this repository with Rstudio
2. open `analysis.Rmd` from Rstudio
3. Click on `Run` > `run all` (or Ctrl+Alt+R)

#### Without Rstudio

Run

```
R -e "rmarkdown::render('analysis.Rmd')"
```


## Citation

```
Paper in preparation
```
