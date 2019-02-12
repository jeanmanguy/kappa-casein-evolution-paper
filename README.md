
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Implications of kappa-casein evolutionary diversity for the self-assembly and stabilization of casein micelles

![](https://img.shields.io/badge/R%20%3E=-3.4.4-green.svg)
![](https://img.shields.io/badge/Python%20%3E=-3-green.svg)

This repository contains the code and data for the paper: “Implications
of kappa-casein evolutionary diversity for the self-assembly and
stabilization of casein micelles”.

The code for the figures and the output can be viewed in a web browser,
download and open the file `figures_and_tables_paper.html`.

## To run the analysis

We use ProjectTemplate to manage the analysis. The code was tested only
on Linux.

### Install dependencies

1.  download IUPred2A from <https://iupred2a.elte.hu/download> and unzip
    the archive in `lib/iupred2a/`
2.  install R packages (see `packages.txt`)

### Generates the figures only

1.  go to root of this repository

2.  run:

<!-- end list -->

``` r
library(ProjectTemplate)
load.project()
```

### Run the rmarkdown file

#### With Rstudio

1.  open this repository with Rstudio
2.  open `figures_and_tables_paper.Rmd` from Rstudio
3.  Click on `Run` \> `run all` (or Ctrl+Alt+R)

#### Without Rstudio

Run

``` sh
R -e "rmarkdown::render('figures_and_tables_paper.Rmd')"
```

## Citation

    Paper in preparation
