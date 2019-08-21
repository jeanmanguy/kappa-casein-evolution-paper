
<!-- README.md is generated from README.Rmd. Please edit that file -->

# kappa-casein-evolution-paper

![](https://img.shields.io/badge/R%20≥-3.5-green.svg)
![](https://img.shields.io/badge/Python%20≥-3-green.svg)
[![](https://img.shields.io/badge/doi-10.5281/zenodo.2587122-blue.svg)](https://doi.org/10.5281/zenodo.2587122)
<!-- [![](https://img.shields.io/badge/doi--red.svg)](https://doi.org/) -->

This repository contains the data and code for our paper:

> Jean Manguy and Denis C. Shields, (2019). *Implications of
> kappa-casein evolutionary diversity for the self-assembly and
> stabilization of casein micelles*. in preparation \<\>

This repository is archived on Zenodo:
[DOI: 10.5281/zenodo.2587122](https://doi.org/10.5281/zenodo.2587122)

## Analysis

The code for the figures and the output can be viewed in a web browser,
download and open the file `figures_and_tables_paper.html`.

We use
[{ProjectTemplate}](https://github.com/KentonWhite/ProjectTemplate) to
manage the analysis. The code was tested only on Linux.

### Install dependencies

1.  download IUPred2A from <https://iupred2a.elte.hu/download> and unzip
    the archive in `lib/iupred2a/` (non-commercial users only)
2.  install R \>= 3.5 and python \>= 3
3.  install R packages (see `session_info.txt`)

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

Run in a unix shell

``` sh
R -e "rmarkdown::render('figures_and_tables_paper.Rmd')"
```
