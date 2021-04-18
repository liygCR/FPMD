# FPMD
`FPMD` is an `R` package for the paper 'On Functional Processes with Multiple Discontinuities'.


## Installation
```r
library(devtools)
devtools::install_github("liygCR/FPMD/FPMD")
```

## Description
This package handles the problem of estimating multiple discontinuities for a functional data process. A half-kernel approach is considered that addresses the inference of the total number, locations, and jump sizes of the changes in the mean function.

The `simulation.R` refers to the simulation settings in the main paper. 
The `Meanfromfdapace.R` refers to the implementations in Zhang and Wang (2016, AOS).   

## Citation
> Li, J., Li, Y., and Hsing, T. (2021). On Functional Processes with Multiple Discontinuities. (Submitted)
