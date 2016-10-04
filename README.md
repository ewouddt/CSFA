
<!-- README.md is generated from README.Rmd. Please edit that file -->

Installing the CSFA Package
---------------------------


``` r
setRepositories(ind=c(1:5))
install.packages("CSFA",repos="http://R-Forge.R-project.org")
```

or

``` r
setRepositories(ind=c(1:5))
install.packages("devtools") # If not yet installed on your R Version
devtools::install_github("hadley/devtools") # Only run this if your currently installed 
                                            # devtools version is <= 1.12 (recursive dependencies bug)

devtools::install_github("ewouddt/CSFA")
```
