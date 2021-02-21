[![R-CMD-check](https://github.com/kolesarm/ManyIV/workflows/R-CMD-check/badge.svg)](https://github.com/kolesarm/ManyIV/actions) [![Coverage status](https://codecov.io/gh/kolesarm/ManyIV/branch/master/graph/badge.svg)](https://codecov.io/github/kolesarm/ManyIV?branch=master)

# ManyIV

Inference in instrumental variables models with many instruments. The package
implements procedures from [Kolesár (2018)](https://doi.org/10.1016/j.jeconom.2018.01.004).

See the [package vignette](doc/ManyIV.pdf) for a description of the package
(available through `vignette("ManyIV")` once package is installed), and
the package [manual](doc/manual.pdf) for documentation of the package functions.

## Installation

You can get the current development version from GitHub:

``` r
if (!requireNamespace("remotes")) {
  install.packages("remotes")
}
remotes::install_github("kolesarm/ManyIV")
```
