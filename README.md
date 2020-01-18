FrenchFISH
=====
`FrenchFISH` is an R package for for correcting spot counts used for quantifying DNA copy-number from fluorescence in situ hybridisation of tissue sections. It can correct either spot counts from manual counting, or from an automatic spot counting software. The underlying models used to correct these counts are Poissonian. Details can be found in this [manuscript](https://doi.org/10.1101/487926).

Installation  
----
```r
require(devtools)
devtools::install_github("markowetzlab/frenchFISH")
```

Getting Started
----
```
vignettes/frenchfish-vignette.Rmd
```

Disclaimer
----
Note that this is prerelease software. Please use accordingly.
