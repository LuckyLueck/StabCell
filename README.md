# StabCell

StabCell allows the generation of stable clusterings and sets of differentially
expressed genes. It applies concepts from stability selection to arbitrary
clustering and marker detection methods and offers an integrated approach for 
modelling.


## Installation

You can install `StabCell` directly from [GitHub](https://github.com/LuckyLueck/StabCell):

```r
if(!require("devtools")){ install.packages("devtools") }
devtools::install_github("LuckyLueck/StabCell")
```

or if you also want to use the package vignette (see below) use

```r
if(!require("devtools")){ install.packages("devtools") }
devtools::install_github("LuckyLueck/StabCell", build_vignettes=TRUE)
```

## Usage

For an example use case see

```r
vignette("StabCell_Example", package="StabCell")
```
