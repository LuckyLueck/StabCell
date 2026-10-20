# StabCell

StabCell allows the generation of stable clusterings and sets of differentially expressed genes. It applies concepts from stability selection to arbitrary clustering and marker detection methods and offers an integrated approach for modelling.

<img src="img/schematic.svg" width="600"/>

## Installation

You can install `StabCell` directly from [GitHub](https://github.com/LuckyLueck/StabCell):

``` r
if(!require("devtools")){ install.packages("devtools") }
devtools::install_github("LuckyLueck/StabCell")
```

## Usage

For an example use case see the [marker detection vignette](vignettes/StabCell_Marker_Detection_Example.Rmd).

To reproduce the results published in the article refer to the [StabCell_paper repository](https://github.com/LuckyLueck/StabCell_paper).

## References

[1] [Lück, Niklas, Rossi, Andrea, Staerk, Christian. "StabCell: Stability selection for clustering and marker detection in single-cell RNA sequencing" bioRxiv (2026)](https://doi.org/10.64898/2026.05.07.720061).

[2] [Meinshausen, Nicolai, Bühlmann, Peter. "Stability Selection" Journal of the Royal Statistical Society Series B: Statistical Methodology (2012)](https://doi.org/10.1111/j.1467-9868.2010.00740.x).

[3] [Shah, Rajen D., Samworth, Richard J. "Variable Selection with Error Control: Another Look at Stability Selection" Journal of the Royal Statistical Society Series B: Statistical Methodology](https://doi.org/10.1111/j.1467-9868.2011.01034.x)
