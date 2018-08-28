bamsignals
==========
R package to quickly obtain count vectors from indexed bam files.

## Installation

Install the `devtools` package to be able to directly install R packages hosted on github :

```R
install.packages("devtools")
```

To install `bamsignals` type:

```R
devtools::install_github("lamortenera/bamsignals")
```

Alternatively, you can install it from Bioconductor (the two repos are synchronized):

```R
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("bamsignals")
```

