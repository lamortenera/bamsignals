#this R script does all the automatic code generation (for package maintainers).
#You need to call this script if you change:
#1. Rcpp functions (i.e. add, remove, change a function tagged with \\ [[Rcpp::export]])
#2. Roxygen documentation (including the roxygen import statements)
#it should be located in the scripts subdirectory of the bamsignals package
#and called like this:
#R --slave < bamsignals/scripts/generateCode.R
#after that you can use R CMD build and install the tar.gz

#wrap C++ functions into R functions
pckg <- "bamsignals"
library(Rcpp)
compileAttributes(pckg)

#create man pages and the NAMESPACE file
library(roxygen2)
roxygenize(pckg)

#register function (needs the Kmisc package from CRAN)
setwd(pckg)
Kmisc:::registerFunctions("")
setwd("../")
