#this R script does all the automatic code generation (for package maintainers).
#You need to call this script if you change:
#1. Rcpp functions (i.e. add, remove, change a function tagged with \\ [[Rcpp::export]])
#2. Roxygen documentation (including the roxygen import statements)
#it should be located in the scripts subdirectory of the bamsignals package
#and called like this:
#R --slave < bamsignals/scripts/generateCode.R
#after that you can use R CMD build and install the tar.gz

load_or_get <- function(pkg){
	if (!require(pkg, character.only=T)) install.packages(pkg, dependencies=TRUE)
	library(pkg, character.only=T)
}

#wrap C++ functions into R functions
pckg <- "bamsignals"
load_or_get("Rcpp")
compileAttributes(pckg)

#create man pages and the NAMESPACE file
load_or_get("devtools")
devtools::install_deps(pckg, dependencies=TRUE)
load_or_get("roxygen2")
roxygenize(pckg)

# automatically register functions. Needs the Kmisc package
# from github, NOT from CRAN (it's behind of many commits).
# To install it, do:
# devtools::install_github("kevinushey/Kmisc")
# the last commit that was working was:
# x=b9a340ecef5bdd61fb09f462c179ccb86aff4e31, so in case
# the newest version doesn't work, a quick fix would be
# devtools::install_github("kevinushey/Kmisc", rev=x)
if (!require("Kmisc")){
	load_or_get("devtools")
	devtools::install_github("kevinushey/Kmisc")
}

setwd(pckg)
Kmisc:::registerFunctions("")
setwd("../")
