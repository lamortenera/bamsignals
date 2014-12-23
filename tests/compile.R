library(roxygen2)
library(Rcpp)

# If interactive, you should uncomment these lines prior to compiling
#detach(package:bamsignals, unload=T)

compileAttributes("../../bamsignals")
roxygenize("../../bamsignals")
system("R CMD build ../../bamsignals")
all( install.packages("bamsignals_1.1.tar.gz") )
unlink( "bamsignals_1.1.tar.gz" )

# Creating bedtools reference counts and profiles
system("cd data && sh README.txt")
