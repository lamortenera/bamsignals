#get a CountSignals object
library(GenomicRanges)
bampath <- 
system.file("extdata", "randomBam.bam", package="bamsignals")
genes <- 
get(load(system.file("extdata", "randomAnnot.Rdata", package="bamsignals")))
csig <- bamProfile(bampath, genes, ss=TRUE)

#show it
show(csig)

#number of contained signals
len <- length(csig)

#width of each signal
w <- width(csig)

#get one element as a vector (or matrix)
v <- csig[1]

#use as if it was a list
tot_per_sig <- sapply(csig, sum)

#convert to a list
siglist <- as.list(csig)

#get regions and signals of the same width
proms <- promoters(genes, upstream=150, downstream=150)
csig <- bamCoverage(bampath, proms)

#convert to matrix
mat <- alignSignals(csig)
