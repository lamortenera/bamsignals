# A naive test script for bamsignals::pileup consistency
#
# Run this script with
# R --slave < pileup.R
#
# This test script requires the Roxygen2 R package. Roxygen2 is available through CRAN. You can install the package with intall.packages("roxygen2").
#
# 2014-04-22
require(bamsignals)
require(GenomicRanges)
bed <- scan("data/regions.bed", what=list(character(), numeric(), numeric(), character(), character(), character()), skip="#")
gr <- GRanges(seqnames=bed[[1]], ranges=IRanges(start=bed[[2]]+1,end=bed[[3]]), strand=bed[[6]])

# Single End data
bampath <- "data/Toy_SE.bam"

# Count with bamsignals
pu.unstranded      <- pileup(gr, bampath, ss=F)
pu.stranded        <- pileup(gr, bampath, ss=T)
pu.stranded.binned <- pileup(gr, bampath, ss=T, binsize=5)
pu.shift           <- pileup(gr, bampath, ss=F, shift=100)
pu.qual            <- pileup(gr, bampath, ss=F, mapqual=40)

#self consistency
idxs <- bamsignals:::getIdxRanges(pu.unstranded, gr, gr[1:2], strand=T)
l <- length(pu.stranded$counts)
r1 <- pu.stranded.binned$counts[pu.stranded.binned$starts[1]:pu.stranded.binned$ends[1]]
a1 <- matrix(r1, ncol=2, byrow=T)
r2 <- pu.stranded$counts[pu.stranded$starts[1]:pu.stranded$ends[1]]
a2temp <- array(r2, dim=c(2,5,length(r1)/2))
a2 <- apply(a2temp, 1, colSums)
all(pu.unstranded$counts - pu.stranded$counts[2*(1:l)] - pu.stranded$counts[2*(1:l)-1]==0)
all(a2-a1==0)
all(!is.na(pu.stranded.binned$counts[pu.stranded.binned$starts[3]:pu.stranded.binned$ends[3]]))
all(as.numeric(getSignal(pu.stranded, 2)) == as.numeric(pu.stranded$counts)[pu.stranded$starts[2]:pu.stranded$ends[2]])

# Compare to correct counts
bedtools.Toy_SE <- read.delim("data/Toy_SE.profile", header=T)

all( pu.unstranded$counts == bedtools.Toy_SE$FivePrime_profile )
all(sapply(1:length(gr), function(i) { all( getSignal(pu.unstranded, i) == bedtools.Toy_SE$FivePrime_profile[ (pu.unstranded$starts[i]):(pu.unstranded$ends[i]) ]) }) )
all( pu.stranded$counts[1,] == bedtools.Toy_SE$FivePrime_profile_sense )
all(sapply(1:length(gr), function(i) { all( getSignal(pu.stranded, i)[1,] == bedtools.Toy_SE$FivePrime_profile_sense[ ((pu.stranded$starts[i]-1)/2 +1):(pu.stranded$ends[i]/2) ]) }) )
all( pu.stranded$counts[2,] == bedtools.Toy_SE$FivePrime_profile_antisense )
all(sapply(1:length(gr), function(i) { all( getSignal(pu.stranded, i)[2,] == bedtools.Toy_SE$FivePrime_profile_antisense[ ((pu.stranded$starts[i]-1)/2 +1):(pu.stranded$ends[i]/2) ]) }) )
all( pu.shift$counts == bedtools.Toy_SE$FivePrime_profile_100_shift )
all(sapply(1:length(gr), function(i) { all( getSignal(pu.shift, i) == bedtools.Toy_SE$FivePrime_profile_100_shift[ (pu.shift$starts[i]):(pu.shift$ends[i]) ]) }) )
all( pu.qual$counts == bedtools.Toy_SE$FivePrime_profile_min40qual )
all(sapply(1:length(gr), function(i) { all( getSignal(pu.qual, i) == bedtools.Toy_SE$FivePrime_profile_min40qual[ (pu.qual$starts[i]):(pu.qual$ends[i]) ]) }) )

# Paired End data
bampath <- "data/Toy_PE.bam"

# Count with bamsignals
pu.stranded        <- pileup(gr, bampath, ss=T, paired.end=T, paired.end.midpoint=F)
pu.stranded.binned <- pileup(gr, bampath, ss=T, binsize=5, paired.end=T, paired.end.midpoint=F)
pu.unstranded      <- pileup(gr, bampath, ss=F, paired.end=T, paired.end.midpoint=F)
pu.shift           <- pileup(gr, bampath, ss=F, shift=100, paired.end=T, paired.end.midpoint=F)
pu.qual            <- pileup(gr, bampath, ss=F, mapqual=40, paired.end=T, paired.end.midpoint=F)
pu.midpoint.stranded        <- pileup(gr, bampath, ss=T, paired.end=T, paired.end.midpoint=T)
pu.midpoint.stranded.binned <- pileup(gr, bampath, ss=T, binsize=5, paired.end=T, paired.end.midpoint=T)
pu.midpoint.unstranded      <- pileup(gr, bampath, ss=F, paired.end=T, paired.end.midpoint=T)
pu.midpoint.shift           <- pileup(gr, bampath, ss=F, shift=100, paired.end=T, paired.end.midpoint=T)
pu.midpoint.qual            <- pileup(gr, bampath, ss=F, mapqual=40, paired.end=T, paired.end.midpoint=T)

#self consistency
idxs <- bamsignals:::getIdxRanges(pu.unstranded, gr, gr[1:2], strand=T)
l <- length(pu.stranded$counts)
r1 <- pu.stranded.binned$counts[pu.stranded.binned$starts[1]:pu.stranded.binned$ends[1]]
a1 <- matrix(r1, ncol=2, byrow=T)
r2 <- pu.stranded$counts[pu.stranded$starts[1]:pu.stranded$ends[1]]
a2temp <- array(r2, dim=c(2,5,length(r1)/2))
a2 <- apply(a2temp, 1, colSums)
all(pu.unstranded$counts - pu.stranded$counts[2*(1:l)] - pu.stranded$counts[2*(1:l)-1]==0)
all(a2-a1==0)
all(!is.na(pu.stranded.binned$counts[pu.stranded.binned$starts[3]:pu.stranded.binned$ends[3]]))
all(as.numeric(getSignal(pu.stranded, 2)) == as.numeric(pu.stranded$counts)[pu.stranded$starts[2]:pu.stranded$ends[2]])

idxs <- bamsignals:::getIdxRanges(pu.midpoint.unstranded, gr, gr[1:2], strand=T)
l <- length(pu.midpoint.stranded$counts)
r1 <- pu.midpoint.stranded.binned$counts[pu.midpoint.stranded.binned$starts[1]:pu.midpoint.stranded.binned$ends[1]]
a1 <- matrix(r1, ncol=2, byrow=T)
r2 <- pu.midpoint.stranded$counts[pu.midpoint.stranded$starts[1]:pu.midpoint.stranded$ends[1]]
a2temp <- array(r2, dim=c(2,5,length(r1)/2))
a2 <- apply(a2temp, 1, colSums)
all(pu.midpoint.unstranded$counts - pu.midpoint.stranded$counts[2*(1:l)] - pu.midpoint.stranded$counts[2*(1:l)-1]==0)
all(a2-a1==0)
all(!is.na(pu.midpoint.stranded.binned$counts[pu.midpoint.stranded.binned$starts[3]:pu.midpoint.stranded.binned$ends[3]]))
all(as.numeric(getSignal(pu.midpoint.stranded, 2)) == as.numeric(pu.midpoint.stranded$counts)[pu.midpoint.stranded$starts[2]:pu.midpoint.stranded$ends[2]])

# Compare to correct counts
bedtools.Toy_PE <- read.delim("data/Toy_PE.profile", header=T)

all( pu.unstranded$counts == bedtools.Toy_PE$FivePrime_profile )
all(sapply(1:length(gr), function(i) { all( getSignal(pu.unstranded, i) == bedtools.Toy_PE$FivePrime_profile[ (pu.unstranded$starts[i]):(pu.unstranded$ends[i]) ]) }) )
all( pu.stranded$counts[1,] == bedtools.Toy_PE$FivePrime_profile_sense )
all(sapply(1:length(gr), function(i) { all( getSignal(pu.stranded, i)[1,] == bedtools.Toy_PE$FivePrime_profile_sense[ ((pu.stranded$starts[i]-1)/2 +1):(pu.stranded$ends[i]/2) ]) }) )
all( pu.stranded$counts[2,] == bedtools.Toy_PE$FivePrime_profile_antisense )
all(sapply(1:length(gr), function(i) { all( getSignal(pu.stranded, i)[2,] == bedtools.Toy_PE$FivePrime_profile_antisense[ ((pu.stranded$starts[i]-1)/2 +1):(pu.stranded$ends[i]/2) ]) }) )
all( pu.shift$counts == bedtools.Toy_PE$FivePrime_profile_100_shift )
all(sapply(1:length(gr), function(i) { all( getSignal(pu.shift, i) == bedtools.Toy_PE$FivePrime_profile_100_shift[ (pu.shift$starts[i]):(pu.shift$ends[i]) ]) }) )
all( pu.qual$counts == bedtools.Toy_PE$FivePrime_profile_min40qual )
all(sapply(1:length(gr), function(i) { all( getSignal(pu.qual, i) == bedtools.Toy_PE$FivePrime_profile_min40qual[ (pu.qual$starts[i]):(pu.qual$ends[i]) ]) }) )
all( pu.midpoint.unstranded$counts == bedtools.Toy_PE$MidPoint_profile )
all(sapply(1:length(gr), function(i) { all( getSignal(pu.midpoint.unstranded, i) == bedtools.Toy_PE$MidPoint_profile[ (pu.midpoint.unstranded$starts[i]):(pu.midpoint.unstranded$ends[i]) ]) }) )
all( pu.midpoint.stranded$counts[1,] == bedtools.Toy_PE$MidPoint_profile_sense )
all(sapply(1:length(gr), function(i) { all( getSignal(pu.midpoint.stranded, i)[1,] == bedtools.Toy_PE$MidPoint_profile_sense[ ((pu.midpoint.stranded$starts[i]-1)/2 +1):(pu.midpoint.stranded$ends[i]/2) ]) }) )
all( pu.midpoint.stranded$counts[2,] == bedtools.Toy_PE$MidPoint_profile_antisense )
all(sapply(1:length(gr), function(i) { all( getSignal(pu.midpoint.stranded, i)[2,] == bedtools.Toy_PE$MidPoint_profile_antisense[ ((pu.midpoint.stranded$starts[i]-1)/2 +1):(pu.midpoint.stranded$ends[i]/2) ]) }) )
all( pu.midpoint.shift$counts == bedtools.Toy_PE$MidPoint_profile_100_shift )
all(sapply(1:length(gr), function(i) { all( getSignal(pu.midpoint.shift, i) == bedtools.Toy_PE$MidPoint_profile_100_shift[ (pu.midpoint.shift$starts[i]):(pu.midpoint.shift$ends[i]) ]) }) )
all( pu.midpoint.qual$counts == bedtools.Toy_PE$MidPoint_profile_min40qual )
all(sapply(1:length(gr), function(i) { all( getSignal(pu.midpoint.qual, i) == bedtools.Toy_PE$MidPoint_profile_min40qual[ (pu.midpoint.qual$starts[i]):(pu.midpoint.qual$ends[i]) ]) }) )
