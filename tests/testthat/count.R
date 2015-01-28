# A naive test script for bamsignals::count consistency
#
# Run this script with
# R --slave < count.R
#
# This test script requires the Roxygen2 R package. Roxygen2 is available through CRAN. You can install the package with intall.packages("roxygen2").
#
# 2014-04-15
require(bamsignals)
require(GenomicRanges)
bed <- scan("data/regions.bed", what=list(character(), numeric(), numeric(), character(), character(), character()), skip="#")
gr <- GRanges(seqnames=bed[[1]], ranges=IRanges(start=bed[[2]]+1,end=bed[[3]]), strand=bed[[6]])

# Single End data
bampath <- "data/Toy_SE.bam"

# Count with bamsignals
count.unstranded      <- count(gr, bampath, ss=F)
count.stranded        <- count(gr, bampath, ss=T)
count.shift           <- count(gr, bampath, shift=100, ss=F)
count.qual            <- count(gr, bampath, mapqual=40, ss=F)

# Compare to truth
bedtools.Toy_SE <- read.delim("data/Toy_SE.counts", header=T)

all( count.unstranded == bedtools.Toy_SE$FivePrime_count )
all( count.stranded[1,] == bedtools.Toy_SE$FivePrime_count_sense )
all( count.stranded[2,] == bedtools.Toy_SE$FivePrime_count_antisense )
all( count.shift == bedtools.Toy_SE$FivePrime_count_100_shift )
all( count.qual == bedtools.Toy_SE$FivePrime_count_min40qual )

# Paired End data
bampath <- "data/Toy_PE.bam"

# Count with bamsignals
count.unstranded      <- count(gr, bampath, ss=F, paired.end=T, paired.end.midpoint=F)
count.stranded        <- count(gr, bampath, ss=T, paired.end=T, paired.end.midpoint=F)
count.shift           <- count(gr, bampath, ss=F, paired.end=T, shift=100, paired.end.midpoint=F)
count.qual            <- count(gr, bampath, ss=F, paired.end=T, mapqual=40, paired.end.midpoint=F)
count.midpoint        <- count(gr, bampath, ss=F, paired.end=T, paired.end.midpoint=T)
count.midpoint.stranded <- count(gr, bampath, ss=T, paired.end=T, paired.end.midpoint=T)
count.midpoint.shift  <- count(gr, bampath, ss=F, paired.end=T, paired.end.midpoint=T, shift=100)
count.midpoint.qual   <- count(gr, bampath, ss=F, paired.end=T, paired.end.midpoint=T, mapqual=40)

# Compare to truth
bedtools.Toy_PE <- read.delim("data/Toy_PE.counts", header=T)

all( count.unstranded == bedtools.Toy_PE$FivePrime_count )
all( count.stranded[1,] == bedtools.Toy_PE$FivePrime_count_sense )
all( count.stranded[2,] == bedtools.Toy_PE$FivePrime_count_antisense )
all( count.shift == bedtools.Toy_PE$FivePrime_count_100_shift )
all( count.qual == bedtools.Toy_PE$FivePrime_count_min40qual )
all( count.midpoint == bedtools.Toy_PE$MidPoint_count )
all( count.midpoint.stranded[1,] == bedtools.Toy_PE$MidPoint_count_sense )
all( count.midpoint.stranded[2,] == bedtools.Toy_PE$MidPoint_count_antisense )
all( count.midpoint.shift == bedtools.Toy_PE$MidPoint_count_100_shift )
all( count.midpoint.qual == bedtools.Toy_PE$MidPoint_count_min40qual )
