# A naive test script for bamsignals::depth consistency
#
# Run this script with
# R --slave < depth.R
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
depth      <- depth(gr, bampath, paired.end=F)
depth.qual <- depth(gr, bampath, mapqual=40, paired.end=F)

#self consistency of the depth function
d1 <- depth(gr[1], bampath)
grtmp <- gr[1]
strand(grtmp[1]) = "-"
d2 <- depth(grtmp, bampath)

#self consistency
all(d1$counts == d2$counts[length(d2$counts):1])

# Compare to correct counts
bedtools.Toy_SE <- read.delim("data/Toy_SE.depth", header=T)

all( depth$counts == bedtools.Toy_SE$Depth )
all(sapply(1:length(gr), function(i) { all( getSignal(depth, i) == bedtools.Toy_SE$Depth[ (depth$starts[i]):(depth$ends[i]) ]) }) )
all( depth$counts == bedtools.Toy_SE$Depth_min40qual )
all(sapply(1:length(gr), function(i) { all( getSignal(depth, i) == bedtools.Toy_SE$Depth_min40qual[ (depth$starts[i]):(depth$ends[i]) ]) }) )

# Paired End data
bampath <- "data/Toy_PE.bam"

# Count with bamsignals
depth      <- depth(gr, bampath, paired.end=T)
depth.qual <- depth(gr, bampath, mapqual=40, paired.end=T)

# Compare to correct counts
bedtools.Toy_PE <- read.delim("data/Toy_PE.depth", header=T)

all( depth$counts == bedtools.Toy_PE$Depth )
all(sapply(1:length(gr), function(i) { all( getSignal(depth, i) == bedtools.Toy_PE$Depth[ (depth$starts[i]):(depth$ends[i]) ]) }) )
cat("Ranges not same are:\n")
which( !sapply(1:length(gr), function(i) { all( getSignal(depth, i) == bedtools.Toy_PE$Depth[ (depth$starts[i]):(depth$ends[i]) ]) }) )
all( depth.qual$counts == bedtools.Toy_PE$Depth_min40qual )
all(sapply(1:length(gr), function(i) { all( getSignal(depth.qual, i) == bedtools.Toy_PE$Depth_min40qual[ (depth.qual$starts[i]):(depth.qual$ends[i]) ]) }) )

#i=10; which( getSignal(depth, i) != bedtools.Toy_PE$Depth[ (depth$starts[i]):(depth$ends[i]) ])
#i=10; cbind( (start(gr)[i]-1):end(gr)[i], getSignal(depth, i), bedtools.Toy_PE$Depth[ (depth$starts[i]):(depth$ends[i]) ])

