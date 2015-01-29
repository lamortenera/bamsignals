# testthat script for bamsignals::count
#
# 2015-01-28

library(bamsignals)
context("bamsignals::count()")


#
# Artificial data preparation
#
binwidth <- 500
regionlength <- 10*binwidth

#1. generate list of reads and sort by starting position
readpos <- matrix( rnbinom( 2 * regionlength, 50, .99), nrow=2 )
reads <- list()
for (i in 1:ncol(readpos)) {
	if (readpos[1,i]!=0) {
		for (j in 1:readpos[1,i]) {
			reads <- append(reads, list(c(i, min(i+50-1, regionlength), 0, floor(runif(1, 0, 50)))))
		}
	}
	if (readpos[2,i]!=0) {
		for (j in 1:readpos[2,i]) {
			reads <- append(reads, list(c(max(i-50+1,0), i, 16, floor(runif(1, 0, 50)))))
		}
	}
}
out <- do.call("rbind", reads)
out <- out[ order(out[,1]),]

#2. write to SAM file
dir.create("data", showWarnings=F)
sam.file.se <- "data/bamsignals_SE_counts.sam"
writeLines("@HD\tVN:1.0\tSO:unsorted\n@SQ\tSN:chr1\tLN:249250621", sam.file.se)
write.table(x=cbind("*", out[,3], "chr1", out[,1], out[,4], "50M", "*", "0", "0", "*", "*"),
			file=sam.file.se,
			sep="\t",
			quote=F,
			row.names=F,
			col.names=F,
			append=T
			)

#3. SAM -> BAM, BAM->BAMi
bamsignals:::writeSamAsBamAndIndex(sam.file.se, "data/bamsignals_SE_counts.bam")

#
# Test cases
#
test_that("Test that count works for Single End Data", {
		  #gr <- GRanges(seqnames="chr1", ranges=IRanges(start=1:300,end=bed[[3]]), strand=bed[[6]])
		  #TODO continue here with counting from bamfile and comparing to readpos infered counts

})

test_that("Test that count works for Paired End Data", {

})

#bed <- scan("data/regions.bed", what=list(character(), numeric(), numeric(), character(), character(), character()), skip="#")
#
## Single End data
#bampath <- "data/Toy_SE.bam"
#
## Count with bamsignals
#count.unstranded      <- count(gr, bampath, ss=F)
#count.stranded        <- count(gr, bampath, ss=T)
#count.shift           <- count(gr, bampath, shift=100, ss=F)
#count.qual            <- count(gr, bampath, mapqual=40, ss=F)
#
## Compare to truth
#bedtools.Toy_SE <- read.delim("data/Toy_SE.counts", header=T)
#
#all( count.unstranded == bedtools.Toy_SE$FivePrime_count )
#all( count.stranded[1,] == bedtools.Toy_SE$FivePrime_count_sense )
#all( count.stranded[2,] == bedtools.Toy_SE$FivePrime_count_antisense )
#all( count.shift == bedtools.Toy_SE$FivePrime_count_100_shift )
#all( count.qual == bedtools.Toy_SE$FivePrime_count_min40qual )
#
## Paired End data
#bampath <- "data/Toy_PE.bam"
#
## Count with bamsignals
#count.unstranded      <- count(gr, bampath, ss=F, paired.end=T, paired.end.midpoint=F)
#count.stranded        <- count(gr, bampath, ss=T, paired.end=T, paired.end.midpoint=F)
#count.shift           <- count(gr, bampath, ss=F, paired.end=T, shift=100, paired.end.midpoint=F)
#count.qual            <- count(gr, bampath, ss=F, paired.end=T, mapqual=40, paired.end.midpoint=F)
#count.midpoint        <- count(gr, bampath, ss=F, paired.end=T, paired.end.midpoint=T)
#count.midpoint.stranded <- count(gr, bampath, ss=T, paired.end=T, paired.end.midpoint=T)
#count.midpoint.shift  <- count(gr, bampath, ss=F, paired.end=T, paired.end.midpoint=T, shift=100)
#count.midpoint.qual   <- count(gr, bampath, ss=F, paired.end=T, paired.end.midpoint=T, mapqual=40)
#
## Compare to truth
#bedtools.Toy_PE <- read.delim("data/Toy_PE.counts", header=T)
#
#all( count.unstranded == bedtools.Toy_PE$FivePrime_count )
#all( count.stranded[1,] == bedtools.Toy_PE$FivePrime_count_sense )
#all( count.stranded[2,] == bedtools.Toy_PE$FivePrime_count_antisense )
#all( count.shift == bedtools.Toy_PE$FivePrime_count_100_shift )
#all( count.qual == bedtools.Toy_PE$FivePrime_count_min40qual )
#all( count.midpoint == bedtools.Toy_PE$MidPoint_count )
#all( count.midpoint.stranded[1,] == bedtools.Toy_PE$MidPoint_count_sense )
#all( count.midpoint.stranded[2,] == bedtools.Toy_PE$MidPoint_count_antisense )
#all( count.midpoint.shift == bedtools.Toy_PE$MidPoint_count_100_shift )
#all( count.midpoint.qual == bedtools.Toy_PE$MidPoint_count_min40qual )
