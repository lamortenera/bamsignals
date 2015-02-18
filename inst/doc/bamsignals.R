## ----style, echo = FALSE, results = 'asis'-------------------------------
BiocStyle::markdown()

## ----message = FALSE-----------------------------------------------------
library(GenomicRanges)
library(Rsamtools)
library(bamsignals)

## ------------------------------------------------------------------------
bampath <- system.file("extdata", "randomBam.bam", package="bamsignals")
genes <- get(load(system.file("extdata", "randomAnnot.Rdata", package="bamsignals")))
genes

## ------------------------------------------------------------------------
#sequence names of the GenomicRanges object
seqinfo(genes)
#sequence names in the bam file
bf <- Rsamtools::BamFile(bampath)
seqinfo(bf)
#checking if there is an index
file.exists(gsub(".bam$", ".bam.bai", bampath))


## ------------------------------------------------------------------------
proms <- GenomicRanges::promoters(genes, upstream=100, downstream=100)
counts <- count(proms, bampath, verbose=FALSE)
str(counts)

## ------------------------------------------------------------------------
counts <- count(proms, bampath, verbose=FALSE, shift=75)
str(counts)

## ------------------------------------------------------------------------
strand(proms)
counts <- count(proms, bampath, verbose=FALSE, ss=TRUE)
str(counts)

## ------------------------------------------------------------------------
sigs <- bamsignals::pileup(genes, bampath, verbose=FALSE)
sigs

## ------------------------------------------------------------------------
#CountSignals is conceptually like a list
lsigs <- as.list(sigs)
stopifnot(length(lsigs[[1]])==length(sigs[1]))
#sapply and lapply can be used as if we were using a list
stopifnot(all(sapply(sigs, sum) == sapply(lsigs, sum)))

## ------------------------------------------------------------------------
stopifnot(all(width(sigs)==width(genes)))

## ------------------------------------------------------------------------
sssigs <- bamsignals::pileup(genes, bampath, verbose=FALSE, ss=TRUE)
sssigs

## ------------------------------------------------------------------------
str(sssigs[1])
#summing up the counts from the two strands is the same as using ss=FALSE
stopifnot(colSums(sssigs[1])==sigs[1])
#the width function takes into account that now the signals are strand-specific
stopifnot(width(sssigs)==width(sigs))
#the length function does not, a strand-specific signal is twice as long as before
stopifnot(length(sssigs[1])==2*length(sigs[1]))

## ------------------------------------------------------------------------
xlab <- "offset from start of the region"
ylab <- "counts per base pair (negative means antisense)"
main <- paste0("read profile of the region ", seqnames(genes)[1], ":", start(genes)[1], "-", end(genes)[1])
plot(sigs[1], ylim=c(-max(sigs[1]), max(sigs[1])), ylab=ylab, xlab=xlab, main=main, type="l")
lines(sssigs[1]["sense",], col="blue")
lines(-sssigs[1]["antisense",], col="red")
legend("topright", c("sense", "antisense", "both"), col=c("blue", "red", "black"), lty=1)

## ------------------------------------------------------------------------
#The promoter regions have all the same width
sigs <- bamsignals::pileup(proms, bampath, ss=FALSE, verbose=FALSE)
sssigs <- bamsignals::pileup(proms, bampath, ss=TRUE, verbose=FALSE)

sigsMat <- alignSignals(sigs)
sigsArr <- alignSignals(sssigs)


## ------------------------------------------------------------------------
#the dimensions are [base pair, region]
str(sigsMat)
#the dimensions are [strand, base pair, region]
str(sigsArr)

stopifnot(all(sigsMat == sigsArr["sense",,] + sigsArr["antisense",,]))

## ------------------------------------------------------------------------
avgSig <- rowMeans(sigsMat)
avgSenseSig <- rowMeans(sigsArr["sense",,])
avgAntisenseSig <- rowMeans(sigsArr["antisense",,])
ylab <- "average counts per base pair"
xlab <- "distance from TSS"
main <- paste0("average profile of ", length(proms), " promoters")
xs <- -99:100
plot(xs, avgSig, ylim=c(0, max(avgSig)), xlab=xlab, ylab=ylab, main=main, type="l")
lines(xs, avgSenseSig, col="blue")
lines(xs, avgAntisenseSig, col="red")
legend("topright", c("sense", "antisense", "both"), col=c("blue", "red", "black"), lty=1)

## ------------------------------------------------------------------------

binsize <- 20
binnedSigs <- bamsignals::pileup(proms, bampath, binsize=binsize)
stopifnot(all(width(binnedSigs)==ceiling(width(sigs)/binsize)))
binnedSigs

## ------------------------------------------------------------------------
avgBinnedSig <- rowMeans(alignSignals(binnedSigs))
#the counts in the bin are the sum of the counts in each base pair
stopifnot(all.equal(colSums(matrix(avgSig, nrow=binsize)),avgBinnedSig))
#let's plot it
ylab <- "average counts per base pair"
plot(xs, avgSig, xlab=xlab, ylab=ylab, main=main, type="l")
lines(xs, rep(avgBinnedSig, each=binsize)/binsize, lty=2)
legend("topright", c("base pair count", "bin count"), lty=c(1, 2))


## ------------------------------------------------------------------------
covSigs <- bamsignals::depth(genes, bampath, verbose=FALSE)
puSigs <- bamsignals::pileup(genes, bampath, verbose=FALSE)
xlab <- "offset from start of the region"
ylab <- "reads per base pair"
main <- paste0("read coverage and profile of the region ", seqnames(genes)[1], ":", start(genes)[1], "-", end(genes)[1])
plot(covSigs[1], ylim=c(0, max(covSigs[1])), ylab=ylab, xlab=xlab, main=main, type="l")
lines(puSigs[1], lty=2)
legend("topright", c("covering the base pair", "5' end maps to the base pair"), lty=c(1,2))

