## TOY DATA ##
library(GenomicRanges)
bampath <- 
system.file("extdata", "randomBam.bam", package="bamsignals")
genes <- 
get(load(system.file("extdata", "randomAnnot.Rdata", package="bamsignals")))


## THE FUNCTION 'count' ##
#count how many reads map in each region (according to 5' end)
v <- bamCount(bampath, genes)
#plot it
labs <- paste0(seqnames(genes), ":", start(genes), "-", end(genes))
par(mar=c(5, 6, 4, 2))
barplot(v, names.arg=labs, main="read counts per region", las=2, 
	horiz=TRUE, cex.names=.6)

#distinguish between strands
v2 <- bamCount(bampath, genes, ss=TRUE)
#plot it
par(mar=c(5, 6, 4, 2))
barplot(v2, names.arg=labs, main="read counts per region", las=2, 
	horiz=TRUE, cex.names=.6, col=c("blue", "red"), legend=TRUE)



## THE FUNCTIONS 'bamProfile' and 'bamCoverage' ##
#count how many reads map to each base pair (according to 5' end)
pu <- bamProfile(bampath, genes)
#count how many reads cover each base pair
du <- bamCoverage(bampath, genes)
#plot it
xlab <- "offset from start of the region"
ylab <- "reads per base pair"
main <- paste0("read coverage and profile of the region ", labs[1])
plot(du[1], ylim=c(0, max(du[1])), ylab=ylab, xlab=xlab, main=main, type="l")
lines(pu[1], lty=2)
llab <- c("covering the base pair", "5' end maps to the base pair")
legend("topright", llab, lty=c(1,2), bg="white")



## REGIONS OF THE SAME SIZE AND OPTIONS FOR 'bamProfile' ##
proms <- promoters(genes, upstream=150, downstream=150)
#pileup according to strand
pu_ss <- bamProfile(bampath, proms, ss=TRUE)
#compute average over regions
avg_ss <- apply(alignSignals(pu_ss), 2, rowMeans)

#profile using a strand-specific shift
pu_shift <- bamProfile(bampath, proms, shift=75)
#compute average over regions
avg_shift <- rowMeans(alignSignals(pu_shift))

#profile using a strand-specific shift and a binning scheme
binsize <- 20
pu_shift_bin <- bamProfile(bampath, proms, shift=75, binsize=binsize)
#compute average over regions
avg_shift_bin <- rowMeans(alignSignals(pu_shift_bin))

#plot it
xs <- -149:150
main <- paste0("average read profile over ", length(proms), " promoters")
xlab <- "distance from TSS"
ylab <- "average reads per base pair"
plot(xs, avg_shift, xlab=xlab, ylab=ylab, main=main, type="l", 
	ylim=c(0, max(avg_shift)))
lines(xs, avg_ss["sense",], col="blue")
lines(xs, avg_ss["antisense",], col="red")
lines(xs, rep(avg_shift_bin/binsize, each=binsize), lty=2)
llabs <- 
c("sense reads", "antisense reads", "with shift", "binned and with shift")
legend("topright", llabs, col=c("blue", "red", "black", "black"),
	lty=c(1,1,1,2), bg="white")

