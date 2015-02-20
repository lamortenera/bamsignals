context("bamsignals methods")
library(GenomicRanges)
library(Rsamtools)
source("utils.R")


## GENERATE TOY DATA ##
##reads
nRef <- 3  #number of chromosomes
lastStart <- 1e3 #last start of a read pair (i.e. chromosome length more or less)
nPairs <- 1e5 #number of read pairs
avgRLen <- 30 #average read length
avgFLen <- 150 #average fragment length or template length
nRemove <- 1e3 #remove random reads from the pairs
bampath <- paste0(tempfile(), ".bam") #path where to save the bam file

##regions
nRegions <- 20		#number of regions
avgRegLen <- 200	#average region length


#generate reads
posReads <- data.frame(
	rname=sample(nRef, nPairs, replace=TRUE),			#chromosome name
	pos=sample(lastStart, nPairs, replace=TRUE),	#start of the read
	qwidth=1+rpois(nPairs, lambda=(avgRLen-1)),		#length of the read
	strand=rep("+", nPairs), 											#strand
	isize=1+rpois(nPairs, lambda=(avgFLen-1)),		#template (or fragment) length
	read1=(sample(2, nPairs, replace=TRUE)==1),		#is it the first read in the pair?
	mapq=sample(254, nPairs, replace=TRUE))				#mapping quality

negReads <- data.frame(
	rname=posReads$rname,
	qwidth=1+rpois(nPairs, lambda=(avgRLen-1)),
	strand=rep("-", nPairs),
	isize=-posReads$isize,
	read1=!posReads$read1,
	mapq=sample(254, nPairs, replace=TRUE))
negReads$pos <- posReads$pos + posReads$isize - negReads$qwidth

reads <- rbind(posReads, negReads)
#remove some of them
reads <- reads[-sample(nrow(reads), nRemove),]

#write to a bam file
readsToBam(reads, bampath)

#generate regions
regions <- GRanges(
	seqnames=sample(nRef, nRegions, replace=TRUE),
	strand=sample(c("+", "-"), nRegions, replace=TRUE),
	ranges=IRanges(
		start=sample(lastStart, nRegions, replace=TRUE), 
		width=1+rpois(nRegions, lambda=(avgRegLen-1))))



getPem <- function(pe){
	if (pe) return(c(TRUE, FALSE))
	FALSE
}

argsToStr <- function(args){
	vals <- sapply(args, get, envir=globalenv())
	paste(collapse=",", sep="=", args, vals)
}

test_that("bamCount function", {
	for (shift in c(0, 100)){
		for (mapq in c(0, 100)){
			for (ss in c(TRUE, FALSE)){
				for (pe in c(TRUE, FALSE)){
					for (pem in getPem(pe)){
						expect_equal(
							label=paste0("bamCount{", 
							paste("shift", shift, "mapq", mapq, "ss", ss, "pe", pe, "pem", pem, sep="="),
							"}"),
							countR(regions, reads, ss=ss, shift=shift, paired.end=pe, paired.end.midpoint=pem, mapqual=mapq), 
							bamCount(regions, bampath, ss=ss, shift=shift, paired.end=pe, paired.end.midpoint=pem, mapqual=mapq, verbose=FALSE))
	}	}	}	}	}
})

test_that("bamProfile function", {
	for (shift in c(0, 100)){
		for (mapq in c(0, 100)){
			for (ss in c(TRUE, FALSE)){
				for (pe in c(TRUE, FALSE)){
					for (pem in getPem(pe)){
						expect_equal(
							label=paste0("bamProfile{", 
							paste("shift", shift, "mapq", mapq, "ss", ss, "pe", pe, "pem", pem, sep="="),
							"}"),
							profileR(regions, reads, ss=ss, shift=shift, paired.end=pe, paired.end.midpoint=pem, mapqual=mapq), 
							as.list(bamProfile(regions, bampath, ss=ss, shift=shift, paired.end=pe, paired.end.midpoint=pem, mapqual=mapq, verbose=FALSE)))
	}	}	}	}	}
})


test_that("bamCoverage function", {
	for (mapq in c(0, 100)){
		for (pe in c(TRUE, FALSE)){
				expect_equal(
					label=paste0("bamCoverage{mapq=",mapq, ", pe=", pe, "}"),
					coverageR(regions, reads, paired.end=pe, mapqual=mapq), 
					as.list(bamCoverage(regions, bampath, paired.end=pe, mapqual=mapq, verbose=FALSE)))
	}	}
})

#remove the temporary files
unlink(paste0(bampath, c("", ".bai")))

