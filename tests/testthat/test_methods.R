context("bamsignals methods")

#add a strand specific shift to each range
ssshift <- function(gr, shift){
	shifts <- rep(shift, length(gr))
	shifts[as.logical(strand(gr)=="-")] <- -shift
	
	GenomicRanges::shift(gr, shifts)
}

#implement the count function in R
countR <- function(genes, reads, shift=0, ss=FALSE){
	#apply the shift to the reads
	reads <- ssshift(reads, shift)
	
	ov <- findOverlaps(genes, reads, select="all", type="any", ignore.strand=TRUE)
	
	ssCounts <- sapply(1:length(genes), function(g){
		gStart <- start(genes)[g]
		gEnd <- end(genes)[g]
		#get the reads overlapping with this gene
		ovReads <- reads[subjectHits(ov)[queryHits(ov)==g]]
		#sum up positive reads
		preads <- sum(start(ovReads)>=gStart & start(ovReads)<=gEnd & strand(ovReads)=="+")
		#sum up negative reads
		nreads <- sum(end(ovReads)>=gStart & end(ovReads)<=gEnd & strand(ovReads)=="-")
		c(preads, nreads)
	})
	
	#reverse columns of negative regions
	toRev <- as.logical(strand(genes)=="-")
	ssCounts[,toRev] <- ssCounts[c(2,1),toRev]
	rownames(ssCounts) <- c("sense", "antisense")
	if (!ss) return(colSums(ssCounts))
	
	ssCounts
}

#implement the pileup function in R
pileupR <- function(genes, reads, shift=0, ss=FALSE){
	#apply the shift to the reads
	reads <- ssshift(reads, shift)
	
	ov <- findOverlaps(genes, reads, select="all", type="any", ignore.strand=TRUE)
	
	isNegGene <- as.logical(strand(genes)=="-")
	lapply(seq_along(genes), function(g){
		gStart <- start(genes)[g]
		gEnd <- end(genes)[g]
		gLen <- gEnd-gStart+1
		#get the reads overlapping with this gene
		ovReads <- reads[subjectHits(ov)[queryHits(ov)==g]]
		isNeg <- as.logical(strand(ovReads)=="-")
		#sum up positive reads
		pStarts <- start(ovReads)[!isNeg]
		preads <- table(factor(pStarts-gStart+1, levels=1:gLen))
		#sum up negative reads
		nEnds <- end(ovReads)[isNeg]
		nreads <- table(factor(nEnds-gStart+1, levels=1:gLen))
		
		mat <- t(cbind(preads, nreads))
		if (isNegGene[g]) {
			#take into account the strand of the region
			mat <- rev(mat)
		}
		
		mat <- as.integer(mat)
		dim(mat) <- c(2, length(preads))
		dimnames(mat) <- NULL
		
		if (ss){#make matrix
			rownames(mat) <- c("sense", "antisense")
			colnames(mat) <- NULL
			return(mat)
		} else return(colSums(mat))
	})
}

depthR <- function(genes, reads){
	isNegGene <- as.logical(strand(genes)=="-")
	seqNames <- as.character(seqnames(genes))
	
	#use the already implemented coverage function
	cvrg <- coverage(reads)
	lapply(seq_along(genes), function(g){
		sig <- cvrg[[seqNames[g]]][start(genes)[g]:end(genes)[g]]
		#take into account strand of the region
		if (isNegGene[g]) sig <- rev(sig)
		as.integer(sig)
	})
}


## TOY DATA ##
library(GenomicRanges)
library(Rsamtools)
bampath <- system.file("extdata", "randomBam.bam", package="bamsignals")
genes <- get(load(system.file("extdata", "randomAnnot.Rdata", package="bamsignals")))

#load the reads with Rsamtools
reads <- scanBam(bampath)[[1]]
#make a GRanges object
reads <- GRanges(seqnames=reads$rname, strand=reads$strand, IRanges(start=reads$pos, width=reads$qwidth))

test_that("count function", {
	shifts <- c(0, 100)
	sss <- c(TRUE, FALSE)
	for (shift in shifts){
		for (ss in sss){
			expect_equal(countR(genes, reads, ss=ss, shift=shift), 
			count(genes, bampath, ss=ss, shift=shift, verbose=FALSE))
		}
	}
})

test_that("pileup function", {
	shifts <- c(0, 100)
	sss <- c(TRUE, FALSE)
	for (shift in shifts){
		for (ss in sss){
			expect_equal(pileupR(genes, reads, ss=ss, shift=shift), 
			as.list(bamsignals::pileup(genes, bampath, ss=ss, shift=shift, verbose=FALSE)))
		}
	}
})

test_that("depth function", {
	expect_equal(depthR(genes, reads), as.list(depth(genes, bampath, verbose=FALSE)))
})

