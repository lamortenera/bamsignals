#' Efficient counting of reads in a bam files
#'
#' Functions for extracting signals from a bam file. Differently than the other packages,
#' this package cannot be used to import reads in R.
#' All the read-processing is done in C/C++ and the only output are read counts. 
#'
#' @name bamsignals
#' @import Rcpp
#' @import IRanges
#' @import GenomicRanges
#' @docType package
#' @author Alessandro Mammana \email{mammana@@molgen.mpg.de}
#' @author Johannes Helmuth \email{helmuth@@molgen.mpg.de}
#' @useDynLib bamsignals
NULL
#' Pileup reads from a bam file.
#'
#' Compute read density in the regions specified by a GenomicRanges object.
#' A read position is always specified by its 5' end, so a read mapping to the reference strand
#' is positioned at its leftmost coordinate, a read mapping to the alternative strand
#' is positioned at its rightmost coordinate. To change that use the \code{shift} parameter
#' or the \code{coverage} function.
#' @param gr GenomicRanges object used to specify the regions
#' @param bampath path to the bam file storing the read. The file must be indexed. 
#' If a range is on the negative strand the profile will be reverse-complemented.
#' @param binsize If the value is set to 1, the method will return basepair-resolution read densities,
#' for bigger values the density profiles will be binned (and the memory requirements
#' will scale accordingly). 
#' @param mapqual discard reads with mapping quality strictly lower than this parameter.
#' The value 0 ensures that no read will be discarded, the value 254 that only reads
#' with the highest possible mapping quality will be considered.
#' @param shift shift the read position by a user-defined number of basepairs. This can
#' be handy in the analysis of chip-seq data.
#' @param ss produce a strand-specific profile or ignore the strand of the read. This option
#' does not have any effect on the strand of the region. Strand-specific profiles are
#' twice as long then strand-independent profiles.
#' @param format attempts to find a suitable matrix/array format for the count vector. 
#' if the profile is strand-specific one dimension will correspond to sense
#' and anti-sense strand, if the ranges have all the same width one dimension
#' will correspond to the range number.
#' @param paired.end a logical value indicating whether the bampath contains paired-end 
#' sequencing output. In this case, only first reads in proper mapped pairs are considered 
#' (FLAG 66).
#' @param paired.end.midpoint a logical value indicating whether the fragment middle 
#' points of each fragment should be counted. Therefore it uses the tlen information from
#' the given bam file (MidPoint = fragment_start + int( abs(tlen) / 2) )). For even tlen, 
#' the given midpoint is the round half down real midpoint.
#' @param paired.end.max.frag.length an integer indicating which fragments should be 
#' considered in paired-end sequencing data. Default value of 1,000 bases is generally
#' a good pick.
#' @param verbose a logical value indicating whether verbose output is desired
#' @return a list with the following arguments:
#' 	\item{counts}{the vector containing the read counts. This will be formatted
#' 	into a matrix or an array depending on whether the profile is strand-specific
#' 	and whether the ranges have all the same length.}
#' 	\item{starts, ends}{Vectors defining the boundaries of the count vector. 
#' 	To extract counts relative to the i-th range, use 
#'		\code{as.numeric(counts)[starts[i]:ends[i]]}, 
#' 	or the \code{getSignal} function to preserve the formatting.}
#'		\item{format}{This element is present if pu$counts is formatted
#' 	differently than a simple vector and it describes the formatting.}
#' @export
pileup <- function(gr, bampath, binsize=1, mapqual=0, shift=0, ss=F, format=T, paired.end=F, paired.end.midpoint=F, paired.end.max.frag.length=1000, verbose=T){
	if (verbose) {
		cat( "Processing ", bampath, " and ") 
		printStupidSentence()
	}
	if (binsize < 1){
		stop("provide a binsize greater or equal to 1")
	} else if (binsize > 1 && any((width(gr) %% binsize) != 0)){
		warning(paste("some ranges' widths are not a multiple of the selected binsize,",
		"some bins will correspond to less than binsize basepairs"))
	}

	pu <- pileup_core(gr, bampath, mapqual, binsize, shift, ss, paired.end, paired.end.midpoint, paired.end.max.frag.length)
	new("ConcatCounts", counts=pu$counts, starts=pu$starts, ends=pu$ends, ss=ss)
}

#' Compute read depth (or read coverage) from a bam file.
#'
#' Compute read coverage in the regions specified by a GenomicRanges object.
#' @param gr GenomicRanges object used to specify the regions. If a range
#' is on the negative strand its coverage will be reverse-complemented.
#' @param bampath path to the bam file storing the read. The file must be indexed.
#' @param mapqual discard reads with mapping quality strictly lower than this parameter.
#' The value 0 ensures that no read will be discarded, the value 254 that only reads
#' with the highest possible mapping quality will be considered.
#' @param format if the ranges have all the same width this method
#' will return a matrix.
#' @param paired.end a logical value indicating whether the bampath contains paired-end 
#' sequencing output. In this case, only first reads in proper mapped pairs are considered 
#' (FLAG 66).
#' @param paired.end.max.frag.length an integer indicating which fragments should be 
#' considered in paired-end sequencing data. Default value of 1,000 bases is generally
#' a good pick.
#' @param verbose a logical value indicating whether verbose output is desired
#' @return a list with the following arguments:
#' 	\item{counts}{the vector containing the read counts. This will be formatted
#' 	into a matrix depending on whether the ranges have all the same length.}
#' 	\item{starts, ends}{Vectors defining the boundaries of the count vector. 
#' 	To extract counts relative to the i-th range, use 
#'		\code{as.numeric(counts)[starts[i]:ends[i]]}, 
#' 	or the \code{getSignal} function to preserve the formatting.}
#'		\item{format} {This element is present if counts is formatted
#' 	differently than a simple vector and it describes the formatting.}
#' @export
depth <- function(gr, bampath, mapqual=0, format=T, paired.end=F, paired.end.max.frag.length=1000, verbose=T){
	if (verbose) {
		cat( "Processing ", bampath, " and ") 
		printStupidSentence()
	}
	pu <- coverage_core(gr, bampath, mapqual, paired.end, paired.end.max.frag.length)
	new("ConcatCounts", counts=pu$counts, starts=pu$starts, ends=pu$ends, ss=F)
}

#' Count reads from a bam file.
#'
#' Count reads in the bins specified by a GenomicRanges object.
#' A read position is always specified by its 5' end, so a read mapping to the reference strand
#' is positioned at its leftmost coordinate, a read mapping to the alternative strand
#' is positioned at its rightmost coordinate. To change that use the \code{shift} parameter.
#' @param gr GenomicRanges object used to specify the regions
#' @param bampath path to the bam file storing the read. The file must be indexed.
#' @param mapqual discard reads with mapping quality strictly lower than this parameter.
#' The value 0 ensures that no read will be discarded, the value 254 that only reads
#' with the highest possible mapping quality will be considered.
#' @param shift shift the read position by a user-defined number of basepairs. This can
#' be handy in the analysis of chip-seq data.
#' @param ss produce a strand-specific count or ignore the strand of the read. Strand-specific counts
#' will be returned in a 2*length(gr) matrix.
#' @param paired.end a logical value indicating whether the bampath contains paired-end 
#' sequencing output. In this case, only first reads in proper mapped pairs are considered 
#' (FLAG 66).
#' @param paired.end.midpoint a logical value indicating whether the fragment middle 
#' points of each fragment should be counted. Therefore it uses the tlen information from
#' the given bam file (MidPoint = fragment_start + int( abs(tlen) / 2) )). For even tlen, 
#' the given midpoint is the round half down real midpoint.
#' @param paired.end.max.frag.length an integer indicating which fragments should be 
#' considered in paired-end sequencing data. Default value of 1,000 bases is generally
#' a good pick.
#' @param verbose a logical value indicating whether verbose output is desired
#' @return a vector or a matrix with the counts
#' @export
count <- function(gr, bampath, mapqual=0, shift=0, ss=F, paired.end=F, paired.end.midpoint=F, paired.end.max.frag.length=1000, verbose=T){
	if (verbose) {
		cat( "Processing ", bampath, " and ") 
		printStupidSentence()
	}
	pu <- pileup_core(gr, bampath, mapqual, max(width(gr)), shift, ss, paired.end, paired.end.midpoint, paired.end.max.frag.length)
	if (ss) {
		dim(pu$counts) <- c(2, length(gr))
		rownames(pu$counts) <- c("sense", "antisense")
	}
	return(pu$counts)
}



printStupidSentence <- function(){
	sentences <- c(
	"tHaT'S tHa fAStEsT pIlE-uP bAm iN tHe SoUth!!!\n",
	"yOu cAn'T pIlE-Up FaStEr!!!\n",
	"I'M gOnNa cHaSe'em and PiLe'em aLl up!!!\n",
	"fOr brOoMmHiLdA!!!\n",
	"tHe lEgEnD said, hE cOuLd PiLe uP fAsTeR thAn LiGht\n",
	"I gEt gOoSeBuMPs wHen I seE yOu pilEuPpiNg...\n")
	cat(sample(sentences, 1))
}


strandToBool <- function(s){
	l <- levels(s)
	h <- logical(length(l))
	h[match(c("+","*","-"), l)] <- c(TRUE,TRUE,-FALSE)
	h[as.vector(s, mode="integer")]
}

xnor <- function(x, y) {(x & y) | !(x | y)}
 
getIdxRanges <- function(pu, pugr, subsetgr, strand=F){
	countsPerBase = (pu$end[1] - pu$start[1] + 1)/width(pugr[1])
	if (	length(pu$end) != length(pugr) || 
			any(pu$end-pu$start+1 != width(pugr)*countsPerBase)) stop("pu and pugr object do not match")
	
	ov <- findOverlaps(subsetgr, pugr, type="within", select="arbitrary")
	
	if (any(is.na(ov))) stop("the subset is not entirely contained in the given genomic range")
	if (any((start(subsetgr)-start(pugr)[ov]) %% (1/countsPerBase) != 0)) stop("the subset of the grange is not a subset of the bins") 
	cutstarts <- pu$start[ov] + (start(subsetgr) - start(pugr)[ov])*countsPerBase
	newwidths <- width(subsetgr)*countsPerBase
	
	if (strand) return( list(start=cutstarts, width=newwidths, strand=xnor(strandToBool(strand(subsetgr)) , strandToBool(strand(pugr)[ov]))))
	else return( list(start=cutstarts, width=newwidths) )
}


#these two functions are not exported and not documented.
#they basically query the pileup object in a certain index range.
subsetPileup <- function(pu, idx, start, width, strand=NULL){
	if (length(width)==1) width <- rep(width, length(idx))
	if (is.null(strand)) strand <- rep(TRUE, length(idx))
	subsetCounts(pu$counts, start + pu$starts[idx], width, strand)
}

countInPileup <- function(pu, idx, start, width){
	if (length(width)==1) width <- rep(width, length(idx))
	countInSubset(pu$counts, start + pu$starts[idx], width)
}

