#' Efficient counting of reads in a bam files
#'
#' Functions and classes for extracting signals from a bam file. Differently 
#' than the other packages, this package cannot be used to import reads in R.
#' All the read-processing is done in C/C++ and the only output are read counts. 
#'
#' @name bamsignals
#' @import Rcpp
#' @import IRanges
#' @import GenomicRanges
#' @import methods
#' @import BiocGenerics
#' @docType package
#' @author Alessandro Mammana \email{mammana@@molgen.mpg.de}
#' @author Johannes Helmuth \email{helmuth@@molgen.mpg.de}
#' @seealso \code{\link{bamsignals-methods}} for the functions and 
#' \code{\link{CountSignals}} for the class.
#' @useDynLib bamsignals, .registration=TRUE
NULL

#' Functions for extracting count signals from a bam file.
#'
#' @param gr GenomicRanges object used to specify the regions
#' @param bampath path to the bam file storing the read. The file must be
#' indexed. If a range is on the negative strand the profile will be 
#' reverse-complemented.
#' @param binsize If the value is set to 1, the method will return 
#' basepair-resolution read densities, for bigger values the density profiles
#' will be binned (and the memory requirements will scale accordingly). 
#' @param mapqual discard reads with mapping quality strictly lower than this 
#' parameter. The value 0 ensures that no read will be discarded, the value 254 
#' that only reads with the highest possible mapping quality will be considered.
#' @param shift shift the read position by a user-defined number of basepairs. 
#' This can be handy in the analysis of chip-seq data.
#' @param ss produce a strand-specific profile or ignore the strand of the read.
#' This option does not have any effect on the strand of the region. 
#' Strand-specific profiles are twice as long then strand-independent profiles.
#' @param paired.end a logical value indicating whether the bampath contains 
#' paired-end sequencing output. In this case, only first reads in proper mapped
#' pairs are considered. 
#' @param paired.end.midpoint a logical value indicating whether the fragment 
#' middle points of each fragment should be counted. Therefore it uses the tlen 
#' information from the given bam file 
#' (MidPoint = fragment_start + int( abs(tlen) / 2) )). For even tlen, 
#' the given midpoint will be moved of 0.5 basepairs in the 3' direction.
#' @param paired.end.max.frag.length an integer indicating which fragments 
#' should be considered in paired-end sequencing data. Default value of 1,000 
#' bases is generally a good pick.
#' @param verbose a logical value indicating whether verbose output is desired
#' @return \itemize{
#'   \item \code{bamProfile} and \code{bamCoverage}: a CountSignals object with a signal 
#'     per region
#'   \item \code{bamCount}: a vector with one element per region or, 
#'     if \code{ss==TRUE}, a matrix with one column per region and two rows 
#'     (sense and antisense).
#' }
#' @details A read position is always specified by its 5' end, so a read mapping
#' to the reference strand is positioned at its leftmost coordinate, a 
#' read mapping to the alternative strand is positioned at its rightmost 
#' coordinate. This can be changed using the \code{shift} parameter.
#' @seealso \code{\link{CountSignals}} for handling the return value of 
#' \code{bamProfile} and \code{bamCoverage}
#' @name bamsignals-methods
#' @example inst/examples/methods_example.R
NULL

#' @export
setGeneric("bamCount", function(gr, bampath, ...) standardGeneric("bamCount"))
#' \code{bamCount}: for each range, count the reads whose 5' end map in it.
#' @aliases bamCount
#' @rdname bamsignals-methods
#' @export
setMethod("bamCount", c("GenomicRanges", "character"), 
	function(gr, bampath, mapqual=0, shift=0, ss=FALSE, 
	paired.end=FALSE, paired.end.midpoint=FALSE, 
	paired.end.max.frag.length=1000, verbose=TRUE){
		if (verbose) printStupidSentence(bampath)
		
		pu <- pileup_core(gr, bampath, mapqual, max(width(gr)), shift, ss, 
		paired.end, paired.end.midpoint, paired.end.max.frag.length)
		if (ss) {
			dim(pu$counts) <- c(2, length(gr))
			rownames(pu$counts) <- c("sense", "antisense")
		}
		return(pu$counts)
	}
)


#' @export
setGeneric("bamProfile", function(gr, bampath, ...) standardGeneric("bamProfile"))
#' \code{bamProfile}: for each base pair in the ranges, compute the number of reads
#' whose 5' end maps there.
#' @aliases bamProfile
#' @rdname bamsignals-methods
#' @export
setMethod("bamProfile", c("GenomicRanges", "character"), 
	function(gr, bampath, binsize=1, mapqual=0, shift=0, ss=FALSE, 
	paired.end=FALSE, paired.end.midpoint=FALSE,
	 paired.end.max.frag.length=1000, verbose=TRUE){
		if (verbose) printStupidSentence(bampath)
		
		if (binsize < 1){
			stop("provide a binsize greater or equal to 1")
		} else if (binsize > 1 && any((width(gr) %% binsize) != 0)){
			warning("some ranges' widths are not a multiple of the selected
			 binsize, some bins will correspond to less than binsize basepairs")
		}

		pu <- pileup_core(gr, bampath, mapqual, binsize, shift, ss, 
		paired.end, paired.end.midpoint, paired.end.max.frag.length)
		new("CountSignals", counts=pu$counts, breaks=pu$breaks, ss=pu$ss)
	}
)

#' @export
setGeneric("bamCoverage", function(gr, bampath, ...) standardGeneric("bamCoverage"))
#' \code{bamCoverage}: for each base pair in the ranges, compute the number of reads
#' covering it.
#' @aliases bamCoverage
#' @rdname bamsignals-methods
#' @export
setMethod("bamCoverage", c("GenomicRanges", "character"), 
	function(gr, bampath, mapqual=0, paired.end=FALSE,
	paired.end.max.frag.length=1000, verbose=TRUE){
		if (verbose) printStupidSentence(bampath)
		
		pu <- coverage_core(gr, bampath, mapqual, paired.end, 
		paired.end.max.frag.length)
		new("CountSignals", counts=pu$counts, breaks=pu$breaks, ss=pu$ss)
	}
)

sentences <- c(
	"tHaT'S tHa fAStEsT pIlE-uP bAm iN tHe SoUth!!!",
	"yOu cAn'T pIlE-Up FaStEr!!!",
	"I'M gOnNa cHaSe'em and PiLe'em aLl up!!!",
	"fOr brOoMmHiLdA!!!",
	"tHe lEgEnD said, hE cOuLd PiLe uP fAsTeR thAn LiGht",
	"I gEt gOoSeBuMPs wHen I seE yOu pilEuPpiNg...")
printStupidSentence <- function(path){
	cat("Processing ", path, ": ", sample(sentences, 1), "\n", sep="")
}
