#' Efficient counting of reads in a bam file
#'
#' Functions and classes for extracting signals from a bam file. Differently 
#' than the other packages, this package is not supposed to import reads in R.
#' All the read-processing is done in C/C++ and the only output are read counts. 
#'
#' @name bamsignals
#' @import Rcpp
#' @import IRanges
#' @import GenomicRanges
#' @import zlibbioc
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
#' @param bampath path to the bam file storing the read. The file must be
#' indexed. 
#' @param gr GenomicRanges object used to specify the regions. If a range is on
#' the negative strand the profile will be reverse-complemented.
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
#' @param paired.end a character string indicating how to handle paired-end 
#' reads. If \code{paired.end!="ignore"} then only first reads in proper mapped
#' pairs will be considered (SAMFLAG 66, i.e. in the flag of the read, the bits 
#' in the mask 66 must be all ones). \cr
#' If \code{paired.end=="midpoint"} then the midpoint of a filtered fragment is 
#' considered, where \code{mid = fragment_start + int(abs(tlen)/2)}, and where 
#' tlen is the template length stored in the bam file. For even tlen, the given 
#' midpoint will be moved of 0.5 basepairs in the 3' direction (bamCount, 
#' bamProfile). \cr
#' If \code{paired.end=="extend"} then the whole fragment is treated 
#' as a single read (bamCoverage).
#' @param tlen.filter A filter on fragment length as estimated from alignment 
#' in paired end experiments (TLEN). If set to \code{c(min,max)} only reads are 
#' considered where \code{min <= TLEN <= max}. If \code{paired.end=="ignore"}, 
#' this argument is set to \code{NULL} and no filtering is done. If 
#' \code{paired.end!="ignore"}, this argument defaults to \code{c(0,1000)}.
#' @param filteredFlag Filter out reads with a certain flag set, e.g. "1024" to
#' filter out PCR or optical duplicates.
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
#' coordinate. This can be changed using the \code{shift} parameter. \cr
#'  
#' @seealso \code{\link{CountSignals}} for handling the return value of 
#' \code{bamProfile} and \code{bamCoverage}
#' @name bamsignals-methods
#' @example inst/examples/methods_example.R
NULL

#helper functions to set the right flag
flagMask <- function(paired.end){
  #considers only the first read in a properly mapped pair
  if (paired.end != "ignore") return(66L)
  #considers all reads
  0L
}

#helper functions to set the right flag
tlenFilter <- function(tlen.filter, paired.end){
  #No filtering if paired end is ignored
  if (paired.end == "ignore") return(integer())
  #Default filter if tlen.filter is unset
  if (is.null(tlen.filter)) return(c(0,1000))
  #Check supplied tlen.filter
  if (length(tlen.filter) != 2 || tlen.filter[1] < 0 ||
      tlen.filter[2] < 0) {
    stop("tlen.filter must be NULL or vector of 2 positive integers")
  } 
  if (tlen.filter[1] > tlen.filter[2]) {
    stop("tlen.filter[1] must be smaller or equal to tlen.filter[2]")
  }
  tlen.filter
}

#' @export
setGeneric("bamCount", function(bampath, gr, ...) standardGeneric("bamCount"))
#' \code{bamCount}: for each range, count the reads whose 5' end map in it.
#' @aliases bamCount
#' @rdname bamsignals-methods
#' @export
setMethod("bamCount", c("character", "GenomicRanges"), 
    function(bampath, gr, mapqual=0, shift=0, ss=FALSE, 
    paired.end=c("ignore", "filter", "midpoint"), 
    tlen.filter=NULL, filteredFlag=-1, verbose=TRUE){
        if (verbose) printStupidSentence(bampath)
        
        pe <- match.arg(paired.end)
        
        bampath <- path.expand(bampath)
        pu <- pileup_core(bampath, gr, tlenFilter(tlen.filter, pe), mapqual, 
        -1, shift, ss, flagMask(pe), filteredFlag, (pe=="midpoint"))
        
        pu[[1]]
    }
)


#' @export
setGeneric("bamProfile", function(bampath, gr, ...) standardGeneric("bamProfile"))
#' \code{bamProfile}: for each base pair in the ranges, compute the number of reads
#' whose 5' end maps there.
#' @aliases bamProfile
#' @rdname bamsignals-methods
#' @export
setMethod("bamProfile", c("character", "GenomicRanges"), 
    function(bampath, gr, binsize=1, mapqual=0, shift=0, ss=FALSE, 
    paired.end=c("ignore", "filter", "midpoint"),
    tlen.filter=NULL, filteredFlag=-1, verbose=TRUE){
        if (verbose) printStupidSentence(bampath)
        
        if (binsize < 1){
            stop("provide a binsize greater or equal to 1")
        } else if (binsize > 1 && any((width(gr) %% binsize) != 0)){
            warning("some ranges' widths are not a multiple of the selected
             binsize, some bins will correspond to less than binsize basepairs")
        }

        pe <- match.arg(paired.end)
        
        bampath <- path.expand(bampath)
        pu <- pileup_core(bampath, gr, tlenFilter(tlen.filter, pe), mapqual, binsize, 
        shift, ss, flagMask(pe), filteredFlag, (pe=="midpoint"))
        
         new("CountSignals", signals=pu, ss=ss)
    } 
)

#' @export
setGeneric("bamCoverage", function(bampath, gr, ...) standardGeneric("bamCoverage"))
#' \code{bamCoverage}: for each base pair in the ranges, compute the number of reads
#' covering it.
#' @aliases bamCoverage
#' @rdname bamsignals-methods
#' @export
setMethod("bamCoverage", c("character", "GenomicRanges"), 
    function(bampath, gr, mapqual=0, paired.end=c("ignore", "extend"),
    tlen.filter=NULL, filteredFlag=-1, verbose=TRUE){
        if (verbose) printStupidSentence(bampath)
        
        pe <- match.arg(paired.end)
        
        bampath <- path.expand(bampath)
        pu <- coverage_core(bampath, gr, tlenFilter(tlen.filter, pe), mapqual, 
        flagMask(pe), filteredFlag, (pe=="extend"))
        
        new("CountSignals", signals=pu, ss=FALSE)
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
    message("Processing ", path, ": ", sample(sentences, 1))
}
