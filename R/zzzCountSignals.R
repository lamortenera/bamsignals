#' Container for count signals
#'
#' This s4 class is a tiny wrapper around a normal list (stored in the
#' \code{signals} slot) and it is the output of the methods in the
#' bamsignals package. 
#' Among other things the container provides an accessor method, 
#' that returns single signals as vectors and matrices, and the 
#' methods \code{as.list} and \code{alignSignals}, that convert the 
#' container to a list or an array/matrix respectively. A CountSignals
#' object is read-only, i.e. it cannot be modified.
#'
#' @param x A CountSignals object
#' @param i Index for subsetting. It can be a single index as well as
#' a vector of indices.
#' @param drop In case \code{i} is a vector of length 1, after subsetting, 
#' collapse the CountSignal object to a single signal or not.
#' @slot ss A single boolean value indicating whether all
#'     signals are strand-specific or not
#' @slot signals A list of integer vectors (if \code{ss==TRUE}) or of integer
#'     matrices, representing each signal
#' @aliases CountSignals
#' @seealso \code{\link{bamsignals-methods}} for the functions that produce 
#' this object
#' @return return values are described in the Methods section.
#' @example inst/examples/class_example.R
#' @export
setClass( "CountSignals", 
    representation = representation(signals = "list", ss = "logical" )
)


setValidity( "CountSignals", 
    function(object) {
        x <- object
        #check ss slot
        if (length(x@ss)!=1 || is.na(x@ss)) return("invalid ss slot")
        #check signals slot (must be written in C)
        if (!checkList(x@signals, x@ss)) return("invalid list")
        TRUE
    }
)

#' @describeIn CountSignals Number of contained signals
#' @aliases length
#' @export
setMethod("length", "CountSignals", function(x) length(x@signals))

#' @describeIn CountSignals Width of each signal. If the CountSignals
#' object \code{csig} is strand-specific then
#' \code{width(csig)[i] == ncol(csig[i])}, otherwise 
#' \code{width(csig)[i] = length(csig[i])}.
#' @aliases width
#' @export
setMethod("width", "CountSignals", function(x) {
    fastWidth(x@signals, x@ss)
})

#' @describeIn CountSignals Access single signals or subset the CountSignals 
#' object.
#' If \code{i} is a single index and \code{drop==TRUE} then the accessor returns
#' a single signal. If \code{x} is strand-specific then a single signal is a 
#' matrix with two rows, the first for the sense, the second for the antisense 
#' strand. Otherwise a signle signal is simply a vector of integers. If \code{i}
#' is a vector of length different than 1, then the acessor returns a subset of 
#' the CountSignals object. Invalid indices result into errors.
#' @aliases [,CountSignals,ANY-method
#' @export
setMethod("[", "CountSignals", 
    function(x, i, drop=TRUE){
        if (length(i)==1 && drop){
            return(x@signals[[i]])
        } else {
            subsigs <- x@signals[i]
            return(new("CountSignals", signals=subsigs, ss=x@ss))
        }
    }
)

#' @describeIn CountSignals Converts the container to a list \code{l} such that
#' \code{l[[i]]} is the i-th signal.
#' @aliases as.list
#' @export
setMethod("as.list", "CountSignals",
    #it should be using the generic defined in BiocGenerics
    function(x) x@signals
)

#' Converts the container to a list \code{l} such that
#' \code{l[[i]]} is the i-th signal.
#' @method as.list CountSignals
#' @param x A CountSignals object
#' @param ... not used
#' @return a list \code{l} such that
#' \code{l[[i]]} is \code{x[i]}.
#' @export
as.list.CountSignals <- function(x, ...) {x@signals}

#' @export
setGeneric("alignSignals", function(x) standardGeneric("alignSignals"))
#' @describeIn CountSignals Convert to a matrix or to an array. This is only
#' possible if all signals have the same width \code{w}. If the CountSignals
#' object \code{csig} is strand-specific, the result is an array of dimensions 
#' \code{[2, w, length(csig)]}, otherwise it will be a matrix of dimensions
#' \code{[w, length(csig)]}.
#' @aliases alignSignals
#' @export
setMethod("alignSignals", "CountSignals",
    function(x){
        ws <- width(x)
        if (any(ws != ws[1])) stop("all signals must have the same length")
        simplify2array(x@signals)
    }
)


showCounts <- function(v){
    nprint <- min(length(v), 10)
    res <- paste0(as.character(v[1:nprint]), collapse=" ")
    if (nprint < length(v)) res <- paste(res, "...")
    res
}

showCountSignals <- function(x){
    cat("CountSignals object with ", length(x), 
    ifelse(x@ss, " strand-specific", ""), 
    " signal", ifelse(length(x)!=1, "s", ""), fill=TRUE, sep="")
    
    if (length(x) == 0) return()
    
    for (i in 1:min(5, length(x))){
        el <- x[i]
        npos <- ifelse(x@ss, ncol(el), length(el))
        cat("[",i, "] signal of width ", npos, fill=TRUE, sep="")
        if (x@ss){
            cat("sense      ", showCounts(el[1,]), sep="", fill=TRUE)
            cat("antisense  ", showCounts(el[2,]), sep="", fill=TRUE)
        } else cat(showCounts(el), fill=TRUE)
    }
    if (length(x)>5)
    cat("....", fill=TRUE)
}

setMethod("show", "CountSignals", function(object) showCountSignals(object))
