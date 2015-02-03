setClass( "CountSignals", 
	representation = representation( 
		counts = "integer",
		breaks = "integer",
		ss = "logical" )
)


setValidity( "CountSignals", 
	function(object) {
		x <- object
		if (length(x@ss)!=1 || is.na(x@ss)) return("invalid ss slot")
		len <- length(x@breaks)-1
		if (len < 0 || x@breaks[1]!=0 || x@breaks[len+1] != length(x@counts)) return("breaks first element must be 0 and the last must be length(x@counts)")
		if (len > 0){
			d <- diff(x@breaks)
			if (any(d) < 0) return("breaks must be increasing numbers")
			if (x@ss && any(d %% 2 == 1)) return("odd number of elements in a strand-specific signal...")
		}
		TRUE
	}
)

#' @export
setMethod("length", "CountSignals", function(x) length(x@breaks)-1)

#' @export
setMethod("width", "CountSignals", function(x) fastWidth(x))


#' @export
setMethod("[", "CountSignals", 
	function(x, i, j, ..., drop){
		if (!missing(j) || length(list(...)) > 0L)
			stop("invalid subsetting")
		
		if (length(i)==1 && drop){
			if (x@ss) return(getMatrix(x, i))
			else         return(getVector(x, i))
		} 
		
		largs <- getSubset(x, as.integer(i))
		new("CountSignals", counts=largs$counts, breaks=largs$breaks, ss=largs$ss)
	}
)

#' @export
setMethod("as.list", "CountSignals",
	function(x, ...){
		if (length(list(...)) > 0) stop("unrecognized options")
		asList(x)
	}
)

#' @export
setGeneric("alignSignals", function(x) standardGeneric("alignSignals"))
#' @export
setMethod("alignSignals", "CountSignals",
	function(x){
		lens <- diff(x@breaks)
		if (any(lens[1] != lens)) stop("all signals must have the same length")
		lastDim <- length(x)
		firstDims <- length(x@counts)/lastDim
		if (x@ss){
			ret <- x@counts
			dim(ret) <- c(2, firstDims/2, lastDim)
			dimnames(ret) <- list(c("sense", "antisense"), NULL, NULL)
			ret
		} else {
			ret <- x@counts
			dim(ret) <- c(firstDims, lastDim)
			ret
		}
	}
)


showCounts <- function(v){
	nprint <- min(length(v), 10)
	res <- paste0(as.character(v[1:nprint]), collapse=" ")
	if (nprint < length(v)) res <- paste(res, "...")
	res
}

showCountSignals <- function(x){
	cat("CountSignals object with ", length(x), ifelse(x@ss, " strand-specific", ""), 
	" signal", ifelse(length(x)!=1, "s", ""), fill=T, sep="")
	
	if (length(x) == 0) return()
	
	for (i in 1:min(5, length(x))){
		el <- x[i]
		npos <- ifelse(x@ss, ncol(el), length(el))
		cat("[",i, "] signal of width ", npos, fill=T, sep="")
		if (x@ss){
			cat("sense      ", showCounts(el[1,]), sep="", fill=T)
			cat("antisense  ", showCounts(el[2,]), sep="", fill=T)
		} else cat(showCounts(el), fill=T)
	}
	if (length(x)>5)
	cat("....", fill=T)
}

#' @export
setMethod("show", "CountSignals", function(object) showCountSignals(object))
