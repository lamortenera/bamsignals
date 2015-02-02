setClass( "CountSignals", 
	representation = representation( 
		counts = "integer",
		breaks = "integer",
		ss = "logical" )
)

setValidity( "CountSignals", function( x ) {
	if (length(x@ss)!=1 || is.na(x@ss)) return("invalid ss slot")
	len <- length(x@breaks)-1
	if (len < 0 || x@breaks[1]!=0 || x@breaks[len+1] != length(x@counts)) return("breaks first element must be 0 and the last must be length(x@counts)")
	if (len > 0){
		d <- diff(x@breaks)
		if (any(d) < 0) return("breaks must be increasing numbers")
		if (x@ss && any(d %% 2 == 1)) return("odd number of elements in a strand-specific signal...")
	}
	TRUE
})

setMethod("length", "CountSignals", function(x) length(x@breaks)-1)

setMethod(
	f = "[", 
	signature = "CountSignals",
	definition = function(x, i, j, ..., drop){
		if (!missing(j) || length(list(...)) > 0L)
			stop("invalid subsetting")
		
		#needs to be done in C
		if (any(i > length(x) || i <= 0)) stop("index out of bounds")
		
		nlen <- length(i)
		nbreaks <- integer(nlen+1)
		#figure out starts and ends in the new vector
		#and total amount of memory to allocate
		acc <- 0L; nbreaks[1] <- acc
		for (ii in seq_along(i)){
			idx <- i[ii]
			acc <- acc + x@breaks[idx+1] - x@breaks[idx]
			nbreaks[ii+1] <- acc
		}
		
		#copy the subsets in a new vector
		ncounts <- integer(acc)
		for (ii in seq_along(i)){
			idx <- i[ii]
			ncounts[(nbreaks[ii]+1):nbreaks[ii+1]] <- x@counts[(x@breaks[idx]+1):x@breaks[idx+1]]
		}
		
		#deal with the drop option
		if (length(i)==1 && drop){
			if (x@ss){
				dim(ncounts) <- c(2, length(ncounts)/2)
				rownames(ncounts) <- c("sense", "antisense")
			}
			return(ncounts)
		}
		
		#return new s4 object
		new("CountSignals", counts=ncounts, breaks=nbreaks, ss=x@ss)
	}
)


setMethod("as.list", "CountSignals",
	function(x, ...){
		if (length(list(...)) > 0) stop("unrecognized options")
		l <- list()
		for (i in seq_along(x)) l[[i]] <- x[i]
		l
	}
)

setGeneric("alignSignals", function(x) standardGeneric("alignSignals"))
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
})


showCounts <- function(v){
	nprint <- min(length(v), 10)
	res <- paste0(as.character(v[1:nprint]), collapse=" ")
	if (nprint < length(v)) res <- paste(res, "...")
	res
}

setMethod("show", "CountSignals", function(x){
	cat("CountSignals object with ", length(x), ifelse(x@ss, " strand-specific", ""), 
	" signal", ifelse(length(x)!=1, "s", ""), fill=T, sep="")
	
	if (length(x) == 0) return()
	
	for (i in 1:min(5, length(x))){
		el <- x[i]
		npos <- ifelse(x@ss, ncol(el), length(el))
		cat("[",i, "] signal of length ", npos, fill=T, sep="")
		if (x@ss){
			cat("sense      ", showCounts(el[1,]), sep="", fill=T)
			cat("antisense  ", showCounts(el[2,]), sep="", fill=T)
		} else cat(showCounts(el), fill=T)
	}
	if (length(x)>5)
	cat("....", fill=T)
})
