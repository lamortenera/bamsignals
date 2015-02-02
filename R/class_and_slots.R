setClass( "ConcatCounts", 
	representation = representation( 
		counts = "integer",
		starts = "integer",
		ends = "integer",
		ss = "logical" )
)

setValidity( "ConcatCounts", function( x ) {
	if (length(x@ss)!=1 || is.na(x@ss)) return("invalid ss slot")
	if (length(x@starts) != length(x@ends)) return("starts and ends don't match")
	if (any(x@ends - x@starts < 0)) return("negative ranges are not allowed")
	len <- length(x@starts)
	if (len > 0){
		if (any(x@starts[2:len] != x@ends[1:(len-1)]+1)) return("starts and ends don't match")
		if (x@ss && any((x@ends-x@starts) %% 2 == 0)) return("odd number of elements in a strand-specific signal...")
		if (x@starts[1] != 1 || x@ends[len] != length(x@counts)) return("invalid starts or ends slot")
	}
	TRUE
})

setMethod("length", "ConcatCounts", function(x) length(x@starts))

setMethod(
	f = "[", 
	signature = "ConcatCounts",
	definition = function(x, i, j, ..., drop){
		if (!missing(j) || length(list(...)) > 0L)
			stop("invalid subsetting")
		
		#needs to be done in C
		if (any(i > length(x))) stop("index out of bounds")
		
		nlen <- length(i)
		nstarts <- integer(nlen)
		nends <- integer(nlen)
		#figure out starts and ends in the new vector
		#and total amount of memory to allocate
		acc <- 0L
		for (ii in seq_along(i)){
			idx <- i[ii]
			nstarts[ii] <- acc + 1L
			acc <- acc + x@ends[idx] - x@starts[idx] + 1L
			nends[ii] <- acc
		}
		
		#copy the subsets in a new vector
		ncounts <- integer(acc)
		for (ii in seq_along(i)){
			idx <- i[ii]
			ncounts[nstarts[ii]:nends[ii]] <- x@counts[x@starts[idx]:x@ends[idx]]
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
		new("ConcatCounts", counts=ncounts, starts=nstarts, ends=nends, ss=x@ss)
	}
)


setMethod("as.list", "ConcatCounts",
	function(x, ...){
		if (length(list(...)) > 0) stop("unrecognized options")
		l <- list()
		for (i in seq_along(x)) l[[i]] <- x[i]
		l
	}
)

setGeneric("alignSignals", function(x) standardGeneric("alignSignals"))
setMethod("alignSignals", "ConcatCounts",
	function(x){
		lens <- x@ends - x@starts
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

setMethod("show", "ConcatCounts", function(x){
	cat("ConcatCounts object with ", length(x), ifelse(x@ss, " strand-specific", ""), 
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
