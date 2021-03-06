##CONVERT READS FROM A DATA-FRAME-LIKE FORMAT TO BAM

#to print integers as characters without weird stuff
int2chr <- function(nums){
    format(nums, scientific=FALSE, trim=TRUE, justify="none")
}

isPEData <- function(reads){
    all(c("isize", "read1", "pnext") %in% names(reads))
}

#converts reads in R format to an indexed bamfile
#the format for the reads is similar to the one used by
#scanBam in Rsamtools, except that only the fields
#qname, rname*, strand*, pos*, qwidth*, mapq* are used 
#('*' means required)
#and for paired end reads also the fields read1 and isize.
#The field read1 is not in the Rsamtools format, but we need
#it to compute the flag. It means "is this read here the read 1
#in the pair?", and it defines the orientation of the pair.
#another difference, the reads object must inherit from "data.frame"
#another difference, there is no flag, every read is supposed to
#be mapped properly, as well as its paired end (if it is there)
#refs is a list with the fields 'refname' and 'reflen'
readsToBam <- function(reads, bampath, refs=NULL){
    #checking and completing fields
    if (!inherits(reads,"data.frame")) stop("reads must be a data.frame")
    #check that all required fields are there
    reqFields <- c("rname", "strand", "pos", "qwidth", "mapq")
    if (any(! reqFields %in% names(reads))) stop("missing required fields")
    #check that all fields have the right length
    nreads <- nrow(reads)
    #check that the format for the strand is correct
    if (any(!reads$strand %in% c("+", "-"))) stop("strand can only be '-' or '+'")
    #deal with paired end reads
    pe <- isPEData(reads)
    if (!pe){
        if (any(c("isize", "read1", "pnext") %in% names(reads))){ 
            warning("incomplete paired-end specification, pair information will be ignored")
        }
        reads$isize <- rep(0, nreads)
        reads$pnext <- reads$isize
        reads$read1 <- rep(TRUE, nreads)
    } else {
        if (!is.logical(reads$read1)) stop("the 'read1' field must be of type logical")
        if (any(reads$isize==0)) stop("an 'isize' of 0 is not allowed")
        if (any((reads$strand=="+")!=(reads$isize>0))) {
            warning("the sign of 'isize' should be the same as the strand of the read
            for properly mapped read pairs. Fixing that.")
            reads$isize <- abs(reads$isize)*(2*as.integer(reads$strand=="+") - 1)
        }
    }
    #complete missing fields
    if (!"qname" %in% names(reads)) reads$qname <- 1:nreads
    reads$mapq[is.na(reads$mapq)] <- 255
    
    
    #compute the flag
    #reads$flag <- 0x10*(reads$strand=="-")
    #if (pe) reads$flag <- reads$flag + 0x1 + 0x2 + 0x80 + (0x40-0x80)*reads$read1 + 0x20*(reads$strand=="+")
    
    
    #deal with the references
    needRefs <- unique(reads$rname)
    if (!is.null(refs)){
        if (any(! c("refname", "reflen") %in% names(regs))) stop("missing required fields in 'refs'")
        if (length(refs$refname) != length(refs$reflen)) stop("all fields must have the same length in 'refs'")
        if (any(! needRefs %in% refs$refname)) stop("missing some chromosome names")
    } else {
        refs <- list(
            refname=needRefs,
            reflen=sapply(needRefs, function(chr){
                v <- reads$rname == chr
                max(reads$pos[v] + reads$qwidth[v]) + 1
            }))
    }
    #sort the references
    refs <- data.frame(refs)
    o <- order(refs$refname)
    refs <- refs[o,]
    
    #build the final data frame
    df <- data.frame(
        qname=reads$qname,
        flag=reads$flag,
        rname=reads$rname,
        pos=reads$pos,
        mapq=reads$mapq,
        cigar=paste0(reads$qwidth, "M"),
        rnext="=",
        pnext=reads$pnext,
        tlen=reads$isize,
        seq="*",
        qual="*")
    
    #sort the reads
    o <- order(df$rname, df$pos)
    df <- df[o,]
    
    #avoid scientific notation
    df$pos <- int2chr(df$pos)
    df$qname <- int2chr(df$qname)
    
    #create temporary samfile
    tmpsam <- paste0(tempfile(), ".sam")
    #write header
    header <- c(
        "@HD\tVN:1.0\tSO:coordinate",
        paste0("@SQ\tSN:", refs$refname, "\tLN:", int2chr(refs$reflen)))
    writeLines(header, tmpsam)

    #write reads
    write.table(x=df, file=tmpsam, sep="\t", quote=F, row.names=F, col.names=F, append=T)
    
    #write to bam and index
    bamsignals:::writeSamAsBamAndIndex(tmpsam, bampath)
    
    #remove temporary sam file
    file.remove(tmpsam)
}

generateToyData <- function() {
  nRef <- c("chr1", "chr2", "chr3")  #number of chromosomes
  lastStart <- 1e4 #last start of a read pair (i.e. chromosome length more or less)
  nPairs <- 5e4 #number of read pairs
  avgRLen <- 30 #average read length
  avgFLen <- 150 #average fragment length or template length
  nRemove <- 1e3 #remove random reads from the pairs
  bampath <- paste0(tempfile(), ".bam") #path where to save the bam file

  #generate reads
  posReads <- data.frame(
    rname=sample(nRef, nPairs, replace=TRUE),      #chromosome name
    pos=sample(sample(100:lastStart, 1e3, replace=T), nPairs, replace=T) 
        + as.integer(rnorm(nPairs,0,30)), #start of the read
    qwidth=1+rpois(nPairs, lambda=(avgRLen-1)),    #length of the read
    strand=rep("+", nPairs),                       #strand
    isize=1+rnbinom(nPairs, mu=(avgFLen-1), size=10),#template (or fragment) length
    read1=(sample(2, nPairs, replace=TRUE)==1),    #is it the first read in the pair?
    mapq=rnbinom(nPairs,size=20,prob=.5),          #mapping quality
    flag=0x0)                                      #flag 0 for positive read

  negReads <- data.frame(
    rname=posReads$rname,
    qwidth=1+rpois(nPairs, lambda=(avgRLen-1)),
    strand=rep("-", nPairs),
    isize=-posReads$isize,
    read1=!posReads$read1,
    mapq=rnbinom(nPairs,size=20,prob=.5),
    flag=0x10)
  negReads$pos <- posReads$pos + posReads$isize - negReads$qwidth
  posReads$pnext <- negReads$pos
  negReads$pnext <- posReads$pos

  #mark duplicates 
  dups <- which(duplicated(cbind(posReads[,1:2], negReads[,1:2])))
  posReads$flag[dups] <- posReads$flag[dups] + 0x400
  negReads$flag[dups] <- negReads$flag[dups] + 0x400

  ##merge pos and neg reads and remove some of them
  reads <- rbind(posReads, negReads)
  reads <- reads[-sample(nrow(reads), nRemove),]

  #set paired end flags
  reads$flag <- reads$flag + 0x1 + 0x2 + 0x80 + (0x40-0x80)*reads$read1 + 0x20*(reads$strand=="+")

  #write to a bam file
  bampath <- "inst/extdata/randomBam.bam"
  readsToBam(reads, bampath)
  save(reads, file="inst/extdata/randomReads.RData")
}


##REWRITE BAMSIGNALS FUNCTIONS IN R
##READS ARE IN A DATA-FRAME-LIKE FORMAT

#convert the reads to a GRanges object
df2gr <- function(reads, paired.end=FALSE, paired.end.midpoint=FALSE, shift=0, mapqual=0,
                  tlenFilter=NULL){
  if (!paired.end %in% c("ignore", "filter", "midpoint", "extend"))
    stop("invalid paired.end option")

  #filter based on mapqual
  reads <- reads[reads$mapq>=mapqual,]

  #if paired.end, discard reads that are not the first read in the pair
  if (paired.end != "ignore") {
    reads <- reads[reads$read1,]
    #filter reads on supplied maximum and minimum fragment lengths
    if (is.null(tlenFilter)) { #defaults to c(0,1000)
      reads <- reads[abs(reads$isize) >= 0 & abs(reads$isize) <= 1000, ]
    } else {
      reads <- reads[abs(reads$isize) >= tlenFilter[1] & abs(reads$isize) <= tlenFilter[2], ]
    }
  }

  #do as if the pair was a single read
  if (paired.end %in% c("extend", "midpoint")){
    isNeg <- reads$strand=="-"
    #shift back negative reads
    reads[isNeg,]$pos <- reads[isNeg,]$pos - abs(reads[isNeg,]$isize) + reads[isNeg,]$qwidth
    #extend width to the template width
    reads$qwidth <- abs(reads$isize)
  }

  #convert to GRanges
  gr <- GRanges(seqnames=reads$rname, strand=reads$strand, IRanges(start=reads$pos, width=reads$qwidth))

  #if paired.end.midpoint, consider only the midpoint
  if (paired.end == "midpoint"){
    mids <- (start(gr) + end(gr))/2
    #take care of the rounding
    signedMids <- mids*(2*as.integer(strand(gr)=="+")-1)
    mids <- abs(ceiling(signedMids))

    start(gr) <- mids
    end(gr) <- mids
  }

  #consider shift
  shifts <- rep(shift, length(gr))
  shifts[as.logical(strand(gr)=="-")] <- -shift

  GenomicRanges::shift(gr, shifts)
}


countR <- function(reads, genes, ss=FALSE, ...){
    reads <- df2gr(reads, ...)
    
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

profileR <- function(reads, genes, ss=FALSE, ...){
    reads <- df2gr(reads, ...)
    
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

coverageR <- function(reads, genes, ...){
    reads <- df2gr(reads, ...)
    
    isNegGene <- as.logical(strand(genes)=="-")
    seqNames <- as.character(seqnames(genes))
    
    #use the already implemented coverage function
    #but first make sure that seqinfo of the reads is set properly
    rng <- range(c(genes, reads), ignore.strand=TRUE)
    seqlevels(reads) <- as.character(seqnames(rng))
    seqlengths(reads) <- end(rng)
    
    cvrg <- coverage(reads)
    lapply(seq_along(genes), function(g){
        sig <- cvrg[[seqNames[g]]][start(genes)[g]:end(genes)[g]]
        #take into account strand of the region
        if (isNegGene[g]) sig <- rev(sig)
        as.integer(sig)
    })
}

