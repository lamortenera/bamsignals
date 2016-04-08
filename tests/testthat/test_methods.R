context("bamsignals methods")
library(GenomicRanges)
library(Rsamtools)
source("utils.R")

#toy data
reads <- get(load("randomReads.RData"))
bampath <- system.file("extdata", "randomBam.bam", package="bamsignals")

#generate regions
nRef <- c("chr1", "chr2", "chr3")  #number of chromosomes
nRegions <- 20    #number of regions
lastStart <- 1e3 #last start of a read pair (i.e. chromosome length more or less)
avgRegLen <- 200  #average region length
regions <- GRanges(
  seqnames=sample(nRef, nRegions, replace=TRUE),
  strand=sample(c("+", "-"), nRegions, replace=TRUE),
  ranges=IRanges(
    start=sample( lastStart, nRegions, replace=TRUE), 
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
      for (ss in c(FALSE, TRUE)){
        for (pe in c("ignore", "filter", "midpoint")){
          for (tlenFilter in list(NULL, c(50, 200))){
              expect_equal(label=paste0("bamCount{", 
                                        paste("shift", shift, "mapq", mapq, "ss", ss, "pe", pe, 
                                              "tlenFilter", paste0(tlenFilter, collapse=","), 
                                              sep="="),
                                        "}"),
                           countR(reads, regions, ss=ss, shift=shift, paired.end=pe, mapqual=mapq,
                                  tlenFilter=tlenFilter), 
                           bamCount(bampath, regions, ss=ss, shift=shift, paired.end=pe, mapqual=mapq, 
                                    tlenFilter=tlenFilter, verbose=FALSE))
  }  }  }  }  }
})

test_that("bamProfile function", {
   for (shift in c(0, 100)){
     for (mapq in c(0, 100)){
      for (ss in c(FALSE, TRUE)){
        for (pe in c("ignore", "filter", "midpoint")){
          for (tlenFilter in list(NULL, c(50, 200))){
              expect_equal(label=paste0("bamProfile{", 
                                        paste("shift", shift, "mapq", mapq, "ss", ss, "pe", pe, 
                                              "tlenFilter", paste0(tlenFilter, collapse=","), 
                                              sep="="),
                                        "}"),
                           profileR(reads, regions, ss=ss, shift=shift, paired.end=pe, mapqual=mapq,
                                    tlenFilter=tlenFilter), 
                           as.list(bamProfile(bampath, regions, ss=ss, shift=shift, paired.end=pe, mapqual=mapq, 
                                    tlenFilter=tlenFilter, verbose=FALSE)))
  }  }  }  }  }
})


test_that("bamCoverage function", {
  for (mapq in c(0, 100)){
    for (pe in c("ignore", "extend")){
      for (tlenFilter in list(NULL, c(50, 200))){
          expect_equal(label=paste0("bamCoverage{", 
                                    paste("mapq", mapq, "pe", pe, 
                                          "tlenFilter", paste0(tlenFilter, collapse=","), 
                                          sep="="),
                                    "}"),
                       coverageR(reads, regions, paired.end=pe, mapqual=mapq, tlenFilter=tlenFilter), 
                       as.list(bamCoverage(bampath, regions, paired.end=pe, mapqual=mapq, tlenFilter=tlenFilter, 
                                verbose=FALSE)))
  }  }  }
})

test_that("filtering on SAMFLAGS", {
  strand(regions) <- "+"#need to set this for below test
  for (shift in c(0, 100)){
    for (mapq in c(0, 100)){
      for (pe in c("ignore", "filter", "midpoint")){
        for (tlenFilter in list(NULL, c(50, 200))){
          expect_equal(label=paste0("filterFlag Test bamCount{", 
                                    paste("shift", shift, "mapq", mapq, "pe", pe, 
                                          "tlenFilter", paste0(tlenFilter, collapse=","), 
                                          sep="="),
                                    "}"),
                       countR(reads, regions, ss=T, shift=shift, paired.end=pe, mapqual=mapq,
                          tlenFilter=tlenFilter)[1,], 
                       bamCount(bampath, regions, ss=F, shift=shift, paired.end=pe, mapqual=mapq, 
                                tlenFilter=tlenFilter, filteredFlag=16, verbose=FALSE))
        }
      }
    }
  }
})

