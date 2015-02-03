# testthat script for artificial bam generation
#
# 2015-01-29

library(bamsignals)
context("bamsignals: Artifical data generation")

#
# Artificial data preparation
#
test_that("Artifical Data generation for single end data on chr1", {
		  binwidth <- 500
		  regionlength <- 50*binwidth
		  offs <- 3000000
		  readlen <- 50

		  #1. generate list of reads and sort by starting position
		  readpos <- matrix(rnbinom(2 * regionlength, 50, .99), nrow=2)
		  reads <- list()
		  for (i in 1:ncol(readpos)) {
			  mapq <- floor(runif(1, 0, 50))
			  if (readpos[1,i]!=0) { 
				  for (j in 1:readpos [1,i]) {
					  reads <- append(reads, list(c(0,           #FLAG
													i+offs,      #POS
													mapq         #MAPQ
													)))
				  }
			  }
			  if (readpos[2,i]!=0) {
				  for (j in 1:readpos[ 2,i]) {
					  reads <- append(reads, list(c(16,              #FLAG
													i+offs-readlen+1,#POS
													mapq             #MAPQ
													)))
				  }
			  }
		  }
		  out <- do.call("rbind", reads)
		  out <- out[ order(out[,2]),] #sort by POS

		  #2. write to SAM file
		  options("scipen"=100)#prevents exponential notation of integers
		  dir.create("data", showWarnings=F)
		  sam.file.se <- "data/bamsignals_SE_counts.sam"
		  writeLines("@HD\tVN:1.0\tSO:coordinate\n@SQ\tSN:chr1\tLN:249250621", sam.file.se)
		  write.table(x=cbind("*", out[,1], "chr1", out[,2], out[,3], paste(readlen,"M",sep=""), "*", "0", "0", "*", "*"),
					  file=sam.file.se,
					  sep="\t",
					  quote=F,
					  row.names=F,
					  col.names=F,
					  append=T
					  )

		  #3. SAM -> BAM, BAM->BAMi
		  bamsignals:::writeSamAsBamAndIndex(sam.file.se, "data/bamsignals_SE_counts.bam")

		  expect_true( file.exists( "data/bamsignals_SE_counts.sam" ), label="Artifical single end sam file successfully created")
		  expect_true( file.exists( "data/bamsignals_SE_counts.bam" ), label="Artifical single end bam file successfully created")
		  expect_true( file.exists( "data/bamsignals_SE_counts.bam.bai" ), label="Artifical single end bam file successfully indexed")
})

test_that("Artifical Data generation for paired end data on chr1", {
		  binwidth <- 500
		  regionlength <- 50*binwidth
		  offs <- 3000000
		  readlen <- 50

		  #1. generate list of reads and sort by starting position
		  readpos <- matrix(rnbinom(2 * regionlength, 50, .99), nrow=2)
		  reads <- list()
		  for (i in 1:ncol(readpos)) {
			  tlen <- floor(rnorm(1, 300, 50))
			  mapq <- floor(runif(1, 0, 50))
			  if (readpos[1,i]!=0) { 
				  for (j in 1:readpos [1,i]) {
					  reads <- append(reads, list(c(83,             #FLAG
													i+offs,         #POS
													mapq,           #MAPQ
													i+offs+tlen-50, #PNEXT
													tlen            #TLEN
													)))
					  #paired read
					  reads <- append(reads, list(c(163,            #FLAG
													i+offs+tlen-50, #POS
													mapq,           #MAPQ
													i+offs,         #PNEXT
													-tlen           #TLEN
													)))
				  }
			  }
			  if (readpos[2,i]!=0) {
				  for (j in 1:readpos[ 2,i]) {
					  reads <- append(reads, list(c(99,                 #FLAG
													i+offs-readlen+1,   #5'POS
													mapq,               #MAPQ
													i+offs+readlen-tlen,#PNEXT
													-tlen
													)))
					  #paired read
					  reads <- append(reads, list(c(147,                #FLAG
													i+offs+readlen-tlen,#POS
													mapq,               #MAPQ
													i+offs-readlen+1,   #PNEXT
													tlen                #TLEN
													)))
				  }
			  }
		  }
		  out <- do.call("rbind", reads)
		  out <- out[ order(out[,2]),] #sort by POS

		  #2. write to SAM file
		  options("scipen"=100)#prevents exponential notation of integers
		  dir.create("data", showWarnings=F)
		  sam.file.se <- "data/bamsignals_PE_counts.sam"
		  writeLines("@HD\tVN:1.0\tSO:coordinate\n@SQ\tSN:chr1\tLN:249250621", sam.file.se)
		  write.table(x=cbind("*", out[,1], "chr1", out[,2], out[,3], "50M", "=", out[,4], out[,5], "*", "*"),
					  file=sam.file.se,
					  sep="\t",
					  quote=F,
					  row.names=F,
					  col.names=F,
					  append=T
					  )

		  #3. SAM -> BAM, BAM->BAMi
		  bamsignals:::writeSamAsBamAndIndex(sam.file.se, "data/bamsignals_PE_counts.bam")

		  expect_true( file.exists( "data/bamsignals_PE_counts.sam" ), label="Artifical paired end sam file successfully created")
		  expect_true( file.exists( "data/bamsignals_PE_counts.bam" ), label="Artifical paired end bam file successfully created")
		  expect_true( file.exists( "data/bamsignals_PE_counts.bam.bai" ), label="Artifical paired end bam file successfully indexed")
})
