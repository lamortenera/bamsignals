# testthat script for artificial bam generation
#
# 2015-01-29

library(bamsignals)
context("bamsignals test case artifical data generation")

#
# Artificial data preparation
#
test_that("Artifical Data generation is done correctly", {
#TODO generate also PE bam in here
		  binwidth <- 500
		  regionlength <- 20*binwidth

		  #1. generate list of reads and sort by starting position
		  readpos <- matrix( rnbinom( 2 * regionlength, 50, .99), nrow=2 )
		  reads <- list()
		  for (i in 1:ncol(readpos)) {
			  if (readpos[1,i]!=0) { 
				  for (j in 1:readpos [1,i]) {
					  reads <- append(reads, list(c(i, min(i+50-1, regionlength), 0, floor(runif(1, 0, 50)))))
				  }
			  }
			  if (readpos[2,i]!=0) {
				  for (j in 1:readpos[ 2,i]) {
					  reads <- append(reads, list(c(max(i-50+1,0), i, 16, floor(runif(1, 0, 50)))))
				  }
			  }
		  }
		  out <- do.call("rbind", reads)
		  out <- out[ order(out[,1]),]

		  #2. write to SAM file
		  dir.create("data", showWarnings=F)
		  sam.file.se <- "data/bamsignals_SE_counts.sam"
		  writeLines("@HD\tVN:1.0\tSO:coordinate\n@SQ\tSN:chr1\tLN:249250621", sam.file.se)
		  write.table(x=cbind("*", out[,3], "chr1", out[,1], out[,4], "50M", "*", "0", "0", "*", "*"),
					  file=sam.file.se,
					  sep="\t",
					  quote=F,
					  row.names=F,
					  col.names=F,
					  append=T
					  )

		  #3. SAM -> BAM, BAM->BAMi
		  bamsignals:::writeSamAsBamAndIndex(sam.file.se, "data/bamsignals_SE_counts.bam")

		  expect_true( file.exists( "data/samsignals_SE_counts.sam" ), label="Artifical sam file successfully created")
		  expect_true( file.exists( "data/bamsignals_SE_counts.bam" ), label="Artifical bam file successfully created")
		  expect_true( file.exists( "data/bamsignals_SE_counts.bam.bai" ), label="Artifical bam file successfully indexed")

		  expect_true( file.exists( "data/samsignals_PE_counts.sam" ), label="Artifical sam file successfully created")
		  expect_true( file.exists( "data/bamsignals_PE_counts.bam" ), label="Artifical bam file successfully created")
		  expect_true( file.exists( "data/bamsignals_PE_counts.bam.bai" ), label="Artifical bam file successfully indexed")
})
