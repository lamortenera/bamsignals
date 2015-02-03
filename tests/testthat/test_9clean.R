# testthat script for artificial data cleanup
#
# 2015-01-29

library(bamsignals)
context("bamsignals: Artifical data clean up")

#
# Artificial data preparation
#
test_that("Artifical Data is removed correctly", {
		  expect_true( file.remove( "data/bamsignals_SE_counts.sam" ), label="data/bamsignals_SE_counts.sam successfully removed")
		  expect_true( file.remove( "data/bamsignals_SE_counts.bam" ), label="data/bamsignals_SE_counts.bam file successfully removed")
		  expect_true( file.remove( "data/bamsignals_SE_counts.bam.bai" ), label="data/bamsignals_SE_counts.bam.bai successfully removed")

		  expect_true( file.remove( "data/bamsignals_PE_counts.sam" ), label="data/bamsignals_PE_counts.sam successfully removed")
		  expect_true( file.remove( "data/bamsignals_PE_counts.bam" ), label="data/bamsignals_PE_counts.bam file successfully removed")
		  expect_true( file.remove( "data/bamsignals_PE_counts.bam.bai" ), label="data/bamsignals_PE_counts.bam.bai successfully removed")
})

