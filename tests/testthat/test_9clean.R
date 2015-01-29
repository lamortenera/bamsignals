# testthat script for artificial data cleanup
#
# 2015-01-29

library(bamsignals)
context("bamsignals test case artifical data clean up")

#
# Artificial data preparation
#
test_that("Artifical Data generation is done correctly", {
		  expect_true( file.remove( "data/Artifical_SE_counts.sam" ), label="data/Artifical_SE_counts.sam successfully removed")
		  expect_true( file.remove( "data/Artifical_SE_counts.bam" ), label="data/Artifical_SE_counts.bam file successfully removed")
		  expect_true( file.remove( "data/Artifical_SE_counts.bam.bai" ), label="data/Artifical_SE_counts.bam.bai successfully removed")

		  expect_true( file.remove( "data/Artifical_PE_counts.sam" ), label="data/Artifical_PE_counts.sam successfully removed")
		  expect_true( file.remove( "data/Artifical_PE_counts.bam" ), label="data/Artifical_PE_counts.bam file successfully removed")
		  expect_true( file.remove( "data/Artifical_PE_counts.bam.bai" ), label="data/Artifical_PE_counts.bam.bai successfully removed")
})

