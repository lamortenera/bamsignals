context("CountSignals class and methods")

getSig <- function(i, ss){
	nums <- as.integer( ((i-1)*4 + 1) : (i*4) )
	if (ss) matrix(nums, nrow=2, dimnames=list(c("sense", "antisense"), NULL))
	else nums
}

runs <- function(expr){
	res <- try(force(expr), TRUE)
	msg1 <- "code did not generate an error"
	msg2 <- "code generated an error"
	expectation("error", msg2, msg1)
}

expect_runs <- function(expr){
	expect_that(expr, runs, label=testthat:::find_expr("expr"))
}

test_that("Test CountSignals class and methods", {
	for (ss in c(TRUE,FALSE)){
		n <- new("CountSignals", signals=lapply(1:4, getSig, ss=ss), ss=ss)

		#see if method length works
		expect_equal(length(n), 4)

		#see if width works
		expect_equal(width(n), rep(ifelse(ss, 2, 4), 4))
		
		#see if method show works
		expect_runs(capture.output(show(n)))

		#subsetting with zero elements
		expect_runs(n[c()])
		
		#negative indices
		expect_error(n[-1])
		
		#too large indices
		expect_error(n[5])
		
		#check if accessor works
		for (i in 1:length(n)){
			expect_equal(n[i], getSig(i, n@ss))
		}

		#check subsetting
		sn <- n[c(2,4)]
		expect_equal(sn[1], n[2])
		expect_equal(sn[2], n[4])
		
		#check list
		l <- as.list(n)
		for (i in seq_along(l)) expect_equal(l[[i]], n[i])

		#check alignSignals
		a <- alignSignals(n)
		for (i in seq_along(l)) {
			if (ss) sa <- a[,,i]
			else sa <- a[,i]
			expect_equal(sa, n[i])
		}
	}}
)

