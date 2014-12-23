default: compile.log alltests

alltests: pileup.log count.log depth.log

compile.log: compile.R
	R --slave < compile.R 2>&1 | tee compile.log

pileup.log: pileup.R
	R --slave < pileup.R 2>&1 | tee pileup.log

count.log: count.R
	R --slave < count.R 2>&1 | tee count.log

depth.log: depth.R
	R --slave < depth.R 2>&1 | tee depth.log

clean:
	rm -f *.log data/Toy_SE.counts data/Toy_SE.profile data/Toy_SE.depth data/Toy_PE.counts data/Toy_PE.profile data/Toy_PE.depth

.PHONY: default alltest *.log
