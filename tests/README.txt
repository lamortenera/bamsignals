Tests for bamsignals
--------------------------------------------------------------------------------------------------
2014-03-31

This folder contains various tests for the bamsignals R package. Therefore it uses bedtools, 
different R scripts and data contained in ./data. Results are compared to outputs of
cmbr.R::countBamInGRangesFast and cmbr::coverageBamInGRangesFast.

compilation.R:
  compilation
  checks for system availability of R packages used in tests
  
deploy.R: compilation.R
  bamsignals::count
  bamsignals::depth
  bamsignals::pileup

#EOF
