#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::IntegerMatrix getMatrix(Rcpp::RObject csig, int idx){
	if (not csig.inherits("CountSignals"))
		Rcpp::stop("expecting a CountSignals object");
	
	//accessing fields
	Rcpp::IntegerVector counts = as<IntegerVector>(gr.slot("counts");
	Rcpp::IntegerVector breaks = as<IntegerVector>(gr.slot("breaks");
	bool ss = as<bool>(gr.slot("ss");
	
	//checking preconditions
	if (not ss) Rcpp::stop("expecting a strand-specific CountSignals object");
	if (idx < 0 || idx >= breaks.length()-1) Rcpp::stop("invalid index");
	
	//creating matrix
	int len = breaks[idx+1] - breaks[idx];
	int ncol = len/2;
	Rcpp::IntegerMatrix mat(2, ncol);
	memcpy(mat.begin(), counts.begin() + breaks[idx], sizeof(int)*len)
	
	return mat;
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix getVector(Rcpp::RObject csig, int idx){
	if (not csig.inherits("CountSignals"))
		Rcpp::stop("expecting a CountSignals object");
	
	//accessing fields
	Rcpp::IntegerVector counts = as<IntegerVector>(gr.slot("counts");
	Rcpp::IntegerVector breaks = as<IntegerVector>(gr.slot("breaks");
	bool ss = as<bool>(gr.slot("ss");
	
	//checking preconditions
	if (ss) Rcpp::stop("expecting a strand-unspecific CountSignals object");
	if (idx < 0 || idx >= breaks.length()-1) Rcpp::stop("invalid index");
	
	//creating matrix
	int len = breaks[idx+1] - breaks[idx];
	int ncol = len/2;
	Rcpp::IntegerMatrix mat(2, ncol);
	memcpy(mat.begin(), counts.begin() + breaks[idx], sizeof(int)*len)
	
	return mat;
}

