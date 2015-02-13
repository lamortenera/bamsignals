#include <Rcpp.h>
//that's the equivalent of the s4 class with the same name
struct CountSignals {
	Rcpp::IntegerVector counts;
	Rcpp::IntegerVector breaks;
	bool ss;
	
	CountSignals(SEXP sexp) {
		Rcpp::RObject csig = Rcpp::as<Rcpp::RObject>(sexp);
		if (not csig.inherits("CountSignals")) Rcpp::stop("expecting a CountSignals object");
		counts = Rcpp::as<Rcpp::IntegerVector>(csig.slot("counts"));
		breaks = Rcpp::as<Rcpp::IntegerVector>(csig.slot("breaks"));
		//if I don't put Rcpp::as<bool> some weid bugs show up
		ss = Rcpp::as<bool>(csig.slot("ss"));
	}
};
