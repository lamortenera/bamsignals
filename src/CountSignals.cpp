#include <Rcpp.h>

// [[Rcpp::export]]
bool checkList(Rcpp::List l, bool ss){
    int nsig = l.length();
    for (int i = 0; i < nsig; ++i){
        //check that we have integers
        if (TYPEOF(l[i]) != INTSXP) return false;
        if (ss){
        //check that we have a matrix with 2 rows
        Rcpp::IntegerVector dims = Rf_getAttrib(l[i], R_DimSymbol);
        if (dims.length() != 2 || dims[0] != 2) return false;
        }
    }
    return true;
}

// [[Rcpp::export]]
Rcpp::IntegerVector fastWidth(Rcpp::List l, bool ss){
    int nsig = l.length();
    int div = ss?2:1;
    Rcpp::IntegerVector w(nsig);
    for (int i = 0; i < nsig; ++i){
        Rcpp::IntegerVector sig = l[i];
        w[i] = sig.length()/div;
    }
    
    return w;
}
