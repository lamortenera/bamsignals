#include <Rcpp.h>
#include <bamsignals.h>

inline void checkIndex(int idx, Rcpp::IntegerVector& breaks){
    if (idx < 0 || idx >= breaks.length()-1) Rcpp::stop("index out of range");
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix getMatrix(CountSignals x, int idx){
    //R idx to C idx
    --idx;
    
    //checking preconditions
    if (!x.ss) Rcpp::stop("expecting a strand-specific CountSignals object");
    checkIndex(idx, x.breaks);
    
    //creating matrix
    int len = x.breaks[idx+1] - x.breaks[idx];
    int ncol = len/2;
    Rcpp::IntegerMatrix mat(2, ncol);
    memcpy(mat.begin(), x.counts.begin() + x.breaks[idx], sizeof(int)*len);
    
    //set the dimnames
    Rcpp::List dnames(2);
    dnames[0] = Rcpp::CharacterVector::create("sense", "antisense");
    mat.attr("dimnames") = dnames;
    
    return mat;
}

// [[Rcpp::export]]
Rcpp::IntegerVector getVector(CountSignals x, int idx){
    //R idx to C idx
    --idx;
    
    //checking preconditions
    if (x.ss) Rcpp::stop("expecting a strand-unspecific CountSignals object");
    checkIndex(idx, x.breaks);
    
    //creating vector
    int len = x.breaks[idx+1] - x.breaks[idx];
    Rcpp::IntegerVector vec(len);
    memcpy(vec.begin(), x.counts.begin() + x.breaks[idx], sizeof(int)*len);
    
    return vec;
}


// [[Rcpp::export]]
Rcpp::List getSubset(CountSignals x, Rcpp::IntegerVector idxs){
    int nlen = idxs.length();
    Rcpp::IntegerVector nbreaks(nlen+1);
    //figure out starts and ends in the new vector
    //and total amount of memory to allocate
    int acc = 0; nbreaks[0]=0;
    for (int i = 0; i < nlen; ++i){
        //R idx to C idx
        int idx = idxs[i] - 1;
        checkIndex(idx, x.breaks);
        acc += x.breaks[idx+1] - x.breaks[idx];
        nbreaks[i+1] = acc;
    }
    
    //copy the subsets in a new vector
    Rcpp::IntegerVector ncounts(acc);
    for (int i = 0; i < nlen; ++i){
        //R idx to C idx
        int idx = idxs[i] - 1;
        int len = nbreaks[i+1] - nbreaks[i];
        memcpy(ncounts.begin() + nbreaks[i], x.counts.begin() + x.breaks[idx], sizeof(int)*len);
    }
    
    return Rcpp::List::create(Rcpp::Named("counts")=ncounts, Rcpp::Named("breaks")=nbreaks, Rcpp::Named("ss")=x.ss);
}

// [[Rcpp::export]]
Rcpp::List asList(CountSignals x){
    int nsig = x.breaks.length()-1;
    Rcpp::List list(nsig);
    
    for (int i = 0; i < nsig; ++i){
        if (x.ss) list[i] = getMatrix(x, i+1);
        else list[i]    = getVector(x, i+1);
    }
    
    return list;
}

// [[Rcpp::export]]
Rcpp::IntegerVector fastWidth(CountSignals x){
    int div = x.ss?2:1;
    int nsig = x.breaks.length()-1;

    Rcpp::IntegerVector w(nsig);
    int last = x.breaks[0];
    for (int i = 0; i < nsig; ++i){
        int next = x.breaks[i+1];
        w[i] = (next-last)/div;
        last = next;
    }
    
    return w;
}
