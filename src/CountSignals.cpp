#include <Rcpp.h>

//allocates the counts, breaks and ss slots
#define unwrapCountSignals(csig) \
	if (not csig.inherits("CountSignals")) Rcpp::stop("expecting a CountSignals object");\
	Rcpp::IntegerVector counts = csig.slot("counts");\
	Rcpp::IntegerVector breaks = csig.slot("breaks");\
	bool ss = csig.slot("ss");\


inline void checkIndex(int idx, Rcpp::IntegerVector& breaks){
	if (idx < 0 || idx >= breaks.length()-1) Rcpp::stop("index out of range");
}

// [[Rcpp::export]]
Rcpp::IntegerMatrix getMatrix(Rcpp::RObject csig, int idx){
	unwrapCountSignals(csig);
	
	//R idx to C idx
	--idx;
	
	//checking preconditions
	if (!ss) Rcpp::stop("expecting a strand-specific CountSignals object");
	checkIndex(idx, breaks);
	
	//creating matrix
	int len = breaks[idx+1] - breaks[idx];
	int ncol = len/2;
	Rcpp::IntegerMatrix mat(2, ncol);
	memcpy(mat.begin(), counts.begin() + breaks[idx], sizeof(int)*len);
	
	//set the dimnames
	Rcpp::List dnames(2);
	dnames[0] = Rcpp::CharacterVector::create("sense", "antisense");
	mat.attr("dimnames") = dnames;
	
	return mat;
}

// [[Rcpp::export]]
Rcpp::IntegerVector getVector(Rcpp::RObject csig, int idx){
	unwrapCountSignals(csig);
	
	//R idx to C idx
	--idx;
	
	//checking preconditions
	if (ss) Rcpp::stop("expecting a strand-unspecific CountSignals object");
	checkIndex(idx, breaks);
	
	//creating vector
	int len = breaks[idx+1] - breaks[idx];
	Rcpp::IntegerVector vec(len);
	memcpy(vec.begin(), counts.begin() + breaks[idx], sizeof(int)*len);
	
	return vec;
}


// [[Rcpp::export]]
Rcpp::List getSubset(Rcpp::RObject csig, Rcpp::IntegerVector idxs){
	unwrapCountSignals(csig);
	
	int nlen = idxs.length();
	Rcpp::IntegerVector nbreaks(nlen+1);
	//figure out starts and ends in the new vector
	//and total amount of memory to allocate
	int acc = 0; nbreaks[0]=0;
	for (int i = 0; i < nlen; ++i){
		//R idx to C idx
		int idx = idxs[i] - 1;
		checkIndex(idx, breaks);
		acc += breaks[idx+1] - breaks[idx];
		nbreaks[i+1] = acc;
	}
	
	//copy the subsets in a new vector
	Rcpp::IntegerVector ncounts(acc);
	for (int i = 0; i < nlen; ++i){
		//R idx to C idx
		int idx = idxs[i] - 1;
		int len = nbreaks[i+1] - nbreaks[i];
		memcpy(ncounts.begin() + nbreaks[i], counts.begin() + breaks[idx], sizeof(int)*len);
	}
	
	return Rcpp::List::create(Rcpp::Named("counts")=ncounts, Rcpp::Named("breaks")=nbreaks, Rcpp::Named("ss")=ss);
}

// [[Rcpp::export]]
Rcpp::List asList(Rcpp::RObject csig){
	unwrapCountSignals(csig);
	
	int nsig = breaks.length()-1;
	Rcpp::List list(nsig);
	
	for (int i = 0; i < nsig; ++i){
		if (ss) list[i] = getMatrix(csig, i+1);
		else list[i]    = getVector(csig, i+1);
	}
	
	return list;
}

// [[Rcpp::export]]
Rcpp::IntegerVector fastWidth(Rcpp::RObject csig){
	if (not csig.inherits("CountSignals")) Rcpp::stop("expecting a CountSignals object");
	Rcpp::IntegerVector breaks = csig.slot("breaks");
	bool ss = csig.slot("ss");
	int div = ss?2:1;
	int nsig = breaks.length()-1;

	Rcpp::IntegerVector w(nsig);
	int last = breaks[0];
	for (int i = 0; i < nsig; ++i){
		int next = breaks[i+1];
		w[i] = (next-last)/div;
		last = next;
	}
	
	return w;
}
