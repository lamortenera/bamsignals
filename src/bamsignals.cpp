#include <Rcpp.h>
#include <stdio.h>  
#include <algorithm>
#include "samtools/sam.h"  

//in C: true==1(or something different than 0), false==0
//coding conventions:
//iterator for type t and variable v: It, variable: i_v,
//end of the iterator, e_v.
//use classes, not structs. Name of a class starts with upper-case letter
using namespace Rcpp;

inline bool isNegStrand(const bam1_t *b){
	return ((b->core).flag & BAM_FREVERSE) != 0;
}

inline int fiveprimepos(const bam1_t *b, int strand){
	if (strand){//read on the negative strand
		return bam_calend(&(b->core), bam1_cigar(b))  - 1;
	} else{//read on the reference strand
		return (b->core).pos;
	}
}

//helmuth 2014-03-31
//helmuth 2014-11-10: Consider only first read in a proper mapped pair (SAM FLAG 66)
inline bool isFirstInProperMappedPair(const bam1_t *b){ 
	return ( ((b->core).flag & BAM_FPROPER_PAIR) && ((b->core).flag & BAM_FREAD1) );
}

typedef NumericVector::iterator Idouble;
typedef IntegerVector::iterator Iint;

//workaround, because the function 
//"bam_init_header_hash" in samtools/bam_aux.c
//is not in any header... this should be equivalent
void bam_init_header_hash(bam_header_t *header){
	int foo = 0;
	bam_parse_region(header, "", &foo, &foo, &foo);
}

inline int getRefId(samfile_t* in, const std::string& refname){
	return bam_get_tid(in->header, refname.c_str());
}

//genomic array: represents a correspondence between a genomic range and
//a memory area.
class GArray {
	//initialized when parsing the GRanges
	public:
	
	int rid;//reference id as in the bamfile
	int loc;//0-indexed
	int len;
	int strand;//-1,0,1 represent alternative, unspecified, reference strand respectively
	inline int end(){ return loc + len; }
	//initialized when allocating memory
	int* array;
	int alen;
	GArray(int _rid, int _loc, int _len, int _strand){
		rid = _rid;
		loc = _loc;
		len = _len;
		strand = _strand;
	}
};

//this object won't be valid anymore when the RObject used for construction
//will be destroyed
class RleIter {
	IntegerVector rlens;
	IntegerVector values;
	CharacterVector names;
	
	public:
	
		int run;
		int rlen;
		int rpos;
		bool valid;

		RleIter(RObject& rle):
			rlens(as<IntegerVector>(rle.slot("lengths"))),
			values(as<IntegerVector>(rle.slot("values"))),
			names(as<CharacterVector>(values.attr("levels"))),
			run(0), rpos(-1)
		{
			next();
		}
		
		bool next(){
			++rpos;
			if (rpos == rlens[run]){ //end of the run, go to the next
				++run; rpos = 0;
				if (run == rlens.length())
					valid = false;
					return valid;
			}
			valid = true;
			return valid;
		}
		
		String getValue(){
			return names[values[run]-1];
		}
};

//parses the GR object.
void parseRegions(std::vector<GArray>& container, RObject& gr, samfile_t* in){
	if (not gr.inherits("GRanges"))
		stop("must provide a GRanges object");
	
	IntegerVector starts = as<IntegerVector>(as<RObject>(gr.slot("ranges")).slot("start"));
	IntegerVector lens =   as<IntegerVector>(as<RObject>(gr.slot("ranges")).slot("width"));
	
	RObject chrsRle = as<RObject>(gr.slot("seqnames"));
	RObject strandsRle = as<RObject>(gr.slot("strand"));
	RleIter chrs(chrsRle);
	RleIter strands(strandsRle);
	container.reserve(container.size() + starts.length());
	Iint e_starts = starts.end(); Iint i_starts = starts.begin(); Iint i_lens = lens.begin();
	
	int lastStrandRun = -1;
	int strand = -1;
	
	int lastChrsRun = -1;
	int rid = -1;
	 
	for (; i_starts < e_starts; ++i_starts, ++i_lens, chrs.next(), strands.next()){
		//if new run, update chromosome
		if (lastChrsRun != chrs.run){
			lastChrsRun = chrs.run;
			rid = getRefId(in, chrs.getValue());
			if (rid == -1)
				stop("chromosome " + (std::string)chrs.getValue() + " not present in the bam file");
		}
		
		//if new run, update strand 
		if (lastStrandRun != strands.run){
			lastStrandRun = strands.run;
			const std::string& s = strands.getValue();
			if (s == "-"){ strand = -1; }
			else if (s == "+"){ strand = +1; }
			else { strand = 0; }
		}
		
		container.push_back(GArray(rid, *i_starts - 1, *i_lens, strand));
	}
}

static List allocateAndDistributeMemory(std::vector<GArray>& ranges, int binsize, bool ss){
	int rnum = ranges.size();//number of ranges
	int mult = ss?2:1;
	
	//range i will be stored in the interval breaks[i]-breaks[i+1]
	IntegerVector breaks(rnum+1);
	
	//compute breaks and total length
	double dbinsize = binsize;
	int acc = 0; 
	breaks[0] = acc; 
	for (int i = 0; i < rnum; ++i){
		acc += mult*ceil(ranges[i].len/dbinsize);
		if (acc < 0) Rcpp::stop("Integer overflow: genomic ranges too large");
		breaks[i+1] = acc;
	}
	
	IntegerVector counts(acc);
	int* countsptr = counts.begin();
	
	for (int i = 0; i < rnum; ++i){
		ranges[i].array = countsptr + breaks[i]; 
		ranges[i].alen = breaks[i+1] - breaks[i];
	}
	
	return List::create(_("counts")=counts, _("breaks")=breaks, _("ss")=ss);
}

//if you forget to close the Bamfile you get a memory leak
class Bamfile {
	public:
		samfile_t* in;
		bam_index_t* idx;
		//allocate memory
		Bamfile(const std::string& bampath, int cache_size=10*BGZF_MAX_BLOCK_SIZE){
			const char* cbampath = bampath.c_str();
			in = samopen(cbampath, "rb", 0);  
			if (in == 0) {  
				stop("Fail to open BAM file " + bampath);  
			}
			
			idx = bam_index_load(cbampath); // load BAM index  
			if (idx == 0) {  
				stop("BAM indexing file is not available for file " + bampath);
			}  
			bam_init_header_hash(in->header);
			if (cache_size > 0){
				bgzf_set_cache_size(in->x.bam, cache_size);
			}
		}
		//deallocate
		void close(){
			bam_index_destroy(idx);  
			samclose(in);
		}
};

template <class TRegion>
static bool sortByStart(const TRegion& a, const TRegion& b){
	int ret = b.rid - a.rid;
	if (ret==0){ return b.loc > a.loc; }
	return ret > 0;
}

//interface for TRegion:
//int rid;
//int loc; start position
//int end(); end position
//interface for TPileup:
//void pileup(TRegion&, const bam1_t*, int, int)
template <class TRegion, class TPileup>
static void overlapAndPileup(Bamfile& bfile, std::vector<TRegion>& ranges, int mapqual, int shift, TPileup& pileupper, int maxgap){

	//sorting intervals according to start coordinate (and ref id of course)
	std::sort(ranges.begin(), ranges.end(), sortByStart<TRegion>);
	
	//trade-off between querying a new region and processing unnecessary reads
	const int MAX_GAP = maxgap;

	//variables processed, chunk_start, chunk_end, curr_range and range are indices for the vector ranges
	unsigned int processed = 0;
	bam1_t* read; read = bam_init1();
	//process one chunk of nearby ranges at a time
	while (processed < ranges.size()){
		unsigned int chunk_start = processed;
		int rid = ranges[chunk_start].rid;
		int start = ranges[chunk_start].loc - shift;
		//find out how many regions to process together
		unsigned int chunk_end = chunk_start+1;
		for (; chunk_end < ranges.size(); ++chunk_end){
			if (ranges[chunk_end].rid != rid || ranges[chunk_end].loc - ranges[chunk_end-1].end() - 2*shift > MAX_GAP){
				break;
			}
		}
		int end = ranges[chunk_end-1].end() + shift;
		//perform query
		bam_iter_t iter = bam_iter_query(bfile.idx, rid, start, end);
		//all ranges behind curr_range should not overlap with the next reads anymore
		unsigned int curr_range = chunk_start;
		//loop through the reads
		while (bam_iter_read((bfile.in)->x.bam, iter, read) >= 0){
			if ((read->core).qual >= mapqual){
				int r_start = (read->core).pos;
				//skip non-overlapping regions at the beginning
				while (curr_range < chunk_end && r_start >= ranges[curr_range].end() + shift) ++curr_range;
				//should never happen, unless the last reads returned by the iterator do not overlap with the queried interval
				if (curr_range == chunk_end) break; 
				
				int r_end = bam_calend(&(read->core), bam1_cigar(read)) -1;
				//go through the regions that overlap this read
				for (unsigned int range = curr_range; range < chunk_end && ranges[range].loc - shift <= r_end; ++range){
					pileupper.pileup(ranges[range], read, r_start, r_end);
				}
			}
			
		}
		
		bam_iter_destroy(iter);
		processed = chunk_end;
	}
	bam_destroy1(read);
}

//helmuth 2014-04-08: PairedEnd implementation of overlapAndPileupPairedEnd. This is used by pileup_core() and coverage_core()
//- only considers first read in a proper mapped pair (MAPQ 66, can be positive or negative strand)
template <class TRegion, class TPileup>
static void overlapAndPileupPairedEnd(Bamfile& bfile, std::vector<TRegion>& ranges, int mapqual, int shift, TPileup& pileupper, int maxgap, bool pe_mid, int maxfraglength){

	//sorting intervals according to start coordinate (and ref id of course)
	std::sort(ranges.begin(), ranges.end(), sortByStart<TRegion>);
	
	//trade-off between querying a new region and processing unnecessary reads
	const int MAX_GAP = maxgap;

	int window = shift;
	//if we do paired end midpoint counting we enlarge the range by maxfraglength to get all fragments desired
	if (pe_mid)
	 window = shift + maxfraglength;

	//variables processed, chunk_start, chunk_end, curr_range and range are indices for the vector ranges
	unsigned int processed = 0;
	bam1_t* read; read = bam_init1();
	//process one chunk of nearby ranges at a time
	while (processed < ranges.size()){
		unsigned int chunk_start = processed;
		int rid = ranges[chunk_start].rid;
		int start = ranges[chunk_start].loc - window;
		//find out how many regions to process together
		unsigned int chunk_end = chunk_start+1;
		for (; chunk_end < ranges.size(); ++chunk_end){
			if (ranges[chunk_end].rid != rid || ranges[chunk_end].loc - ranges[chunk_end-1].end() - 2*window > MAX_GAP){
				break;
			}
		}
		int end = ranges[chunk_end-1].end() + window;
		//perform query
		bam_iter_t iter = bam_iter_query(bfile.idx, rid, start, end); 
		//all ranges behind curr_range should not overlap with the next reads anymore
		unsigned int curr_range = chunk_start;
		//loop through the reads
		while (bam_iter_read((bfile.in)->x.bam, iter, read) >= 0){
			//only take first read in proper pair mapping + mapq threshold
			if ( isFirstInProperMappedPair( read ) && ( (read->core).qual >= mapqual) ){ 
				int r_start = (read->core).pos;
				//skip non-overlapping regions at the beginning
				while (curr_range < chunk_end && r_start >= ranges[curr_range].end() + window) ++curr_range;
				//should never happen, unless the last reads returned by the iterator do not overlap with the queried interval
				if (curr_range == chunk_end) break; 

				int r_end = bam_calend(&(read->core), bam1_cigar(read)) -1;
				int r_shift = shift;

				//if we want to count midpoints of the fragments we have to change r_end
				if (pe_mid) {
					//move counting position relative to fragment middle point
					r_shift = r_shift + abs((read->core).isize)/2; 
				}
				//go through the regions that overlap this read
				for (unsigned int range = curr_range; range < chunk_end && ranges[range].loc - window <= r_end; ++range){
					pileupper.pileupPairedEnd(ranges[range], read, r_start, r_end, r_shift);
				}
			}
		}
		
		bam_iter_destroy(iter);
		processed = chunk_end;
	}
	bam_destroy1(read);
}

class Coverager{
	public:
	
	void pileup(GArray& range, const bam1_t* read, int start, int end){
		//check if the read really overlaps
		if (start < range.end() && end >= range.loc){
		    	//Rcout << "-> took read " << std::endl;
			if (range.strand >= 0){
				//range on the reference strand
				int pos = start-range.loc;
				++range.array[pos>0?pos:0];
				pos = end + 1 - range.loc;
				if (pos < range.len){
					--range.array[pos];
				}
			} else {
				//range on the negative strand
				int pos = range.end() - 1 - end;
				++range.array[pos>0?pos:0];
				pos = range.end() - start;
				if (pos < range.len){
					--range.array[pos];
				}
			}
		}
	}

	//helmuth 2014-04-08: "shift" argument added. Paired End extension for varying shift sizes. 
	//                     It's not used in this method but incorporated for consistency with 
	//                     TPileup::Pileupper.
	void pileupPairedEnd(GArray& range, const bam1_t* read, int start, int end, int shift){
		//Construct the region of the fragment overlap
		int isize = (read->core).isize;
		bool negstrand = isNegStrand(read);
		if (negstrand && isize < 0) {        //-strand read: only calculate if isize is meaningful, otherwise fall back to given start
		    start = end + isize;
		    //Rcout << "-FRAG: " << start << "-" << end << " ( pos = " << (read->core).pos << ", read_end = " << bam_calend(&(read->core), bam1_cigar(read)) << " )" <<  std::endl;
		} else if (!negstrand && isize > 0) { //+strand read: only calculate if isize is meaningful, otherwise fall back to given end (i.e. bam.c::bam_calend output)
		    end   = start + isize - 1;
		    //Rcout << "+FRAG: " << start << "-" << end << " ( pos = " << (read->core).pos << ", read_end = " << bam_calend(&(read->core), bam1_cigar(read)) << " )" <<  std::endl;
		}
		pileup( range, read, start, end);
	}
};

class Pileupper{
	public:
	
	int binsize;
	int shift; 
	bool ss;
	
	Pileupper(int abinsize, int ashift, bool ass){
		binsize = abinsize;
		shift = ashift;
		ss = ass;
	}
	
	//can be accelerated by rewriting a similar one for special cases:
	//ss=false, ignoreRegionStrand, binsize=1
	//and storing in range extra information: range.loc - shift and range.loc + shift
	//this could be probably done by templating this class (so everything would be written just once).
	//
	void pileup(GArray& range, const bam1_t* read, int start, int end){
		bool negstrand = isNegStrand(read);
		//relative position from the start of the range
		int pos = negstrand ? end - range.loc - shift : start + shift - range.loc;
		//check overlap with the region
		if (pos < 0 || pos >= range.len) return;
		//strand of the region defines direction and sense and antisense
		if (range.strand < 0){
			pos = range.len - pos - 1;
			negstrand = !negstrand;
		}
		//even positions of the array are sense-reads, odd are antisense
		if (ss){ ++range.array[2*(pos/binsize) + (negstrand?1:0)]; }
		else {++range.array[pos/binsize]; } 
	}

	//helmuth 2014-04-08: "shift" argument added. Paired End extension for varying shift sizes. 
	//                     It's not used in this method but incorporated for consistency with 
	//                     TPileup.
	void pileupPairedEnd(GArray& range, const bam1_t* read, int start, int end, int r_shift){
		bool negstrand = isNegStrand(read);
		//relative position from the start of the range
		int pos = negstrand ? end - range.loc - r_shift : start + r_shift - range.loc;
		//check overlap with the region
		if (pos < 0 || pos >= range.len) return;
		//strand of the region defines direction and sense and antisense
		if (range.strand < 0){
			pos = range.len - pos - 1;
			negstrand = !negstrand;
		}
		//even positions of the array are sense-reads, odd are antisense
		if (ss){ ++range.array[2*(pos/binsize) + (negstrand?1:0)]; }
		else {++range.array[pos/binsize]; } 
	}
};

// [[Rcpp::export]]
List pileup_core(RObject gr, std::string bampath, int mapqual=0, int binsize=1, int shift=0, bool ss=false, 
	bool pe=false, bool pe_mid=false, int maxfraglength=1000, int maxgap=16385){
	std::vector<GArray> ranges;
	//opening bamfile and index
	Bamfile bfile(bampath);
	//adding regions to the vector
	parseRegions(ranges, gr, bfile.in);
	//allocate memory
	List ret = allocateAndDistributeMemory(ranges, binsize, ss);
	//pileup
	Pileupper p(binsize, shift, ss);
	if (pe) //use PairedEnd routine
	 overlapAndPileupPairedEnd(bfile, ranges, mapqual, shift, p, maxgap, pe_mid, maxfraglength);
	else //use SingleEnd routine
	 overlapAndPileup(bfile, ranges, mapqual, shift, p, maxgap);
	//close bamfile and index
	bfile.close();
	
	return ret;
}

// [[Rcpp::export]]
List coverage_core(RObject gr, std::string bampath, int mapqual=0, bool pe=false, int maxfraglength=1000, int maxgap=16385){
	std::vector<GArray> ranges;
	//opening bamfile and index
	Bamfile bfile(bampath);
	//adding regions to the vector
	parseRegions(ranges, gr, bfile.in);
	//allocate memory
	List ret = allocateAndDistributeMemory(ranges, 1, false);
	//pileup
	Coverager c;
	if (pe)
	 overlapAndPileupPairedEnd(bfile, ranges, mapqual, 0, c, maxgap, true, maxfraglength); //we set pe.mid=true to grasp all regions in fragment size
	else
	 overlapAndPileup(bfile, ranges, mapqual, 0, c, maxgap);
	//close bamfile and index
	bfile.close();
	//do cumsum on all the ranges
	typedef std::vector<GArray>::iterator IGArray;
	IGArray i_ranges = ranges.begin(); IGArray e_ranges = ranges.end(); 
	for (; i_ranges < e_ranges; ++i_ranges){
		int acc = i_ranges->array[0];
		int* e_array = i_ranges->array + i_ranges->len;
		for (int* i_array = i_ranges->array + 1; i_array < e_array; ++i_array){
			*i_array = acc += *i_array;
		}
	}
	
	return ret;
}


// [[Rcpp::export]]
bool writeSamAsBamAndIndex(const std::string& sampath, const std::string& bampath) {
	//streams
	samfile_t *in = 0, *out = 0;

	//open samfile for reading
	const char* csampath = sampath.c_str();
	in = samopen(csampath, "r", 0);
	if (in == 0) {  
		stop("Fail to open SAM file " + sampath);  
	}

	//open bamfile for writing
	const char* cbampath = bampath.c_str();
	out = samopen(cbampath, "wb", in->header);
	if (out == 0) {  
		stop("Fail to open BAM file ." + bampath);  
	}

	//read sam and write to bam: adapted from sam_view.c:232-244
	bam1_t *b = bam_init1();
	while (samread(in, b) >= 0) { // read one alignment from `in'
		samwrite(out, b); // write the alignment to `out'
	}
	bam_destroy1(b);

	//close streams
	samclose(in);
	samclose(out);

	//build the index
	bam_index_build(bampath.c_str());

	return 0;
}
