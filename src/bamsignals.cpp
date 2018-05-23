#include <Rcpp.h>
#include <stdio.h>  
#include <algorithm>
#include "htslib/sam.h"
#include "htslib/bgzf.h"

//in C: true==1(or something different than 0), false==0
using namespace Rcpp;

//returns true if the read is in the negative strand
inline bool isNegStrand(const bam1_t *b){
    return ((b->core).flag & BAM_FREVERSE) != 0;
}

//leftmost position of a read (0-based, inclusive)
inline int readEnd(const bam1_t *b){
    return bam_endpos(b) - 1;
}

//it returns true if any of the bits in the mask is not set
inline bool invalidFlag(const bam1_t *b, uint32_t mask){
    return (mask & ~(b->core).flag);
}

//gets the reference id of a chromosome as stored in the bam file
inline int getRefId(bam_hdr_t* header, const std::string& refname){
    return bam_name2id(header, refname.c_str());
}

//genomic array: represents a correspondence between a genomic range and
//a memory area.
class GArray {
    //initialized when parsing the GRanges
    public:
    
    int rid;//reference id as in the bamfile
    int loc;//0-indexed
    int len;//num of bases in the genomic range
    int strand;//-1,0,1 represent alternative, unspecified, reference strand respectively
    inline int end(){ return loc + len; }
    //initialized when allocating memory
    int* array;//array where counts are stored
    int alen;//num of elements of the array
    GArray(int _rid, int _loc, int _len, int _strand){
        rid = _rid;
        loc = _loc;
        len = _len;
        strand = _strand;
    }
};

//iterate through a run-length-encoded character vector (s4 class Rle)
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
        bool valid;//iff(run < rlens.length() && rpos < rlens[run])

        RleIter(RObject& rle):
            rlens(as<IntegerVector>(rle.slot("lengths"))),
            values(as<IntegerVector>(rle.slot("values"))),
            names(as<CharacterVector>(values.attr("levels"))),
            run(0), rpos(-1), valid(true)
        {
            next();
        }
        
        bool next(){
            if (!valid) return false;
            ++rpos;
            if (rpos == rlens[run]){ //end of the run, go to the next
                ++run; rpos = 0;
                valid = run < rlens.length();
            }
            return valid;
        }
        
        String getValue(){
            return names[values[run]-1];
        }
};

//parses the GR object.
void parseRegions(std::vector<GArray>& container, RObject& gr, samFile* in){
    if (not gr.inherits("GRanges"))
        stop("must provide a GRanges object");
    bam_hdr_t *header = sam_hdr_read(in);
    
    IntegerVector starts = as<IntegerVector>(as<RObject>(gr.slot("ranges")).slot("start"));
    IntegerVector lens =   as<IntegerVector>(as<RObject>(gr.slot("ranges")).slot("width"));
    int nranges = starts.length();
    
    RObject chrsRle = as<RObject>(gr.slot("seqnames"));
    RObject strandsRle = as<RObject>(gr.slot("strand"));
    RleIter chrs(chrsRle);
    RleIter strands(strandsRle);
    container.reserve(container.size() + nranges);
    
    int lastStrandRun = -1;
    int strand = -1;
    
    int lastChrsRun = -1;
    int rid = -1;
     
    for (int i = 0; i < nranges; ++i, chrs.next(), strands.next()){
        //if new run, update chromosome
        if (lastChrsRun != chrs.run){
            lastChrsRun = chrs.run;
            rid = getRefId(header, chrs.getValue());
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
        
        container.push_back(GArray(rid, starts[i] - 1, lens[i], strand));
    }
    
    bam_hdr_destroy(header);
}

//allocates the memory and sets the pointer for each GArray object
//it also sets the binsize if binsize <= 0
static List allocateList(std::vector<GArray>& ranges, int* binsize, bool ss){
    int rnum = ranges.size();//number of ranges
    int mult = ss?2:1;
    double dbinsize = *binsize;
    //set the dimnames (this is needed only if ss==TRUE)
    List dnames(2);
    if (ss) dnames[0] = CharacterVector::create("sense", "antisense");
    
    //allocate a list with one vector (or matrix) containing all counts
    if (dbinsize <= 0){
        List sigs(1);
        int* C;
        if (ss){
          IntegerMatrix sig(2, rnum);
          sig.attr("dimnames") = dnames;
          sigs[0] = sig;
          C = sig.begin();
        } else {
          IntegerVector sig(rnum);
          sigs[0] = sig;
          C = sig.begin();
        }
        //assign pointers and compute binsize for the pileupper
        int maxw = -1;
        for (int i = 0; i < rnum; ++i){
            ranges[i].array = C + mult*i;
            maxw = std::max(maxw, ranges[i].len);
        }
        *binsize = maxw;
        return sigs;
    }
    
    //allocate a list where range i will be stored as the i-th element of sigs
    List sigs(rnum);
    for (int i = 0; i < rnum; ++i){
        //the width does not depend on the 'ss' parameter
        int width = ceil(ranges[i].len/dbinsize);
        if (ss){
            //if strand specific, allocate matrix
            IntegerMatrix sig(2, width);
            sig.attr("dimnames") = dnames;
            sigs[i] = sig;
            ranges[i].array = sig.begin(); 
        } else {
            //otherwise, allocate vector
            IntegerVector sig(width);
            sigs[i] = sig;
            ranges[i].array = sig.begin(); 
        }
        //set the pointer of the range
        ranges[i].alen = mult*width;
    }
    return sigs;
}

//if you forget to close the Bamfile you get a memory leak
class Bamfile {
    public:
        samFile* in;
        hts_idx_t* idx;
        //allocate memory
        Bamfile(const std::string& bampath, int cache_size=10*BGZF_MAX_BLOCK_SIZE){
            const char* cbampath = bampath.c_str();
            in = sam_open(cbampath, "rb");
            if (in == NULL) {  
                stop("Fail to open BAM file " + bampath);  
            }
            
            idx = bam_index_load(cbampath); // load BAM index  
            if (idx == 0) {  
                stop("BAM indexing file is not available for file " + bampath);
            }  
            if (cache_size > 0){
                bgzf_set_cache_size(in->fp.bgzf, cache_size);
            }
        }
        //deallocate
        ~Bamfile(){
            hts_idx_destroy(idx);  
            sam_close(in);
        }
};

static inline bool sortByStart(const GArray& a, const GArray& b){
    int ret = b.rid - a.rid;
    if (ret==0){ return b.loc > a.loc; }
    return ret > 0;
}

//this function overlaps each read with the regions it might fall into
//each read is processed only once, and when it is matched to a region
//the function TPileup.pileup is called

//interface for TPileup:
//1. int setRead(const bam1_t*): check compatibility of the read, if compatible
//set the internal variables of the pileupper and return the end of the read 
//(0-based, inclusive, we don't want to compute it twice), otherwise -1
//2. void pileup(GArray&): pileup the last set read with the given interval

//the parameter 'ext' is an extension added to beginning and end of a read
//when computing the ranges that it overlaps
template <class TPileup>
static void overlapAndPileup(Bamfile& bfile, std::vector<GArray>& ranges, 
                            int ext, TPileup& pileupper, int maxgap){
    if (ext < 0) Rcpp::stop("negative 'ext' values don't make sense");
    
    //sorting intervals according to start coordinate (and ref id of course)
    std::sort(ranges.begin(), ranges.end(), sortByStart);
    
    //variables processed, chunk_start, chunk_end, curr_range and range are indices for the vector ranges
    unsigned int processed = 0;
    bam1_t* read; read = bam_init1();
    //process one chunk of nearby ranges at a time
    while (processed < ranges.size()){
        unsigned int chunk_start = processed;
        int rid = ranges[chunk_start].rid;
        int start = ranges[chunk_start].loc - ext;
        int end = ranges[chunk_start].end() + ext;
        //find out how many regions to process together
        unsigned int chunk_end = chunk_start+1;
        for (; chunk_end < ranges.size(); ++chunk_end){
            int next_start = ranges[chunk_end].loc - ext;
            if (ranges[chunk_end].rid != rid || next_start - end > maxgap){
                break;
            }
            end = std::max(end, ranges[chunk_end].end() + ext);
        }
        //perform query
        hts_itr_t* iter = bam_itr_queryi(bfile.idx, rid, start, end);
        //all ranges behind curr_range should not overlap with the next reads anymore
        unsigned int curr_range = chunk_start;
        //loop through the reads
        while (bam_itr_next(bfile.in, iter, read) >= 0){
            int read_end = pileupper.setRead(read);
            if (read_end < 0) continue;
            //overlap start end end (extremes included)
            int ov_start = (read->core).pos - ext;
            int ov_end = read_end + ext;
            //skip non-overlapping regions at the beginning
            while (curr_range < chunk_end && ov_start >= ranges[curr_range].end()) ++curr_range;
            //should never happen
            if (curr_range == chunk_end) break; 
            //go through the regions that overlap this read (with extension)
            for (unsigned range = curr_range; 
                range < chunk_end && ranges[range].loc <= ov_end; ++range){
                pileupper.pileup(ranges[range]);
            }
        }
        bam_itr_destroy(iter);
        processed = chunk_end;
    }
    bam_destroy1(read);
}

class Pileupper{
    public:
    
    const int binsize;
    const int shift; 
    const bool ss;
    const int mapqual;
    const int requiredF;//a required flag
    const bool midpoint;//consider the midpoint or not
    const int filteredF;//a filtering flag
    const int* tlen_filter; //filter on sam TLEN
    
    //these values refer to the last read and are set using "setRead"
    
    //0-based position of the 5' end base pair from the start of the chromosome
    int pos;
    //true if is on the reference strand, false otherwise
    bool negstrand;
    
    Pileupper(int abinsize, int ashift, bool ass, int amapqual, 
              int arequiredF, bool amidpoint, int afilteredF, 
              int* atlen_filter) : 
        binsize(abinsize), shift(ashift), ss(ass), mapqual(amapqual), 
        requiredF(arequiredF), midpoint(amidpoint), filteredF(afilteredF),
        tlen_filter(atlen_filter)
    {
        pos = -1;
        negstrand = false;
    }
    
    //it filter the reads if they are not ok
    //and sets the relevant variables if they are ok
    //it returns the end of the read if everything is ok, -1 otherwise
    inline int setRead(const bam1_t* read){
        //filter the read on MAPQ, SAMFLAG and ISIZE
        if ((read->core).qual < mapqual || invalidFlag(read, requiredF)
            || !invalidFlag(read, filteredF)
            || (tlen_filter != 0 && 
                (abs((read->core).isize) < tlen_filter[0] || 
                 abs((read->core).isize) > tlen_filter[1])))
             return -1;
        //compute read end
        int read_end = readEnd(read);
        //set 'negstrand' and 'pos'
        negstrand = isNegStrand(read);
        //shift in the 5' direction
        int offset = midpoint?(abs((read->core).isize)/2 + shift):shift; 
        if (negstrand) {
            pos = read_end - offset;
        } else {
            pos = (read->core).pos + offset;
        }
        return read_end;
    }
    
    //it accumulates the last set read with a given range
    inline void pileup(GArray& range){
        //relative position from the start of the range
        int relpos = pos - range.loc;
        //check overlap with the region
        if (relpos < 0 || relpos >= range.len) return;
        //strand of the region defines direction and sense and antisense
        int antisense = negstrand?1:0;
        if (range.strand < 0){
            relpos = range.len - relpos - 1;
            antisense = 1-antisense;
        }
        //even positions of the array are sense-reads, odd are antisense
        if (ss){ ++range.array[2*(relpos/binsize) + antisense]; }
        else {++range.array[relpos/binsize]; } 
    }
};

class Coverager{
    public:
    
    const int mapqual;
    const int requiredF;//a required flag
    const int filteredF;//a filtering flag
    const bool tspan;//consider the span of the whole read pair or not
    const int* tlen_filter; //filter on sam TLEN
    
    //these values refer to the last read and are set using "setRead"
    //'start' and 'end' (0-based, inclusive) 
    int start;
    int end;
    
    Coverager(int amapqual, int arequiredF, int afilteredF, bool atspan, 
              int* atlen_filter) :
        mapqual(amapqual), requiredF(arequiredF), filteredF(afilteredF), 
        tspan(atspan), tlen_filter(atlen_filter)
    {
        start = end = -1;
    }
    
    
    //it filter the reads if they are not ok
    //and sets the relevant variables if they are ok
    //it returns the end of the read if everything is ok, -1 otherwise
    inline int setRead(const bam1_t* read){
        //filter the read on MAPQ, SAMFLAG and ISIZE
        if ((read->core).qual < mapqual || invalidFlag(read, requiredF) 
            || !invalidFlag(read, filteredF)
            || (tlen_filter != 0 && 
               (abs((read->core).isize) < tlen_filter[0] || 
                abs((read->core).isize) > tlen_filter[1])))
             return -1;
        //set 'start' and 'end' (0-based, inclusive)
        int read_end = readEnd(read);
        start = (read->core).pos;
        end = read_end;
        if (tspan){
            bool negstrand = isNegStrand(read);
            int isize = (read->core).isize;
            //only calculate if isize is meaningful, otherwise fall back to given start
            if (negstrand && isize < 0) {
                start = end + isize + 1;
            } else if (!negstrand && isize > 0) {
                end   = start + isize - 1;
            }
        }
        return read_end;
    }
    
    //it accumulates the last set read with a given range
    inline void pileup(GArray& range){
        //check if the read really overlaps
        if (start >= range.end() || end < range.loc) return;
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
};

//this function is called by bamProfile and bamCount.
//when it is bamCount, then binsize <= 0
// [[Rcpp::export]]
List pileup_core(std::string bampath, RObject& gr, IntegerVector& tlen_filter, 
    int mapqual=0, int binsize=1, int shift=0, bool ss=false, int requiredF=0, 
    int filteredF=-1, bool pe_mid=false, int maxgap=16385){
    std::vector<GArray> ranges;
    //opening bamfile and index
    Bamfile bfile(bampath);
    //adding regions to the vector
    parseRegions(ranges, gr, bfile.in);
    //allocate memory (if binsize <= 0, binsize is set to max(width(gr)))
    List ret = allocateList(ranges, &binsize, ss);
    //pileup
    Pileupper p(binsize, shift, ss, mapqual, requiredF, pe_mid, filteredF,
                (tlen_filter.size()==0?0:&tlen_filter[0]));
    int ext = abs(shift) + (pe_mid?tlen_filter[1]:0);
    overlapAndPileup(bfile, ranges, ext, p, maxgap);
    
    return ret; 
}

//do cumulative sum
inline void cumsum(int* C, int len){
    if (len < 2) return;
    int acc = C[0];
    for (int i = 1; i < len; ++i){
      C[i] = (acc += C[i]);
    }
}

//this function is called by bamCoverage
// [[Rcpp::export]]
List coverage_core(std::string bampath, RObject& gr, IntegerVector& tlen_filter, 
    int mapqual=0, int requiredF=0, int filteredF=-1, bool tspan=false, 
    int maxgap=16385){
    std::vector<GArray> ranges;
    //opening bamfile and index
    Bamfile bfile(bampath);
    //adding regions to the vector
    parseRegions(ranges, gr, bfile.in);
    //allocate memory
    int binsize = 1;
    List ret = allocateList(ranges, &binsize, false);
    Coverager c(mapqual, requiredF, filteredF, tspan, 
                (tlen_filter.size()==0?0:&tlen_filter[0]));
    int ext = tspan?tlen_filter[1]:0;
    overlapAndPileup(bfile, ranges, ext, c, maxgap);
    //do cumsum on all the ranges
    for (int i = 0, e = ranges.size(); i < e; ++i){
        cumsum(ranges[i].array, ranges[i].len);
    }
    return ret;
}

//this is used only for the tests right now...
// [[Rcpp::export]]
bool writeSamAsBamAndIndex(const std::string& sampath, const std::string& bampath) {
    //streams
    samFile *in = 0, *out = 0;

    //open samfile for reading
    const char* csampath = sampath.c_str();
    in = sam_open(csampath, "r");
    if (in == 0) {  
        stop("Fail to open SAM file " + sampath);  
    }
    bam_hdr_t *header = sam_hdr_read(in);

    //open bamfile for writing
    const char* cbampath = bampath.c_str();
    out = sam_open(cbampath, "wb");
    if (out == 0) {  
        stop("Fail to open BAM file ." + bampath);  
    }
    sam_hdr_write(out, header);

    //read sam and write to bam: adapted from sam_view.c:232-244
    bam1_t *b = bam_init1();
    while (sam_read1(in, header, b) >= 0) { // read one alignment from `in'
        bam_write1(out->fp.bgzf, b); // write the alignment to `out'
    }
    bam_destroy1(b);

    //close streams
    sam_close(in);
    sam_close(out);

    //build the index
    bam_index_build(bampath.c_str(), 0);

    return true;
}
