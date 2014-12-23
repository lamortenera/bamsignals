# This script is executed by ../compile.R to create reference counts to test bamsignals functionality
#
#*NOTE*
#- BEDFORMAT is 	0-based 	+	end-exclusive
#- bamsignals is	0-based		+	end-exclusive
#- BAMFORMAT is		0-based		+	end-exclusive
#
#- GRanges is		1-based		+	end-inclusive
#
#- SAMFORMAT is		1-based		+	end-exclusive
#
# Create reference counts with bedtools
# Toy_SE.bam 5'-end counting
samtools view Toy_SE.bam | awk '$1==0 {printf "%s\t%s\t%s\t.\t%s\t+\n", $2, $3-1, $3, $4}; $1==16 {printf "%s\t%s\t%s\t.\t%s\t-\n", $2, $3-1+$5-1, $3+$5-1, $4}' > Toy_SE.5p.bed 
samtools view Toy_SE.bam | awk '$1==0 {printf "%s\t%s\t%s\t.\t%s\t+\n", $2, $3-1+100, $3+100, $4}; $1==16 {printf "%s\t%s\t%s\t.\t%s\t-\n", $2, $3-1+$5-1-100, $3+$5-1-100, $4}' > Toy_SE.5p+100.bed 
samtools view -q 40 Toy_SE.bam | awk '$1==0 {printf "%s\t%s\t%s\t.\t%s\t+\n", $2, $3-1, $3, $4}; $1==16 {printf "%s\t%s\t%s\t.\t%s\t-\n", $2, $3-1+$5-1, $3+$5-1, $4}' > Toy_SE.5p+MAPQ40.bed 
bedtools coverage -a Toy_SE.5p.bed -b regions.bed | cut -f 7 > 5p.all
bedtools coverage -s -a Toy_SE.5p.bed -b regions.bed | cut -f 7 > 5p.sense
bedtools coverage -S -a Toy_SE.5p.bed -b regions.bed | cut -f 7 > 5p.antisense
bedtools coverage -a Toy_SE.5p+100.bed -b regions.bed | cut -f 7 > 5p.shift
bedtools coverage -a Toy_SE.5p+MAPQ40.bed -b regions.bed | cut -f 7 > 5p.qual
echo -e "#chrom\tstart\tend\tname\tscore\tstrand\tFivePrime_count\tFivePrime_count_sense\tFivePrime_count_antisense\tFivePrime_count_100_shift\tFivePrime_count_min40qual" > Toy_SE.counts
paste regions.bed 5p.all 5p.sense 5p.antisense 5p.shift 5p.qual >> Toy_SE.counts

# Toy_SE.bam profiles
bedtools coverage -d -a Toy_SE.5p.bed -b regions.bed | cut -f1,2,3,4,5,6,8 > Toy_SE.5p.profile
bedtools coverage -d -s -a Toy_SE.5p.bed -b regions.bed | cut -f 8 > Toy_SE.5p.sense.profile
bedtools coverage -d -S -a Toy_SE.5p.bed -b regions.bed | cut -f 8 > Toy_SE.5p.antisense.profile
bedtools coverage -d -a Toy_SE.5p+100.bed -b regions.bed | cut -f 8 > Toy_SE.5p+100.profile
bedtools coverage -d -a Toy_SE.5p+MAPQ40.bed -b regions.bed | cut -f 8 > Toy_SE.5p+MAPQ40.profile
echo -e "#chrom\tstart\tend\tname\tscore\tstrand\tFivePrime_profile\tFivePrime_profile_sense\tFivePrime_profile_antisense\tFivePrime_profile_100_shift\tFivePrime_profile_min40qual" > Toy_SE.profile
paste Toy_SE.5p.profile Toy_SE.5p.sense.profile Toy_SE.5p.antisense.profile Toy_SE.5p+100.profile Toy_SE.5p+MAPQ40.profile >> Toy_SE.profile

# Toy_SE.bam depth
samtools view Toy_SE.bam | awk '$1==0 {printf "%s\t%s\t%s\t.\t%s\t+\n", $2, $3-1, $3+$5-1, $4}; $1==16 {printf "%s\t%s\t%s\t.\t%s\t-\n", $2, $3-1, $3+$5-1, $4}' > Toy_SE.Fragments.bed 
samtools view -q 40 Toy_SE.bam | awk '$1==0 {printf "%s\t%s\t%s\t.\t%s\t+\n", $2, $3-1, $3+$5-1, $4}; $1==16 {printf "%s\t%s\t%s\t.\t%s\t-\n", $2, $3-1, $3+$5-1, $4}' > Toy_SE.Fragments+MAPQ40.bed 
bedtools coverage -d -a Toy_SE.Fragments.bed -b regions.bed | cut -f1,2,3,4,5,6,8 > Toy_SE.Fragments.profile
bedtools coverage -d -a Toy_SE.Fragments+MAPQ40.bed -b regions.bed | cut -f 8 > Toy_SE.Fragments+MAPQ40.profile
echo -e "#chrom\tstart\tend\tname\tscore\tstrand\tDepth\tDepth_min40qual" > Toy_SE.depth
paste Toy_SE.Fragments.profile Toy_SE.Fragments+MAPQ40.profile >> Toy_SE.depth

# Toy_PE.bam 5'-end counting
# use only first in properly paired pair "-f 66"
samtools view -f 66 Toy_PE.bam | awk '$1==99 {printf "%s\t%s\t%s\t.\t%s\t+\n", $2, $3-1, $3, $4}; $1==83 {printf "%s\t%s\t%s\t.\t%s\t-\n", $2, $3-1+$5-1, $3+$5-1, $4}' > Toy_PE.5p.bed 
samtools view -f 66 Toy_PE.bam | awk '$1==99 {printf "%s\t%s\t%s\t.\t%s\t+\n", $2, $3-1+100, $3+100, $4}; $1==83 {printf "%s\t%s\t%s\t.\t%s\t-\n", $2, $3-1+$5-1-100, $3+$5-1-100, $4}' > Toy_PE.5p+100.bed 
samtools view -f 66 -q 40 Toy_PE.bam | awk '$1==99 {printf "%s\t%s\t%s\t.\t%s\t+\n", $2, $3-1, $3, $4}; $1==83 {printf "%s\t%s\t%s\t.\t%s\t-\n", $2, $3-1+$5-1, $3+$5-1, $4}' > Toy_PE.5p+MAPQ40.bed 
bedtools coverage -a Toy_PE.5p.bed -b regions.bed | cut -f 7 > 5p.all
bedtools coverage -s -a Toy_PE.5p.bed -b regions.bed | cut -f 7 > 5p.sense
bedtools coverage -S -a Toy_PE.5p.bed -b regions.bed | cut -f 7 > 5p.antisense
bedtools coverage -a Toy_PE.5p+100.bed -b regions.bed | cut -f 7 > 5p.shift
bedtools coverage -a Toy_PE.5p+MAPQ40.bed -b regions.bed | cut -f 7 > 5p.qual

# Toy_PE.bam Fragment middle point counting
# use only first in properly paired pair "-f 66"
samtools view -f 66 Toy_PE.bam | awk '$1==99 {printf "%s\t%s\t%s\t.\t%s\t+\n", $2, $3-1+int($8/2), $3+int($8/2), $4}; $1==83 {printf "%s\t%s\t%s\t.\t%s\t-\n", $2, $3-1+$5-1+int($8/2), $3+$5-1+int($8/2), $4}' > Toy_PE.MidPoint.bed 
samtools view -f 66 Toy_PE.bam | awk '$1==99 {printf "%s\t%s\t%s\t.\t%s\t+\n", $2, $3-1+int($8/2)+100, $3+int($8/2)+100, $4}; $1==83 {printf "%s\t%s\t%s\t.\t%s\t-\n", $2, $3-1+$5-1+int($8/2)-100, $3+$5-1+int($8/2)-100, $4}' > Toy_PE.MidPoint+100.bed 
samtools view -f 66 -q 40 Toy_PE.bam | awk '$1==99 {printf "%s\t%s\t%s\t.\t%s\t+\n", $2, $3-1+int($8/2), $3+int($8/2), $4}; $1==83 {printf "%s\t%s\t%s\t.\t%s\t-\n", $2, $3-1+$5-1+int($8/2), $3+$5-1+int($8/2), $4}' > Toy_PE.MidPoint+MAPQ40.bed 
bedtools coverage -a Toy_PE.MidPoint.bed -b regions.bed | cut -f 7 > MidPoint.all
bedtools coverage -s -a Toy_PE.MidPoint.bed -b regions.bed | cut -f 7 > MidPoint.sense
bedtools coverage -S -a Toy_PE.MidPoint.bed -b regions.bed | cut -f 7 > MidPoint.antisense
bedtools coverage -a Toy_PE.MidPoint+100.bed -b regions.bed | cut -f 7 > MidPoint.shift
bedtools coverage -a Toy_PE.MidPoint+MAPQ40.bed -b regions.bed | cut -f 7 > MidPoint.qual
echo -e "#chrom\tstart\tend\tname\tscore\tstrand\tFivePrime_count\tFivePrime_count_sense\tFivePrime_count_antisense\tFivePrime_count_100_shift\tFivePrime_count_min40qual\tMidPoint_count\tMidPoint_count_sense\tMidPoint_count_antisense\tMidPoint_count_100_shift\tMidPoint_count_min40qual" > Toy_PE.counts
paste regions.bed 5p.all 5p.sense 5p.antisense 5p.shift 5p.qual MidPoint.all MidPoint.sense MidPoint.antisense MidPoint.shift MidPoint.qual >> Toy_PE.counts

# Toy_PE.bam profiles
bedtools coverage -d -a Toy_PE.5p.bed -b regions.bed | cut -f1,2,3,4,5,6,8 > Toy_PE.5p.profile
bedtools coverage -d -s -a Toy_PE.5p.bed -b regions.bed | cut -f 8 > Toy_PE.5p.sense.profile
bedtools coverage -d -S -a Toy_PE.5p.bed -b regions.bed | cut -f 8 > Toy_PE.5p.antisense.profile
bedtools coverage -d -a Toy_PE.5p+100.bed -b regions.bed | cut -f 8 > Toy_PE.5p+100.profile
bedtools coverage -d -a Toy_PE.5p+MAPQ40.bed -b regions.bed | cut -f 8 > Toy_PE.5p+MAPQ40.profile
bedtools coverage -d -a Toy_PE.MidPoint.bed -b regions.bed | cut -f 8 > Toy_PE.MidPoint.profile
bedtools coverage -d -s -a Toy_PE.MidPoint.bed -b regions.bed | cut -f 8 > Toy_PE.MidPoint.sense.profile
bedtools coverage -d -S -a Toy_PE.MidPoint.bed -b regions.bed | cut -f 8 > Toy_PE.MidPoint.antisense.profile
bedtools coverage -d -a Toy_PE.MidPoint+100.bed -b regions.bed | cut -f 8 > Toy_PE.MidPoint+100.profile
bedtools coverage -d -a Toy_PE.MidPoint+MAPQ40.bed -b regions.bed | cut -f 8 > Toy_PE.MidPoint+MAPQ40.profile
echo -e "#chrom\tstart\tend\tname\tscore\tstrand\tFivePrime_profile\tFivePrime_profile_sense\tFivePrime_profile_antisense\tFivePrime_profile_100_shift\tFivePrime_profile_min40qual\tMidPoint_profile\tMidPoint_profile_sense\tMidPoint_profile_antisense\tMidPoint_profile_100_shift\tMidPoint_profile_min40qual" > Toy_PE.profile
paste Toy_PE.5p.profile Toy_PE.5p.sense.profile Toy_PE.5p.antisense.profile Toy_PE.5p+100.profile Toy_PE.5p+MAPQ40.profile Toy_PE.MidPoint.profile Toy_PE.MidPoint.sense.profile Toy_PE.MidPoint.antisense.profile Toy_PE.MidPoint+100.profile Toy_PE.MidPoint+MAPQ40.profile >> Toy_PE.profile

# Toy_PE.bam depth
samtools view -f 66 Toy_PE.bam | awk '($1==99) {printf "%s\t%s\t%s\t.\t%s\t+\n", $2, $3-1, $3+($8<0?$8*(-1):$8)-1, $4}; ($1==83&&$8<=0) {printf "%s\t%s\t%s\t.\t%s\t-\n", $2, $3-1+$5-($8<0?$8*(-1):$8)-1, $3+$5-1, $4}' > Toy_PE.Fragments.bed 
samtools view -f 66 -q 40 Toy_PE.bam | awk '($1==99) {printf "%s\t%s\t%s\t.\t%s\t+\n", $2, $3-1, $3+($8<0?$8*(-1):$8)-1, $4}; ($1==83&&$8<=0) {printf "%s\t%s\t%s\t.\t%s\t-\n", $2, $3-1+$5-($8<0?$8*(-1):$8)-1, $3+$5-1, $4}' > Toy_PE.Fragments+MAPQ40.bed 
bedtools coverage -d -a Toy_PE.Fragments.bed -b regions.bed | cut -f1,2,3,4,5,6,8 > Toy_PE.Fragments.profile
bedtools coverage -d -a Toy_PE.Fragments+MAPQ40.bed -b regions.bed | cut -f 8 > Toy_PE.Fragments+MAPQ40.profile
echo -e "#chrom\tstart\tend\tname\tscore\tstrand\tDepth\tDepth_min40qual" > Toy_PE.depth
paste Toy_PE.Fragments.profile Toy_PE.Fragments+MAPQ40.profile >> Toy_PE.depth

# Cleanup
rm -rf Toy_SE.5p* Toy_PE.5p* Toy_PE.MidPoint* 5p.* MidPoint* Toy_PE.Fragments* Toy_SE.Fragments*
