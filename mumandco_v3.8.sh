#!/bin/bash
set -euo pipefail

version="3.8"

##################################################################################
############################## EDIT PATHS AND NAMES ##############################
##################################################################################

### Paths to Additional tools/scripts and check if they are installed and found in path ###
NUCMER=$(which nucmer)
[[ $NUCMER == "" ]] && echo "ERROR: Cannot find nucmer script using 'which nucmer', make sure MUMmer (=>V4) is installed and in path" && exit
DELTAFILTER=$(which delta-filter)
[[ $NUCMER == "" ]] && echo "ERROR: Cannot find delta-filter script using 'which delta-filter', make sure MUMmer (=>V4) is installed and in path" && exit
SHOWCOORDS=$(which show-coords)
[[ $NUCMER == "" ]] && echo "ERROR: Cannot find dnadiff script using 'which dnadiff', make sure MUMmer (=>V4) is installed and in path" && exit

##Genomes for alignments##
#reference_assembly="./yeast.tidy.fa"
#genome_size=12500000
#query_assembly="./yeast_tidy_DEL100.fa"

##output file names for alignment and filtering##
#ref="yeast_tidy"
#query="DEL100"
#alignments_folder="DEL100_alignments"

##output file name for SV detection##
#prefix="DEL100_test_v2.1"

##if you want to remove intemediate files at the end put "yes"##
cleanup="yes"

##subtelomeric regions file. HANGOVER FROM EARLIER VERSIONS. ESSENTIALLY DOES NOTHING NOW.
filter_subtelomeric_region="no"
#subtelo_coords="subtelomeric_regions_file.csv"
rDNA_filter="no"


#default values, unless denoted when running MUM&Co
reference_assembly=""
query_assembly="" 
genome_size=""
prefix="mumandco"
threads="1"
minlen="50"
blast_step="no"



while [[ $# -gt 0 ]]
do
key="$1"

case "$key" in
	-r|--reference_genome)
	reference_assembly="$2"
	shift
	shift
	;;
	-q|--query_genome)
	query_assembly="$2"
	shift
	shift
	;;
	-g|--genome_size)
	genome_size="$2"
	shift
	shift
	;;
	-o|--output)
	prefix="$2"
	shift
	shift
	;;
	-t|--threads)
	threads="$2"
	shift
	shift
	;;
	-ml|--minlen)
	minlen="$2"
	shift
	shift
	;;
	-b|--blast)
	blast_step="yes"
	shift
	;;
	esac
done


#creates error message and exits if these values are not assigned 
[[ $reference_assembly == "" ]] && echo "ERROR: Path to reference genome not found, assign using -r" && exit
[[ $query_assembly == "" ]] && echo "ERROR: Path to query genome not found, assign using -q" && exit
[[ $genome_size == "" ]] && echo "ERROR: genome size not found, assign using -g" && exit
#uses default 'mumandco' prefix for output and creates error and exits if output directory already exists
[[ $prefix == "mumandco" ]] && echo "WARNING: No option for output found, using mumandco"
[[ -d $prefix"_output" ]] && echo "ERROR: Output directory already exists, please remove or set alternate output" && exit

#if blast option is set, check to see if samtools and blast are installed and in path. Exit if not found
if [[ ${blast_step} == "yes" ]]
then
	SAMTOOLS=$(which samtools)
	[[ $SAMTOOLS == "" ]] && echo "ERROR: Cannot find samtools using 'which samtools', make sure samtools is installed and in path" && exit
	BLASTN=$(which blastn)
	[[ $BLASTN == "" ]] && echo "ERROR: Cannot find blastn script using 'which blastn', make sure BLAST is installed and in path" && exit
	BLASTDB=$(which makeblastdb)
	[[ $BLASTDB == "" ]] && echo "ERROR: Cannot find blastdb script using 'which blastdb', make sure BLAST is installed and in path" && exit
fi



###################################################################################
###################################################################################
###################################################################################

####align genomes to one another###
echo ""
echo "Nucmer alignment of genomes, filtering and converting to coordinates"
echo ""

if [ $genome_size -le 100000000 ]
then
	$NUCMER --threads ${threads} --maxmatch --nosimplify -p ""$prefix"_ref" $reference_assembly $query_assembly
	$DELTAFILTER -m ""$prefix"_ref".delta > ""$prefix"_ref".delta_filter
	$SHOWCOORDS -T -r -c -l -d -g ""$prefix"_ref".delta_filter > ""$prefix"_ref".delta_filter.coordsg
	$SHOWCOORDS -T -r -c -l -d ""$prefix"_ref".delta_filter > ""$prefix"_ref".delta_filter.coords

	$NUCMER --threads ${threads} --maxmatch --nosimplify -p ""$prefix"_query" $query_assembly $reference_assembly
	$DELTAFILTER -m ""$prefix"_query".delta > ""$prefix"_query".delta_filter
	$SHOWCOORDS -T -r -c -l -d -g ""$prefix"_query".delta_filter > ""$prefix"_query".delta_filter.coordsg
	$SHOWCOORDS -T -r -c -l -d ""$prefix"_query".delta_filter > ""$prefix"_query".delta_filter.coords
else

	if [ $genome_size -le 500000000 ]
	then
		$NUCMER --threads ${threads} --maxmatch -l 100 -c 500 -p ""$prefix"_ref" $reference_assembly $query_assembly
		$DELTAFILTER -m ""$prefix"_ref".delta > ""$prefix"_ref".delta_filter
		$SHOWCOORDS -T -r -c -l -d -g ""$prefix"_ref".delta_filter > ""$prefix"_ref".delta_filter.coordsg
		$SHOWCOORDS -T -r -c -l -d ""$prefix"_ref".delta_filter > ""$prefix"_ref".delta_filter.coords
	
		$NUCMER --threads ${threads} --maxmatch -l 100 -c 500 -p ""$prefix"_query" $query_assembly $reference_assembly
		$DELTAFILTER -m ""$prefix"_query".delta > ""$prefix"_query".delta_filter
		$SHOWCOORDS -T -r -c -l -d -g ""$prefix"_query".delta_filter > ""$prefix"_query".delta_filter.coordsg
		$SHOWCOORDS -T -r -c -l -d ""$prefix"_query".delta_filter > ""$prefix"_query".delta_filter.coords
	else
		echo "My what a large genome you have, this may take some time"
		$NUCMER --threads ${threads} --maxmatch -l 500 -c 500 -p ""$prefix"_ref" $reference_assembly $query_assembly
		$DELTAFILTER -m ""$prefix"_ref".delta > ""$prefix"_ref".delta_filter
		$SHOWCOORDS -T -r -c -l -d -g ""$prefix"_ref".delta_filter > ""$prefix"_ref".delta_filter.coordsg
		$SHOWCOORDS -T -r -c -l -d ""$prefix"_ref".delta_filter > ""$prefix"_ref".delta_filter.coords
	
		$NUCMER --threads ${threads} --maxmatch -l 500 -c 500 -p ""$prefix"_query" $query_assembly $reference_assembly
		$DELTAFILTER -m ""$prefix"_query".delta > ""$prefix"_query".delta_filter
		$SHOWCOORDS -T -r -c -l -d -g ""$prefix"_query".delta_filter > ""$prefix"_query".delta_filter.coordsg
		$SHOWCOORDS -T -r -c -l -d ""$prefix"_query".delta_filter > ""$prefix"_query".delta_filter.coords
	fi
fi

#!/bin/bash

mkdir "$prefix"_alignments
alignments_folder=""$prefix"_alignments"
mv ""$prefix"_ref".delta* $alignments_folder/
mv ""$prefix"_query".delta* $alignments_folder/

echo "                               #"
echo "                              # #"
echo "                             #   #"
echo "                ###############################"
echo "                #                             #"
echo "                # MUM&Co is open for business #"
echo "                #           version ${version}       #"
echo "                ###############################"

############################################################################################################
################################DELETIONS, INSERTIONS AND TRANSLOCATIONS####################################
############################################################################################################
                         ##USES SHOW-COORDS WITH G OPTION FOR GLOBAL ALIGNMENTS##

echo ""
echo ""
echo "######################################################################################################"
echo "          USING GLOBAL ALIGNMENT COORDINATES FOR DELETIONS, INSERTIONS AND TRANSLOCATIONS"
echo "######################################################################################################"
echo ""
###Begin Filtering###
##filter out chrMT and aligned contigs that don't match the primary matched contigs.
echo ""
echo "Matching query and reference chromosomes"
echo ""

####How to say, from row 1 to row 14 ?? and not having to write it all out seperately####
#tail -n +5 "$alignments_folder/""$prefix"_ref".delta_filter.coordsg" | sed 's/_/\t/g' |\
#													awk '!/chrMT/{for(i=15;i<=NF;++i) if($14==$i) {print $0}}' |\
#													sed 's/\t/_/15g' \
#													>""$prefix"_ref".coordsg_matched

#tail -n +5 "$alignments_folder/""$prefix"_query".delta_filter.coordsg" | \
#													awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$15"\t"$14}' |\
#													sed 's/_/\t/g' |\
#													awk '!/chrMT/{for(i=15; i<=NF; ++i) if($i==$14) {print $0}}' |\
#													sed 's/\t/_/15g' |\
#													awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$15"\t"$14}' \
#													> ""$prefix"_query".coordsg_matched

##get average quality
average_quality=$(tail -n +5 "$alignments_folder/""$prefix"_ref".delta_filter.coordsg" | awk '{sum+=$7} END{print sum/NR}')

### filter for alignment less than 10kb for the chromosome pairing (therefore the limit on translocation size)
### filter for any alignment overlapping another with more than 50bp and with lower quality
tail -n +5 "$alignments_folder/""$prefix"_ref".delta_filter.coordsg" |\
	sort -k15,15d -k3,3n |\
	awk '{if($5>10000) print $0}' |\
		awk 'BEGIN{qual=0};\
			{if(qual==0 && $13=="1") {print $0; qual=$7; start=$3; end=$4; contig=$15}\
			else if(qual==0 && $13=="-1") {print $0; qual =$7; start=$4; end=$3; contig=$15}\
			else if(contig!=$15 && $13=="1") {print $0; qual=$7; start=$3; end=$4; contig=$15}\
			else if(contig!=$15 && $13=="-1") {print $0; qual=$7; start=$4; end=$3; contig=$15}\
			else if(contig==$15 && $13=="1" && $3 >=start+50 && $3 <=end-50 && $7 <=qual || contig==$15 && $13=="1" && $4 >=start+50 && $4 <=end-50 && $7 <=qual) {}\
			else if(contig==$15 && $13=="-1" && $3 >=start+50 && $3 <=end-50 && $7 <=qual || contig==$15 && $13=="-1" && $4 >=start+50 && $4 <=end-50 && $7 <=qual) {}\
			else if($13=="1") {print $0; qual=$7; start=$3; end=$4; contig=$15}\
			else if($13=="-1") {print $0; qual=$7; start=$4; end=$3; contig=$15}}'|\
				sort -k15,15d -k3,3nr |\
					awk 'BEGIN{qual=0};\
						{if(qual==0 && $13=="1") {print $0; qual=$7; start=$3; end=$4; contig=$15}\
						else if(qual==0 && $13=="-1") {print $0; qual =$7; start=$4; end=$3; contig=$15}\
						else if(contig!=$15 && $13=="1") {print $0; qual=$7; start=$3; end=$4; contig=$15}\
						else if(contig!=$15 && $13=="-1") {print $0; qual=$7; start=$4; end=$3; contig=$15}\
						else if(contig==$15 && $13=="1" && $3 >=start+50 && $3 <=end-50 && $7 <=qual || contig==$15 && $13=="1" && $4 >=start+50 && $4 <=end-50 && $7 <=qual) {}\
						else if(contig==$15 && $13=="-1" && $3 >=start+50 && $3 <=end-50 && $7 <=qual || contig==$15 && $13=="-1" && $4 >=start+50 && $4 <=end-50 && $7 <=qual) {}\
						else if($13=="1") {print $0; qual=$7; start=$3; end=$4; contig=$15}\
						else if($13=="-1") {print $0; qual=$7; start=$4; end=$3; contig=$15}}'|\
				sort -k14,14d -k1,1n |\
					awk 'BEGIN{qual=0};\
						{if(qual==0 && $12=="1") {print $0; qual=$7; start=$1; end=$2; contig=$14}\
						else if(qual==0 && $12=="-1") {print $0; qual =$7; start=$2; end=$1; contig=$14}\
						else if(contig!=$14 && $12=="1") {print $0; qual=$7; start=$1; end=$2; contig=$14}\
						else if(contig!=$14 && $12=="-1") {print $0; qual=$7; start=$2; end=$1; contig=$14}\
						else if(contig==$14 && $12=="1" && $1 >=start+50 && $1 <=end-50 && $7 <=qual || contig==$14 && $12=="1" && $2 >=start+50 && $2 <=end-50 && $7 <=qual) {}\
						else if(contig==$14 && $12=="-1" && $1 >=start+50 && $1 <=end-50 && $7 <=qual || contig==$14 && $12=="-1" && $2 >=start+50 && $2 <=end-50 && $7 <=qual) {}\
						else if($12=="1") {print $0; qual=$7; start=$1; end=$2; contig=$14}\
						else if($12=="-1") {print $0; qual=$7; start=$2; end=$1; contig=$14}}' |\
							sort -k14,14d -k1,1nr |\
								awk 'BEGIN{qual=0};\
									{if(qual==0 && $12=="1") {print $0; qual=$7; start=$1; end=$2; contig=$14}\
									else if(qual==0 && $12=="-1") {print $0; qual =$7; start=$2; end=$1; contig=$14}\
									else if(contig!=$14 && $12=="1") {print $0; qual=$7; start=$1; end=$2; contig=$14}\
									else if(contig!=$14 && $12=="-1") {print $0; qual=$7; start=$2; end=$1; contig=$14}\
									else if(contig==$14 && $12=="1" && $1 >=start+50 && $1 <=end-50 && $7 <=qual || contig==$14 && $12=="1" && $2 >=start+50 && $2 <=end-50 && $7 <=qual) {}\
									else if(contig==$14 && $12=="-1" && $1 >=start+50 && $1 <=end-50 && $7 <=qual || contig==$14 && $12=="-1" && $2 >=start+50 && $2 <=end-50 && $7 <=qual) {}\
									else if($12=="1") {print $0; qual=$7; start=$1; end=$2; contig=$14}\
									else if($12=="-1") {print $0; qual=$7; start=$2; end=$1; contig=$14}}' > ""$prefix"_pre_list.txt"

#tail -n +5 "$alignments_folder/""$prefix"_ref".delta_filter.coordsg" |\
#	sort -k15,15d -k3,3n |\
#	awk '{if($5>10000) print $0}' |\
#		awk 'BEGIN{qual=0};\
#			{if(qual==0) {print $0; qual=$7; mid=($3+$4)/2; contig=$15}\
#			else if(contig!=$15) {print $0; qual=$7; mid=($3+$4)/2; contig=$15}\
#			else if(contig==$15 && $3>$4 && $3 >=mid && $4 <=mid && $7 <qual) {}\
#			else if(contig==$15 && $4>$3 && $4 >=mid && $3 <=mid && $7 <qual) {}\
#			else {print $0; qual=$7; mid=($3+$4)/2; contig=$15}}'|\
#				sort -k15,15d -k3,3nr |\
#					awk 'BEGIN{qual=0};\
#						{if(qual==0) {print $0; qual=$7; mid=($1+$2)/2; contig=$15}\
#						else if(contig!=$15) {print $0; qual=$7; mid=($1+$2)/2; contig=$15}\
#						else if(contig==$15 && $1>$2 && $1 >=mid && $2 <=mid && $7 <qual) {}\
#						else if(contig==$15 && $2>$1 && $2 >=mid && $1 <=mid && $7 <qual) {}\
#						else {print $0; qual=$7; mid=($1+$2)/2; contig=$15}}'|\
#							sort -k14,14d -k1,1n |\
#								awk 'BEGIN{qual=0};\
#									{if(qual==0) {print $0; qual=$7; mid=($1+$2)/2; contig=$14}\
#									else if(contig!=$14) {print $0; qual=$7; mid=($1+$2)/2; contig=$14}\
#									else if(contig==$14 && $1>$2 && $1 >=mid && $2 <=mid && $7 <qual) {}\
#									else if(contig==$14 && $2>$1 && $2 >=mid && $1 <=mid && $7 <qual) {}\
#									else {print $0; qual=$7; mid=($1+$2)/2; contig=$14}}'|\
#										sort -k14,14d -k1,1nr |\
#											awk 'BEGIN{qual=0};\
#												{if(qual==0) {print $0; qual=$7; mid=($1+$2)/2; contig=$14}\
#												else if(contig!=$14) {print $0; qual=$7; mid=($1+$2)/2; contig=$14}\
#												else if(contig==$14 && $1>$2 && $1 >=mid && $2 <=mid && $7 <qual) {}\
#												else if(contig==$14 && $2>$1 && $2 >=mid && $1 <=mid && $7 <qual) {}\
#												else {print $0; qual=$7; mid=($1+$2)/2; contig=$14}}' > ""$prefix"_pre_list.txt"
#
#
#

##filter any remaining alignments with less than (1.25 x standard deviation) less than the average of the remaining alignments
##actually filter out alignments that are smaller than 1kb after the chromosome pairing (helps with less contiguous query genomes)
average_quality=$(cat ""$prefix"_pre_list.txt" | awk '{sum+=$7} END{print sum/NR}')

if [[ "$average_quality" == "100" ]]
	then

		cat ""$prefix"_pre_list.txt" | awk  '{print $14"\t"$15}' | sort | uniq > ""$prefix"_chromosome_pairs.txt"
	
		awk 'NR==FNR{a[$1,$2]++;next};a[$14,$15] {print $0}' ""$prefix"_chromosome_pairs.txt" "$alignments_folder/""$prefix"_ref".delta_filter.coordsg" | awk -v minlen="$minlen" '{if($5 > minlen) print}' > ""$prefix"_ref".coordsg_matched
		awk 'NR==FNR{a[$1,$2]++;next};a[$15,$14] {print $0}' ""$prefix"_chromosome_pairs.txt" "$alignments_folder/""$prefix"_query".delta_filter.coordsg" | awk -v minlen="$minlen" '{if($5 > minlen) print}' > ""$prefix"_query".coordsg_matched

	else

		stdev125=$(cat ""$prefix"_pre_list.txt" | awk '{sum=sum+$7 ; sumX2+=(($7)^2)} END{print 1.25*(sqrt(sumX2/(NR) - ((sum/NR)^2)))}')
		quality_filt=$(echo $average_quality - $stdev125 | bc)
		cat ""$prefix"_pre_list.txt" | awk -v quality_filt="$quality_filt" '{if($7 > quality_filt) print $14"\t"$15}' | sort | uniq > ""$prefix"_chromosome_pairs.txt"
	
		awk 'NR==FNR{a[$1,$2]++;next};a[$14,$15] {print $0}' ""$prefix"_chromosome_pairs.txt" "$alignments_folder/""$prefix"_ref".delta_filter.coordsg" | awk -v minlen="$minlen" '{if($5 > minlen) print}' > ""$prefix"_ref".coordsg_matched
		awk 'NR==FNR{a[$1,$2]++;next};a[$15,$14] {print $0}' ""$prefix"_chromosome_pairs.txt" "$alignments_folder/""$prefix"_query".delta_filter.coordsg" | awk -v minlen="$minlen" '{if($5 > minlen) print}' > ""$prefix"_query".coordsg_matched
fi


if [[ $filter_subtelomeric_region == "yes" ]]
then

	##All those calls that start before the end of the subtelomeric region or finish after are filtered
	echo ""
	echo "Removing calls in subtelomeric regions"
	echo ""
	
	awk 'FNR==NR{a[$1]=$2"\t"$3;next} ($14 in a) {print $0"\t"a[$14]}' $subtelo_coords ""$prefix"_ref".coordsg_matched  |\
		awk '{if($2 > $16 && $1 < $17) {print $0}}' |\
			awk -F'\t' 'BEGIN { OFS = FS }; NF{NF-=2};1' > ""$prefix"_ref".coordsg_matched_subtelofilt
	
	awk 'FNR==NR{a[$1]=$2"\t"$3;next} ($15 in a) {print $0"\t"a[$15]}' $subtelo_coords ""$prefix"_query".coordsg_matched  |\
		awk '{if($4 > $16 && $3 < $17) {print $0}}' |\
			awk -F'\t' 'BEGIN { OFS = FS }; NF{NF-=2};1' > ""$prefix"_query".coordsg_matched_subtelofilt

	#################################################################################
	##########getting gaps from subtelomeres to beginning of alignments##############
	#################################################################################
	
	awk 'FNR==NR{a[$1]=$2"\t"$3;next} ($14 in a) {print $0"\t"a[$14]}' $subtelo_coords ""$prefix"_ref".coordsg_matched_subtelofilt |\
		awk 'BEGIN{chrr=0} \
				{if(chrr==0 && $1>$16){print $0; chrr=$14} \
				{if(chrr==0 && $1<$16){chrr=$14} \
				{if(chrr!=$14 && $1>$16) {print $0; chrr=$14}\
				{if(chrr!=$14 && $1<$16) {chrr=$14}}}}}' |\
						awk '{if($13==1){print $14"\t"$15"\t1\t"$1"\t"$1-1"\tdeletion\t1\t"$3} \
							{if($13!=1){print $14"\t"$15"\t1\t"$1"\t"$1-1"\tdeletion\t"$3"\t"$9}}}' \
				> ""$prefix"_ref".coordsg_matched_startgap
	
	
	awk 'FNR==NR{a[$1]=$2"\t"$3;next} ($14 in a) {print $0"\t"a[$14]}' $subtelo_coords ""$prefix"_ref".coordsg_matched_subtelofilt \
			> ""$prefix"_ref".coordsg_matched_endgap_temp
		
	tac ""$prefix"_ref".coordsg_matched_endgap_temp |\
		awk 'BEGIN{chrr=0} \
				{if(chrr==0 && $2<$17){print $0; chrr=$14} \
				{if(chrr==0 && $2>$17){chrr=$14} \
				{if(chrr!=$14 && $2<$17) {print $0; chrr=$14}\
				{if(chrr!=$14 && $2>$17) {chrr=$14}}}}}' |\
					awk '{if($13==1){print $14"\t"$15"\t"$2"\t"$8"\t"$8-$2"\tdeletion\t"$4"\t"$9} \
							{if($13!=1){print $14"\t"$15"\t"$2"\t"$8"\t"$8-$2"\tdeletion\t1\t"$4}}}' \
							> ""$prefix"_ref".coordsg_matched_endgap
	
	cat ""$prefix"_ref".coordsg_matched_startgap ""$prefix"_ref".coordsg_matched_endgap > ""$prefix"_ref".coordsg_matched_subtelogap
	
	
	awk 'FNR==NR{a[$1]=$2"\t"$3;next} ($14 in a) {print $0"\t"a[$14]}' $subtelo_coords ""$prefix"_query".coordsg_matched_subtelofilt |\
		awk 'BEGIN{chrr=0} \
				{if(chrr==0 && $1>$16){print $0; chrr=$14} \
				{if(chrr==0 && $1<$16){chrr=$14} \
				{if(chrr!=$14 && $1>$16) {print $0; chrr=$14}\
				{if(chrr!=$14 && $1<$16) {chrr=$14}}}}}' |\
					awk '{if($13==1){print $15"\t"$14"\t1\t"$3"\t"$3-1"\tinsertion\t1\t"$1} \
							{if($13!=1){print $15"\t"$14"\t1\t"$4"\t"$4-1"\tinsertion\t"$2"\t"$9}}}' \
				> ""$prefix"_query".coordsg_matched_startgap
	
	
	awk 'FNR==NR{a[$1]=$2"\t"$3;next} ($14 in a) {print $0"\t"a[$14]}' $subtelo_coords ""$prefix"_query".coordsg_matched_subtelofilt \
			> ""$prefix"_query".coordsg_matched_endgap_temp
		
	tac ""$prefix"_query".coordsg_matched_endgap_temp |\
		awk 'BEGIN{chrr=0} \
				{if(chrr==0 && $2<$17){print $0; chrr=$14} \
				{if(chrr==0 && $2>$17){chrr=$14} \
				{if(chrr!=$14 && $2<$17) {print $0; chrr=$14}\
				{if(chrr!=$14 && $2>$17) {chrr=$14}}}}}' |\
					awk '{if($13==1){print $15"\t"$14"\t"$4"\t"$8"\t"$8-$4"\tinsertion\t"$2"\t"$9} \
							{if($13!=1){print $15"\t"$14"\t"$3"\t"$8"\t"$8-$3"\tinsertion\t1\t"$4}}}' \
							> ""$prefix"_query".coordsg_matched_endgap
	
	cat ""$prefix"_query".coordsg_matched_startgap ""$prefix"_query".coordsg_matched_endgap > ""$prefix"_query".coordsg_matched_subtelogap
	
	cat ""$prefix"_ref".coordsg_matched_subtelogap ""$prefix"_query".coordsg_matched_subtelogap |\
		awk '{if($5>50) print $0}' > $prefix.coordsg_subtelogaps_50

	cat $prefix.coordsg_subtelogaps_50 |\
		awk '{if($5>1000) print $0}' > $prefix.coordsg_subtelogaps_1000

	rm ""$prefix"_ref".coordsg_matched_endgap_temp
	rm ""$prefix"_ref".coordsg_matched_startgap
	rm ""$prefix"_ref".coordsg_matched_endgap
	rm ""$prefix"_query".coordsg_matched_endgap_temp
	rm ""$prefix"_query".coordsg_matched_startgap
	rm ""$prefix"_query".coordsg_matched_endgap
else 
	cp ""$prefix"_ref".coordsg_matched ""$prefix"_ref".coordsg_matched_subtelofilt
	cp ""$prefix"_query".coordsg_matched ""$prefix"_query".coordsg_matched_subtelofilt
fi

#############################################
###################INDELS####################
#############################################

echo ""
echo "Finding alignment gaps"
echo ""

cat ""$prefix"_ref".coordsg_matched_subtelofilt |\
				awk 'BEGIN{startref=0; chr=""; gap=0; startq=0} \
							{if(chr==$15) {gap=$1-startref; print $0"\t"startref"\t"$1"\t"gap"\t"startq"\t"$3; startref=$2; startq=$4} \
							if(chr!=$15) {print $0; chr=$15; startref=$2; startq=$4}}' \
								> ""$prefix"_ref".coordsg_matched_subtelofilt_gaps

cat ""$prefix"_query".coordsg_matched_subtelofilt |\
				awk 'BEGIN{startq=0; chr=""; gap=0; startref=0} \
							{if(chr==$14) {gap=$1-startq; print $0"\t"startref"\t"$3"\t"gap"\t"startq"\t"$1; startq=$2; startref=$4; startq=$2} \
							if(chr!=$14) {print $0; chr=$14; startq=$2; startref=$4}}' \
								> ""$prefix"_query".coordsg_matched_subtelofilt_gaps


###remove SVS smaller than 50 for 1000bp and remove unecessary columns###
echo ""
echo "Filtering for size labelling SV"
echo ""

cat ""$prefix"_ref".coordsg_matched_subtelofilt_gaps |\
awk '{if($18>50){print $14"\t"$15"\t"$16"\t"$17"\t"$18"\tdeletion\t"$19"\t"$20}}' |\
sort -k1,1 -k3V > ""$prefix"_ref".gaps_50

cat ""$prefix"_ref".gaps_50 | awk '{if($5>1000){print $0}}' > ""$prefix"_ref".gaps_1000

##again for reciprical and also swap chromosome names between ref and query, ref in far left column as in deletions##
cat ""$prefix"_query".coordsg_matched_subtelofilt_gaps |\
awk '{if($18>50){print $15"\t"$14"\t"$16"\t"$17"\t"$18"\tinsertion\t"$19"\t"$20}}' |\
sort -k1,1 -k3V > ""$prefix"_query".gaps_50

cat ""$prefix"_query".gaps_50 | awk '{if($5>1000){print $0}}' > ""$prefix"_query".gaps_1000

#############################################
###############TRANSLOCATIONS################
#############################################

echo ""
echo "Finding translocation fragments"
echo ""

###keep only query contigs that are associated with two or more reference chromosomes
cat ""$prefix"_ref".coordsg_matched_subtelofilt | awk '{print $15}' | sort -u > transloc_list.txt
cat transloc_list.txt | while read chromosome
do
	associated_chromosomes=$(cat ""$prefix"_ref".coordsg_matched_subtelofilt  | awk -v chromosome="$chromosome" '{if($15==chromosome) print $14}' |	sort -u | wc -l)
	if [ $associated_chromosomes -gt 1 ]
	then
		cat ""$prefix"_ref".coordsg_matched_subtelofilt | awk -v chromosome="$chromosome" '{if($15==chromosome) print $0}' >$prefix.$chromosome.transloc_pairing
	fi
done

transloc_true=$( find . -maxdepth 1 -type f -name '*.transloc_pairing' | wc -l ) 
if [ $transloc_true -gt 0 ]
then
	cat *.transloc_pairing | sort -k14,14 -k1V > $prefix.transloc_candidate_alignments
	rm *.transloc_pairing
	rm transloc_list.txt

	###remove smaller false translocations that lie inside the larger true fragment then concatenate again
	cat $prefix.transloc_candidate_alignments | awk 'BEGIN{chr=""; start=""; stop=""}; \
										{if($14==chr && $1>start && $2<stop) {} \
										else if($14==chr) {print $0; chr=$14; start=$1; stop=$2} \
										else {print $0; chr=$14; start=$1; stop=$2}} \
										{if(chr=="") {print $0; chr=$14; start=$1; stop=$2}}' \
										> $prefix.transloc_candidate_alignments2

	###need to fix for INVERTED SEGMENTS!!!!!!!!!!!!

	cat $prefix.transloc_candidate_alignments2 | awk 'BEGIN{chrr=""; chrq=""; startr=0; endr=0; startq=0; endq=0; link=0}\
												{if(chrq=="") {chrr=$14; chrq=$15; startr=$1; endr=$2; startq=$3; endq=$4; orientation=$13}\
												else if(chrr==$14 && chrq==$15) {endr=$2; endq=$4; orientation=$13}\
												else if(chrr==$14 && chrq!=$15) {print $0"\t"chrr"\t"chrq"\t"startr"\t"endr"\t"endr-startr"\ttransloc\t"startq"\t"endq"\t"orientation; startr=$1; endr=$2; startq=$3; endq=$4; chrq=$15; orientation=$13}\
												else if(chrr!=$14){print $0"\t"chrr"\t"chrq"\t"startr"\t"endr"\t"endr-startr"\ttransloc\t"startq"\t"endq"\t"orientation; startr=$1; endr=$2; startq=$3; endq=$4; chrr=$14; chrq=$15; orientation=$13}}\
													END{print $0"\t"chrr"\t"chrq"\t"startr"\t"endr"\t"endr-startr"\ttransloc\t"startq"\t"endq"\t"orientation}'|\
												awk '{print $16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22"\t"$23"\t"$24}' \
												> $prefix.transloc_candidate_alignments3_fragments


	##remove many small fragments and join the large fragments which they broke"
	cat $prefix.transloc_candidate_alignments3_fragments | awk '{if($5< 10000) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\ttransloc_10000"$6"\t"$7"\t"$8}' > ""$prefix"_ref".fragments_less10000
	cat $prefix.transloc_candidate_alignments3_fragments | awk '{if($5>=10000) print $0}'> $prefix.transloc_candidate_alignments4_fragments10000

	##AGAIN remove smaller false translocations that lie inside the larger true fragment then concatenate again
	cat $prefix.transloc_candidate_alignments4_fragments10000 | awk 'BEGIN{chr=""; start=""; stop=""}; \
											{if($1==chr && $3>start && $4<stop) {} \
											else if($1==chr) {print $0; chr=$1; start=$3; stop=$4} \
											else {print $0; chr=$1; start=$3; stop=$4}} \
											{if(chr=="") {print $0; chr=$1; start=$3; stop=$4}}' \
											> $prefix.transloc_candidate_alignments5_fragments10000

	cat $prefix.transloc_candidate_alignments5_fragments10000 | awk 'BEGIN{chrr=""; chrq=""; startr=0; endr=0; startq=0; endq=0; link=0}\
													{if(chrq=="") {chrr=$1; chrq=$2; startr=$3; endr=$4; startq=$7; endq=$8; orientation=$9}\
													else if(chrr==$1 && chrq==$2) {endr=$4; endq=$8}\
													else if(chrr==$1 && chrq!=$2) {print $0"\t"chrr"\t"chrq"\t"startr"\t"endr"\t"endr-startr"\ttransloc\t"startq"\t"endq"\t"orientation; startr=$3; endr=$4; startq=$7; endq=$8; chrq=$2; orientation=$9}\
													else if(chrr!=$1){print $0"\t"chrr"\t"chrq"\t"startr"\t"endr"\t"endr-startr"\ttransloc\t"startq"\t"endq"\t"orientation; startr=$3; endr=$4; startq=$7; endq=$8; chrr=$1; chrq=$2; orientation=$9}}\
														END{print $0"\t"chrr"\t"chrq"\t"startr"\t"endr"\t"endr-startr"\ttransloc\t"startq"\t"endq"\t"orientation}'|\
													awk '{print $10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18}' \
													> $prefix.transloc_candidate_alignments6_fragments10000

	##this makes sure that at the end the orientation of the fragment will reverse the start and stop positions
	##this will be returned to the orignal orientation but only after the correct borders are identified downstream
	cat $prefix.transloc_candidate_alignments6_fragments10000 |\
		awk '{if($9 == "-1") {print $1"\t"$2"\t"$4"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8} else {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}}' \
		> $prefix.transloc_candidate_alignments7_fragments10000
	
	rm $prefix.transloc_candidate_alignments
	rm $prefix.transloc_candidate_alignments2
	rm $prefix.transloc_candidate_alignments3_fragments
	rm $prefix.transloc_candidate_alignments4_fragments10000
	rm $prefix.transloc_candidate_alignments5_fragments10000
	rm $prefix.transloc_candidate_alignments6_fragments10000

else
	touch $prefix.transloc_candidate_alignments7_fragments10000
	touch ${prefix}_ref.fragments_less10000
	rm transloc_list.txt
fi

#
#rows=$(wc -l ""$prefix"_ref".fragments_big | awk '{print $1}')
#cat ""$prefix"_ref".fragments_big | awk -v rowss="$rows" 'BEGIN{all=""; chrr=""; chrq=""; start=0; link=0; startq=0; endq=0} \
#									{if(chrr==$1 && chrq==$2 && rowss!=NR) {end=$4; endq=$8; link=1} \
#									{if(chrr==$1 && chrq!=$2 && link==0 && rowss!=NR) {print all; chrq=$2; start=$3; startq=$7; endq=$8; all=$0} \
#									{if(chrr!="" && chrr!=$1 && link==0 && rowss!=NR) {print all; chrr=$1; chrq=$2; start=$3; startq=$7; all=$0} \
#									{if(chrr==$1 && chrq!=$2 && link==1 && rowss!=NR) {print chrr"\t"chrq"\t"start"\t"end"\t"end-start"\ttransloc\t"startq"\t"endq; link=0; chrq=$2; start=$3; startq=$7; all=$0} \
#									{if(chrr!="" && chrr!=$1 && link==1 && rowss!=NR) {print chrr"\t"chrq"\t"start"\t"end"\t"end-start"\ttransloc\t"startq"\t"endq; link=0; chrr=$1; chrq=$2; start=$3; startq=$7; all=$0} \
#									{if(chrr=="") {chrr=$1; chrq=$2; start=$3; startq=$7; all=$0} \
#									{if(chrr==$1 && chrq!=$2 && link==1 && rowss==NR) {print chrr"\t"chrq"\t"start"\t"end"\t"end-start"\ttransloc\t"startq"\t"endq"\n"$0} \
#									{if(chrr==$1 && chrq==$2 && rowss==NR) {print chrr"\t"chrq"\t"start"\t"$4"\t"$4-start"\ttransloc\t"startq"\t"$8} \
#									{if(chrr==$1 && chrq!=$2 && link==0 && rowss==NR) {print all"\n"$0} \
#									{if(chrr!=$1 && link==0 && rowss==NR) {print all"\n"$0} \
#									{if(chrr!=$1 && link==1 && rowss==NR) {print chrr"\t"chrq"\t"start"\t"end"\t"end-start"\ttransloc\t"startq"\t"endq"\n"$0}}}}}}}}}}}}' \
#									>""$prefix"_ref".fragments_more500concat
#
#cat ${prefix}_ref.fragments_more500concat | awk '{if($2 ~ /_/) print $0}' > ${prefix}_ref.fragments_more500concat2
#cp ${prefix}_ref.fragments_more500concat ${prefix}_ref.fragments_more500concat2
#cat ${prefix}_ref.fragments_more500concat2 | awk '{print $2}' | sed 's/_/\t/g' > ${prefix}_ref.fragments_more500concat2_list
#rm ${prefix}_ref.fragments_more500concat
#cat "$alignments_folder/""$prefix"_ref".delta_filter.coordsg" | awk '{print $14}' | sort | uniq > chromosome_list.txt
#paste ${prefix}_ref.fragments_more500concat2 ${prefix}_ref.fragments_more500concat2_list > ${prefix}_ref.fragments_more500concat3
#awk 'NR==FNR{a[$1]++;next}; a[$10] {print $0}' chromosome_list.txt ${prefix}_ref.fragments_more500concat3 > ${prefix}_ref.fragments_more500concat4
#cat ${prefix}_ref.fragments_more500concat4 | awk '{if($1==$10 || $1==$9) print $0}' | cut -d $'\t' -f1-8 > ${prefix}_ref.fragments_more500concat

#rm ""$prefix"_ref".coordsg_matched_subtelofilt2
#rm ""$prefix"_ref".fragments
#rm ""$prefix"_ref".fragments_less500
#rm ""$prefix"_ref".fragments_more500
#rm ""$prefix"_ref".fragments_big
#rm ""$prefix"_ref".fragments_more500concat2
#rm ""$prefix"_ref".fragments_more500concat2_list
#rm chromosome_list.txt
#rm ""$prefix"_ref".fragments_more500concat3
#rm ""$prefix"_ref".fragments_more500concat4

#############################################
##############LARGE INVERSIONS###############
#############################################

echo ""
echo "Checking alignment sense for inversions involving majority of single chromosome bases"
echo ""

##alignments with a change in the sense compared to the reference##
##if query contigs are not orientated to the reference prior to alignment, this will look like enitre contigs are inverted##
#cat ""$prefix"_ref".coordsg_matched | awk '{if($14==$15 && $12!=$13) print $14"\t"$15"\t"$1"\t"$2"\t"$5"\tinversion\t"$3"\t"$4}' > ""$prefix"_ref".coordsg_matched_largeinv_temp

##merge inverted fragments only split by 500bp##
#cat ""$prefix"_ref".coordsg_matched_largeinv_temp | \
#	awk 'BEGIN{chrr=""; startr=0; endr=0; startq=0; endq=0}; \
#				{if(chrr==$1 && $3-endr<500) {print $0"\tneighbour"; chrr=$1; startr=$3; endr=$4; startq=$7; endq=$8} \
#				else {print $0; chrr=$1; startr=$3; endr=$4; startq=$7; endq=$8}
#				{if(chrr!=$1) {print $0; chrr=$1; startr=$3; endr=$4; startq=$7; endq=$8}}}' \
#				> ""$prefix"_ref".coordsg_matched_largeinv_temp2


#rows1=$(wc -l ""$prefix"_ref".coordsg_matched_largeinv_temp2 | awk '{print $1}')
#cat ""$prefix"_ref".coordsg_matched_largeinv_temp2 | \
#	awk -v rowss="$rows1" 'BEGIN{agg="first"; chrr=""; chrq=""; startr=0; endr=0; startq=0; endq=0} \
#								{if(agg=="" && $9=="" && rowss!=NR) {print all; all=$0; startr=$3; endr=$4; startq=$7; endq=$8; agg=$9} \
#								{if($9=="neighbour" && $4<=endr && rowss!=NR){chrr=$1; chrq=$2; endq=$8; agg=$9} \
#								{if($9=="neighbour" && $4>=endr && rowss!=NR){chrr=$1; chrq=$2; endr=$4; endq=$8; agg=$9} \
#								{if(agg=="neighbour" && $9=="" && rowss!=NR) {print chrr"\t"chrq"\t"startr"\t"endr"\t"endr-startr"\tinversion\t"startq"\t"endq"\taggregated"; all=$0; startr=$3; endr=$4; startq=$7; endq=$8; agg=$9} \
#								{if(agg=="first"){all=$0; startr=$3; endr=$4; startq=$7; endq=$8; agg=$9} \
#								{if(rowss==NR && $9=="neighbour" && $4<=endr) {print $1"\t"$2"\t"startr"\t"endr"\t"endr-startr"\tinversion\t"startq"\t"endq"\taggregated"} \
#								{if(rowss==NR && $9=="neighbour" && $4>endr) {print $1"\t"$2"\t"startr"\t"$4"\t"$4-startr"\tinversion\t"startq"\t"endq"\taggregated"} \
#								{if(rowss==NR && $9=="" && agg=="neighbour") {print chrr"\t"chrq"\t"startr"\t"endr"\t"endr-startr"\tinversion\t"startq"\t"endq"\taggregated\n"$0} \
#								{if(rowss==NR && $9=="" && agg==""){print all"\n"$0}}}}}}}}}}'  > ""$prefix"_ref".coordsg_matched_largeinv_temp3

#cat ""$prefix"_ref".coordsg_matched_largeinv_temp3 | awk '{if($9=="aggregated") {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8} else {print $0}}' > ""$prefix"_ref".coordsg_matched_largeinv

#rm ""$prefix"_ref".coordsg_matched_largeinv_temp
#rm ""$prefix"_ref".coordsg_matched_largeinv_temp2
#rm ""$prefix"_ref".coordsg_matched_largeinv_temp3

echo ""
echo ""
echo "######################################################################################################"
echo "             USING NON-GLOBAL ALIGNMENT FOR INVERSIONS, DUPLICATIONS AND CONTRACTIONS"
echo "######################################################################################################"
echo ""
echo ""
echo "Matching chromosomes based on names, using '_' seperator and filtering chrMT"
echo ""

#tail -n +5 "$alignments_folder/""$prefix"_ref".delta_filter.coords" | sed 's/_/\t/' | awk '!/chrMT/{if($14==$15) {print $0}; if($14==$16) {print $0}}' |\
#						awk '{if($16 == "") {print $0}} \
#						{if($16!="") {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"_"$16}}' \
#						> ""$prefix"_ref".coords_matched

awk 'NR==FNR{a[$1,$2]++;next};a[$14,$15] {print $0}' ""$prefix"_chromosome_pairs.txt" "$alignments_folder/""$prefix"_ref".delta_filter.coords" > ""$prefix"_ref".coords_matched


##All those calls that start before the end of the subtelomeric region or finish after are filtered
if [[ $filter_subtelomeric_region == "yes" ]]
then
	awk 'FNR==NR{a[$1]=$2"\t"$3;next} ($14 in a) {print $0"\t"a[$14]}' $subtelo_coords ""$prefix"_ref".coords_matched  |\
		awk '{if($2 > $16 && $1 < $17) {print $0}}' |\
			awk -F'\t' 'BEGIN { OFS = FS }; NF{NF-=2};1' > ""$prefix"_ref".coords_matched_subtelofilt
else
	cp ""$prefix"_ref".coords_matched ""$prefix"_ref".coords_matched_subtelofilt
fi



awk 'NR==FNR{a[$1,$2]++;next};a[$15,$14] {print $0}' ""$prefix"_chromosome_pairs.txt" "$alignments_folder/""$prefix"_query".delta_filter.coords" > ""$prefix"_query".coords_matched

##All those calls that start before the end of the subtelomeric region or finish after are filtered
if [[ $filter_subtelomeric_region == "yes" ]]
then
	awk 'FNR==NR{a[$1]=$2"\t"$3;next} ($14 in a) {print $0"\t"a[$14]}' $subtelo_coords ""$prefix"_query".coords_matched  |\
		awk '{if($2 > $16 && $1 < $17) {print $0}}' |\
			awk -F'\t' 'BEGIN { OFS = FS }; NF{NF-=2};1' > ""$prefix"_query".coords_matched_subtelofilt
else
	cp ""$prefix"_query".coords_matched ""$prefix"_query".coords_matched_subtelofilt
fi

#############################################
#################INVERSIONS##################
#############################################

echo ""
echo "Looking for interchromosomal changes in alignment sense, filtering for 1kb fragments and merging closely neighbouring calls"
echo ""

##Now looking for region within the same chromosome which swap sense and then take the size of the swapped sense
##Get the sense of the majority of the alignments between any two chromosomes and record as the same (1) or opposite (-1) as done with the alignment##
cat ""$prefix"_ref".coordsg_matched_subtelofilt | \
	awk '{print $0"\t"$13*$6}' |\
	awk '{a[$15]+=$16} END{for(i in a) {print i,a[i]}}' |\
	sed 's/ /\t/g' | while read chr
	do
		chr2=$( echo $chr | awk '{print $1}'  )
		sense=$( echo $chr | awk '{print $2}'  )
		cat ${prefix}_ref.coordsg_matched | awk '{print $14"\t"$15}' | sort -u  | awk -v chr2="$chr2" -v sense="$sense" '{if($2==chr2) print $1"\t"$2"\t"sense}'
	done |\
	awk '{if($3>0) {print $1"\t"$2"\t1"} else {print $1"\t"$2"\t-1"}}' > ""$prefix"_ref".coords_matched.sense

awk 'FNR==NR{c[$1, $2, $3]++;next} (($14, $15, 1) in c) {print $0"\t1"}' ""$prefix"_ref".coords_matched.sense ""$prefix"_ref".coords_matched_subtelofilt > ""$prefix"_ref".coords_matched_subtelofilt_sensetemp1
awk 'FNR==NR{c[$1, $2, $3]++;next} (($14, $15, -1) in c) {print $0"\t-1"}' ""$prefix"_ref".coords_matched.sense ""$prefix"_ref".coords_matched_subtelofilt > ""$prefix"_ref".coords_matched_subtelofilt_sensetemp2
cat ""$prefix"_ref".coords_matched_subtelofilt_sensetemp1 ""$prefix"_ref".coords_matched_subtelofilt_sensetemp2 | \
sort -k14,14 -k1V > ""$prefix"_ref".coords_matched_subtelofilt_sense

##look for change in global alignment sense and local##
cat ""$prefix"_ref".coords_matched_subtelofilt_sense | \
	awk '{if($13!=$16) {print $0"\tinversion"}}' > ""$prefix"_ref".coords_matched_subtelofilt_sense_inv
##trimmed and inverted coordinates in query##
cat ""$prefix"_ref".coords_matched_subtelofilt_sense_inv | \
	awk '{print $14"\t"$15"\t"$1"\t"$2"\t"$5"\tinversion\t"$4"\t"$3}' > ""$prefix"_ref".coords_matched_subtelofilt_sense_inv_trimmed
cat ""$prefix"_ref".coords_matched_subtelofilt_sense_inv_trimmed| \
	awk '$5>=50{print $0}' > ""$prefix"_ref".inv_50
cat ""$prefix"_ref".inv_50 | \
	awk '$5>=1000{print $0}' > ""$prefix"_ref".inv_1000

###label any overlapping inversion calls or within 500bp (can be changed but seems fair for 1kb calls) as neighbours###
cat ""$prefix"_ref".inv_1000 | \
	awk 'BEGIN{chrr=""; startr=0; endr=0; startq=0; endq=0}; \
				{if(chrr==$1 && $3-endr<500) {print $0"\tneighbour"; chrr=$1; startr=$3; endr=$4; startq=$7; endq=$8} \
				else {print $0; chrr=$1; startr=$3; endr=$4; startq=$7; endq=$8}
				{if(chrr!=$1) {print $0; chrr=$1; startr=$3; endr=$4; startq=$7; endq=$8}}}' \
				> ""$prefix"_ref".inv_1000_neighbours

##this is just used to help print the last line during aggregation
rows2=$(wc -l ""$prefix"_ref".inv_1000_neighbours | awk '{print $1}')
###aggregate all neighbours and take the proper start and end in reference, needs to have same name in both ref and query chromosome###
cat ""$prefix"_ref".inv_1000_neighbours | \
	awk -v rowss="$rows2" 'BEGIN{agg="first"; chrr=""; chrq=""; startr=0; endr=0; startq=0; endq=0} \
								{if(agg=="" && $9=="" && rowss!=NR) {print all; all=$0; startr=$3; endr=$4; startq=$7; endq=$8; agg=$9} \
								{if($9=="neighbour" && $4<=endr && rowss!=NR){chrr=$1; chrq=$2; endq=$8; agg=$9} \
								{if($9=="neighbour" && $4>=endr && rowss!=NR){chrr=$1; chrq=$2; endr=$4; endq=$8; agg=$9} \
								{if(agg=="neighbour" && $9=="" && rowss!=NR) {print chrr"\t"chrq"\t"startr"\t"endr"\t"endr-startr"\tinversion\t"startq"\t"endq"\taggregated"; all=$0; startr=$3; endr=$4; startq=$7; endq=$8; agg=$9} \
								{if(agg=="first"){all=$0; startr=$3; endr=$4; startq=$7; endq=$8; agg=$9} \
								{if(rowss==NR && $9=="neighbour" && $4<=endr) {print $1"\t"$2"\t"startr"\t"endr"\t"endr-startr"\tinversion\t"startq"\t"endq"\taggregated"} \
								{if(rowss==NR && $9=="neighbour" && $4>endr) {print $1"\t"$2"\t"startr"\t"$4"\t"$4-startr"\tinversion\t"startq"\t"endq"\taggregated"} \
								{if(rowss==NR && $9=="" && agg=="neighbour") {print chrr"\t"chrq"\t"startr"\t"endr"\t"endr-startr"\tinversion\t"startq"\t"endq"\taggregated\n"$0} \
								{if(rowss==NR && $9=="" && agg==""){print all"\n"$0}}}}}}}}}}' \
								> ""$prefix"_ref".inv_1000_neighbours_agg

rm ""$prefix"_ref".coords_matched.sense
rm ""$prefix"_ref".coords_matched_subtelofilt_sensetemp1
rm ""$prefix"_ref".coords_matched_subtelofilt_sensetemp2
rm ""$prefix"_ref".coords_matched_subtelofilt_sense_inv
rm ""$prefix"_ref".coords_matched_subtelofilt_sense_inv_trimmed
rm ""$prefix"_ref".inv_1000_neighbours

#############################################
################DUPLICATIONS#################
#############################################

echo ""
echo "Locating alignment overlaps for duplication assignment"
echo ""

##Get all alignments that overlap greater than 50bp with the previous alignment, tag those were sense gets inverted## 
##change value of query overlap if both strands in opposite sense to ref## 

cat ""$prefix"_ref".coords_matched_subtelofilt |\
	awk 'BEGIN{chrr=""; chrq=""; endr=0; endq=0} \
		{if(chrr==$14 && chrq==$15 && $13==sense && endr>=$1+50) {print $0"\t"startr"\t"endr"\t"startq"\t"endq"\t"$1-endr"\t"sense"\t"$3-endq"\ttandem_dup"; startr=$1; endr=$2; startq=$3; endq=$4; sense=$13} \
		{if(chrr==$14 && chrq==$15 && $13!=sense && endr>=$1+50) {print $0"\t"startr"\t"endr"\t"startq"\t"endq"\t"$1-endr"\t"sense"\tINVERTED\ttandem_dup_rev"; startr=$1; endr=$2; startq=$3; endq=$4; sense=$13} \
		{if(chrr==$14 && chrq==$15 && endr<$1+50) {startr=$1; endr=$2; startq=$3; endq=$4; sense=$13} \
		{if(chrr==$14 && chrq!=$15) {chrq=$15; startr=$1; endr=$2; startq=$3; endq=$4; sense=$13} \
		{if(chrr!=$14) {chrr=$14; chrq=$15; startr=$1; endr=$2; startq=$3; endq=$4; sense=$13}}}}}}' |\
			awk '{if($22!="INVERTED" && $13==-1 && $21==-1) \
				{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22*-1"\t"$23}\
				else {print $0}}' > ""$prefix"_ref".coords_matched_subtelofilt2

##remove rDNA##	
if [[ $rDNA_filter == "yes" ]]
then
	cat ""$prefix"_ref".coords_matched_subtelofilt2 |\
		awk '{if($14=="chrXII" && $1<=445500) {print$0} \
			{if($14=="chrXII" && $2>=485500) {print $0}\
			{if($14!="chrXII") {print $0}}}}' > ""$prefix"_ref".coords_overlaps
else
	cp ""$prefix"_ref".coords_matched_subtelofilt2 ""$prefix"_ref".coords_overlaps
fi

##list of pure clean, less than 50bp overlaps or gaps in the query and with an overlap in the reference##
cat ""$prefix"_ref".coords_overlaps | awk '{if($22!="INVERTED" && $22<=50 && $22>=-50) {print $0"\tclean"}}' > ""$prefix"_ref".coords_overlaps_clean

##list of positive gaps with and without corresponding deletion##
cat ""$prefix"_ref".coords_overlaps | awk '{if($22!="INVERTED" && $22>50) {print $0}}' > ""$prefix"_ref".coords_overlaps_pos

##Only print those overlaps, with a gap following the overlap, with those that match locations of insertions in the 50bp calls in both senses##
awk 'FNR==NR{c[$2, $3, $4]++;next} (($15, $17, $1) in c) > 0' ""$prefix"_query".gaps_50 ""$prefix"_ref".coords_overlaps_pos > ""$prefix"_ref".coords_overlaps_pos_tempforward
awk 'FNR==NR{c[$2, $3, $4]++;next} (($15, $1, $17) in c) > 0' ""$prefix"_query".gaps_50 ""$prefix"_ref".coords_overlaps_pos > ""$prefix"_ref".coords_overlaps_pos_tempreverse
cat ""$prefix"_ref".coords_overlaps_pos_tempforward ""$prefix"_ref".coords_overlaps_pos_tempreverse > ""$prefix"_ref".coords_overlaps_pos_match

cat ""$prefix"_ref".coords_overlaps_pos_match ""$prefix"_ref".coords_overlaps_clean |\
	sort -k1,1 -k3V |\
		awk '{print $14"\t"$15"\t"$1"\t"$17"\t"$20*-1"\tduplication\t"$19"\t"$3}' > ""$prefix"_ref".coords_dups_50

cat ""$prefix"_ref".coords_dups_50 | awk '$5>=1000{print $0}' > ""$prefix"_ref".coords_dups_1000

rm ""$prefix"_ref".coords_matched_subtelofilt2
rm ""$prefix"_ref".coords_overlaps
rm ""$prefix"_ref".coords_overlaps_clean
rm ""$prefix"_ref".coords_overlaps_pos
rm ""$prefix"_ref".coords_overlaps_pos_tempforward
rm ""$prefix"_ref".coords_overlaps_pos_tempreverse
rm ""$prefix"_ref".coords_overlaps_pos_match

#############################################
################CONTRACTIONS#################
#############################################

echo ""
echo "Locating alignment overlaps for contraction assignment"
echo ""

##Get all alignments that overlap greater than 50bp with the previous alignment, tag those were sense gets inverted## 
##change value of query overlap if both strands in opposite sense to ref## 

cat ""$prefix"_query".coords_matched_subtelofilt |\
	awk 'BEGIN{chrr=""; chrq=""; endr=0; endq=0} \
		{if(chrr==$14 && chrq==$15 && $13==sense && endr>=$1+50) {print $0"\t"startr"\t"endr"\t"startq"\t"endq"\t"$1-endr"\t"sense"\t"$3-endq"\ttandem_dup"; startr=$1; endr=$2; startq=$3; endq=$4; sense=$13} \
		{if(chrr==$14 && chrq==$15 && $13!=sense && endr>=$1+50) {print $0"\t"startr"\t"endr"\t"startq"\t"endq"\t"$1-endr"\t"sense"\tINVERTED\ttandem_dup_rev"; startr=$1; endr=$2; startq=$3; endq=$4; sense=$13} \
		{if(chrr==$14 && chrq==$15 && endr<$1+50) {startr=$1; endr=$2; startq=$3; endq=$4; sense=$13} \
		{if(chrr==$14 && chrq!=$15) {chrq=$15; startr=$1; endr=$2; startq=$3; endq=$4; sense=$13} \
		{if(chrr!=$14) {chrr=$14; chrq=$15; startr=$1; endr=$2; startq=$3; endq=$4; sense=$13}}}}}}' |\
			awk '{if($22!="INVERTED" && $13==-1 && $21==-1) \
				{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18"\t"$19"\t"$20"\t"$21"\t"$22*-1"\t"$23}\
				else {print $0}}' > ""$prefix"_query".coords_matched_subtelofilt2

##remove rDNA##	
if [[ $rDNA_filter == "yes" ]]
then
	cat ""$prefix"_query".coords_matched_subtelofilt2 |\
		awk '{if($14=="chrXII" && $1<=445500) {print$0} \
			{if($14=="chrXII" && $2>=485500) {print $0}\
			{if($14!="chrXII") {print $0}}}}' > ""$prefix"_query".coords_overlaps
else
	cp ""$prefix"_query".coords_matched_subtelofilt2 ""$prefix"_query".coords_overlaps
fi

##list of pure clean, less than 50bp overlaps or gaps in the query and with an overlap in the reference##
cat ""$prefix"_query".coords_overlaps | awk '{if($22!="INVERTED" && $22<=50 && $22>=-50) {print $0"\tclean"}}' > ""$prefix"_query".coords_overlaps_clean

##list of positive gaps with and without corresponding deletion##
cat ""$prefix"_query".coords_overlaps | awk '{if($22!="INVERTED" && $22>50) {print $0}}' > ""$prefix"_query".coords_overlaps_pos

##Only print those overlaps, with a gap following the overlap, with those that match locations of insertions in the 50bp calls in both senses##
awk 'FNR==NR{c[$2, $3, $4]++;next} (($15, $17, $1) in c) > 0' ""$prefix"_ref".gaps_50 ""$prefix"_query".coords_overlaps_pos > ""$prefix"_query".coords_overlaps_pos_tempforward
awk 'FNR==NR{c[$2, $3, $4]++;next} (($15, $1, $17) in c) > 0' ""$prefix"_ref".gaps_50 ""$prefix"_query".coords_overlaps_pos > ""$prefix"_query".coords_overlaps_pos_tempreverse
cat ""$prefix"_query".coords_overlaps_pos_tempforward ""$prefix"_query".coords_overlaps_pos_tempreverse > ""$prefix"_query".coords_overlaps_pos_match

cat ""$prefix"_query".coords_overlaps_pos_match ""$prefix"_query".coords_overlaps_clean |\
	sort -k1,1 -k3V |\
		awk '{print $15"\t"$14"\t"$19"\t"$3"\t"$20*-1"\tcontraction\t"$1"\t"$17}' > ""$prefix"_query".coords_dups_50

cat ""$prefix"_query".coords_dups_50 | awk '$5>=1000{print $0}' > ""$prefix"_query".coords_dups_1000

rm ""$prefix"_query".coords_matched_subtelofilt2
rm ""$prefix"_query".coords_overlaps
rm ""$prefix"_query".coords_overlaps_clean
rm ""$prefix"_query".coords_overlaps_pos
rm ""$prefix"_query".coords_overlaps_pos_tempforward
rm ""$prefix"_query".coords_overlaps_pos_tempreverse
rm ""$prefix"_query".coords_overlaps_pos_match

##############################################
########COMBINING CALLS FOR FILTERING#########
##############################################

echo ""
echo ""
echo "##############################################################################################################"
echo "Combining both deletions, insertion, inversions and duplications and identifying regions of more than one call"
echo "##############################################################################################################"
echo ""
echo ""
echo "Filtering clean inversions using deletion and insertion information"
echo ""

####Inversions will take precedance if previously both an insertion and deletion was labelled at the same location####
####This is due to the global alignment filter removing them  and therefore they are missing only due to this fault####
##only using 1kb for now, joining DELs, INSs, and INVs##
#cat ""$prefix"_ref".gaps_1000 ""$prefix"_query".gaps_1000 ""$prefix"_ref".inv_1000_neighbours_agg ""$prefix"_ref".fragments_more500concat | sort -k1,1 -k3V  > $prefix.1000bp_invfilt
cat ""$prefix"_ref".gaps_1000 ""$prefix"_query".gaps_1000 ""$prefix"_ref".inv_1000_neighbours_agg $prefix.transloc_candidate_alignments7_fragments10000 | sort -k1,1 -k3V  > $prefix.1000bp_invfilt

##filter out rDNA section##
if [[ $rDNA_filter == "yes" ]]
then
	cat $prefix.1000bp_invfilt | awk '{if($1=="chrXII" && $3<=445500) {print$0} \
										{if($1=="chrXII" && $4>=485500) {print $0}\
										{if($1!="chrXII") {print $0}}}}' \
											> $prefix.1000bp_invfilt_r
else
	cp $prefix.1000bp_invfilt $prefix.1000bp_invfilt_r
fi

##index all matching SV calls and label them with numbers, so each group has a unique value, then filter only those with assigned values##
##then in groups with inversion, deletion and insertion, take inversion as the call and deletion as the position##
rows3=$(wc -l $prefix.1000bp_invfilt_r | awk '{print $1}')

cat $prefix.1000bp_invfilt_r | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' | \
							awk -v rowss="$rows3" 'BEGIN{all=""; num=0; chrr=""; chrq=""; startr=0; endr=0; cont=0} \
								{if(rowss!=NR && chrr==$1 && chrq==$2 && startr>=$3-500 && startr<=$3+500 || rowss!=NR && chrr==$1 && chrq==$2 && endr>=$4-500 && endr<=$4+500 || rowss!=NR && chrr==$1 && chrq==$2 && startr<=$3 && endr>=$4 || rowss!=NR && chrr==$1 && chrq==$2 && startr>=$3 && endr<=$4) \
								{print all"\t"num; startr=$3; endr=$4; all=$0; cont=1} \
								{if(rowss!=NR && chrr==$1 && chrq==$2 && cont==0 && startr<$3-500 && endr<$4-500 || rowss!=NR && chrr==$1 && chrq==$2 && cont==0 && startr>=$3+500 && endr<=$4+500) \
								{print all; startr=$3; endr=$4; all=$0} \
								{if(rowss!=NR && chrr==$1 && chrq==$2 && cont==1 && startr<$3-500 && endr<$4-500 || rowss!=NR && chrr==$1 && chrq==$2 && cont==1 && startr>=$3+500 && endr<=$4+500) \
								{print all"\t"num; startr=$3; endr=$4; all=$0; num=num+1; cont=0} \
								{if(chrr=="") {chrr=$1; chrq=$2; startr=$3; endr=$4; all=$0; num=1} \
								{if(rowss!=NR && chrr==$1 && chrq!=$2 && cont==0) \
								{print all; chrr=$1; chrq=$2; startr=$3; endr=$4; all=$0} \
								{if(rowss!=NR && chrr==$1 && chrq!=$2 && cont==1) \
								{print all"\t"num; chrr=$1; chrq=$2; startr=$3; endr=$4; all=$0; num=num+1; cont=0} \
								{if(rowss!=NR && chrr!=$1 && cont==0) \
								{print all; chrr=$1; chrq=$2; startr=$3; endr=$4; all=$0} \
								{if(rowss!=NR && chrr!=$1 && cont==1) \
								{print all"\t"num; chrr=$1; chrq=$2; startr=$3; endr=$4; all=$0; num=num+1} \
								{if(rowss==NR && chrr==$1 && chrq==$2 && startr>=$3-500 && startr<=$3+500 || rowss==NR && chrr==$1 && chrq==$2 && endr>=$4-500 && endr<=$4+500 || rowss==NR && chrr==$1 && chrq==$2 && startr<=$3 && endr>=$4 ||rowss==NR && chrr==$1 && chrq==$2 && startr>=$3 && endr<=$4) \
								{print all"\t"num"\n"$0"\t"num} \
								{if(rowss==NR && chrr==$1 && chrq==$2 && cont==0 && startr<$3-500 && endr<$4-500 || rowss==NR && chrr==$1 && chrq==$2 && cont==0 && startr>=$3+500 && endr<=$4+500) \
								{print all"\n"$0} \
								{if(rowss==NR && chrr==$1 && chrq==$2 && cont==1 && startr<$3-500 && endr<$4-500 || rowss==NR && chrr==$1 && chrq==$2 && cont==1 && startr>=$3+500 && endr<=$4+500) \
								{print all"\t"num"\n"$0} \
								{if(rowss==NR && chrr!=$1 && cont==0) \
								{print all"\n"$0} \
								{if(rowss==NR && chrr!=$1 && cont==1) \
								{print all"\t"num"\n"$0}}}}}}}}}}}}}}' |\
									awk '$9!="" {print $0}' |\
									sort -k9,9V -k6 > $prefix.1000bp_invfilt_r_groups

##remove any value that only had it once, generally for those that have been linked to translocations##	 							
awk -F'\t' 'NR==FNR { dup[$0]; next; } $9 in dup' <(awk -F'\t' '{print $9}' $prefix.1000bp_invfilt_r_groups | sort | uniq -d) $prefix.1000bp_invfilt_r_groups > $prefix.1000bp_invfilt_r_groups2

##if inversion has insertion and deletion alongside it attach the deletion values##
##or if inversion is greater than 30kb use it anyway, it will be translocation associated
cat $prefix.1000bp_invfilt_r_groups2 | awk 'BEGIN{del=0; ins=0; inv=0; dex=1; sv=""; count=0} \
																			{if($9==dex && $6=="deletion" && count==0) {print $0; count=1; alldel=$0} \
																			{if($9==dex && $6=="insertion" && count==1) {print $0; count=2} \
																			{if($9==dex && $6=="insertion" && count==0) {print $0} \
																			{if($9==dex && $6=="inversion" && count==0 && $5 < 30000) {print $0; count=0} \
																			{if($9==dex && $6=="inversion" && count==1 && $5 < 30000) {print $0; count=0} \
																			{if($9==dex && $6=="inversion" && count==0 && $5 >= 30000) {print $0"\t"$0; count=0} \
																			{if($9==dex && $6=="inversion" && count==1 && $5 >= 30000) {print $0"\t"$0; count=0} \
																			{if($9==dex && $6=="inversion" && count==2) {print $0"\t"alldel; count=0} \
																			{if($9!=dex && $6=="deletion" && count==0) {print $0; count=1; alldel=$0; dex=$9} \
																			{if($9!=dex && $6=="deletion" && count==1) {print $0; count=1; alldel=$0; dex=$9} \
																			{if($9!=dex && $6=="deletion" && count==2) {print $0; count=1; alldel=$0; dex=$9} \
																			{if($9!=dex && $6=="inversion" && $5 >= 30000) {print $0"\t"$0; count=1; alldel=$0; dex=$9}}}}}}}}}}}}}'\
																			> $prefix.1000bp_invfilt_r_groups3

##only take those INVs that have the deletion values alongside it##
maxdex=$(sort -k9 -n $prefix.1000bp_invfilt_r_groups3 | tail -1 | awk '{print $9}')
tac $prefix.1000bp_invfilt_r_groups3 | awk -v maxdexx="$maxdex" 'BEGIN{inv=0; dex=maxdexx} \
														{if(dex==$9 && $10=="" && inv==0) {} \
														{if(dex==$9 && $10=="" && inv==1) {} \
														{if(dex==$9 && $10!=""){print $1"\t"$2"\t"$12"\t"$13"\t"$14"\t"$6"\t"$16"\t"$17"\t"$18; inv=1}
														{if(dex!=$9 && $10=="") {dex=$9; inv=0} \
														{if(dex!=$9 && $10!=""){print $1"\t"$2"\t"$12"\t"$13"\t"$14"\t"$6"\t"$16"\t"$17"\t"$18; inv=1; dex=$9}}}}}}' |\
														sort -k9,9V -k6 |\
														awk 'BEGIN{FS=OFS="\t"} NF{--NF};1' > $prefix.1000bp_inversionsclean

rm $prefix.1000bp_invfilt
rm $prefix.1000bp_invfilt_r
rm $prefix.1000bp_invfilt_r_groups
rm $prefix.1000bp_invfilt_r_groups2
rm $prefix.1000bp_invfilt_r_groups3

#######################################################################
########COMBINING CALLS AND REMOVING ALIGNMENT SPECIFIC ISSUES#########
#######################################################################

echo ""
echo "Combining all called SVs"
echo ""


#####50bp-1kb######
##remove insertions and deletion due to mummer3 alignment gaps issue which occur exactly the same in both ref and query yet don't correspond to an inversion##
##therefore remove any deletion or insertion calls with exactly the same coords in ref and query##
tac ""$prefix"_ref".gaps_50 ""$prefix"_query".gaps_50 |\
	sort -k1,1 -k3V |\
		awk 'BEGIN{startr=""; endr=""; startq=""; endq=""} {if(startr==$3 && endr==$4 && startq==$7 && endq==$8 || startr==$3 && endr==$4 && startq==$8 && endq==$7 || startr==$4 && endr==$3 && startq==$7 && endq==$8 || startr==$4 && endr==$3 && startq==$8 && endq==$7) {print $0"\texact\t"size} else {print $0; startr=$3; endr=$4; startq=$7; endq=$8; size=$5}}' >  $prefix.gaps_50filter_temp


tac $prefix.gaps_50filter_temp |\
	awk 'BEGIN{startr=""; endr=""; startq=""; endq=""} {if(startr==$3 && endr==$4 && startq==$7 && endq==$8 || startr==$3 && endr==$4 && startq==$8 && endq==$7 || startr==$4 && endr==$3 && startq==$7 && endq==$8 || startr==$4 && endr==$3 && startq==$8 && endq==$7) {print $0"\texact\t"size} else {print $0; startr=$3; endr=$4; startq=$7; endq=$8; size=$5}}' |\
		awk '{if($9=="exact" && $5 > $10) {print $0"\t"$5-$10} else if($9=="exact" && $10 > $5) {print $0"\t"$10-$5} else {print $0}}' |\
			awk '{if($11>5000) {print $0} else if($9!="exact"){print $0}}' |\
				awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' |\
					sort -k1,1 -k3V > $prefix.gaps_50filtered


##50bp calls do not have inversions YET, need to filter them using INDEL info as before MAYBE....##
cat $prefix.gaps_50filtered ${prefix}_ref.coords_dups_50 ${prefix}_query.coords_dups_50 | sort -V -k1 -k3 | sort -V -k1  |\
	awk '{if($5<1000) print $0}'> $prefix.50bp

##filter out rDNA section##
if [[ $rDNA_filter == "yes" ]]
then
	cat $prefix.50bp | awk '{if($1=="chrXII" && $3<=445500) {print$0} \
							{if($1=="chrXII" && $4>=485500) {print $0}\
							{if($1!="chrXII") {print $0}}}}' \
								>$prefix.50bp_r
else
	cp $prefix.50bp $prefix.50bp_r
fi

##removes doubles, more important for 50bp as it appears more, possibly marker for inversions....Could do as is done with 1kb inversions##
cat $prefix.50bp_r | awk 'BEGIN{start=0; end=0} {if(start==$3 || end==$4) {print $0"\tdouble" ; start=$3; end=$4} \
														{if(start==0 && end ==0 || start!=$3 || end!=$4) {print $0; start=$3; end=$4}}}' > $prefix.50bp_doubles

#####1kb->######
##remove insertions and deletions due to mummer3 alignment gaps issue which occur exactly the same in both ref and query yet don't correspond to an inversion##
##therefore remove any deletion or insertion calls with exactly the same coords in ref and query##
tac ""$prefix"_ref".gaps_1000 ""$prefix"_query".gaps_1000 |\
	sort -k1,1 -k3V |\
		awk 'BEGIN{startr=""; endr=""; startq=""; endq=""} {if(startr==$3 && endr==$4 && startq==$7 && endq==$8 || startr==$3 && endr==$4 && startq==$8 && endq==$7 || startr==$4 && endr==$3 && startq==$7 && endq==$8 || startr==$4 && endr==$3 && startq==$8 && endq==$7) {print $0"\texact\t"size} else {print $0; startr=$3; endr=$4; startq=$7; endq=$8; size=$5}}' >  $prefix.gaps_1000filter_temp

tac $prefix.gaps_1000filter_temp |\
	awk 'BEGIN{startr=""; endr=""; startq=""; endq=""} {if(startr==$3 && endr==$4 && startq==$7 && endq==$8 || startr==$3 && endr==$4 && startq==$8 && endq==$7 || startr==$4 && endr==$3 && startq==$7 && endq==$8 || startr==$4 && endr==$3 && startq==$8 && endq==$7) {print $0"\texact\t"size} else {print $0; startr=$3; endr=$4; startq=$7; endq=$8; size=$5}}' |\
		awk '{if($9=="exact" && $5 > $10) {print $0"\t"$5-$10} else if($9=="exact" && $10 > $5) {print $0"\t"$10-$5} else {print $0}}' |\
			awk '{if($11>5000) {print $0} else if($9!="exact"){print $0}}' |\
				awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' |\
					sort -k1,1 -k3V > $prefix.gaps_1000filtered

cat $prefix.gaps_1000filtered ""$prefix"_ref".coords_dups_1000 ""$prefix"_query".coords_dups_1000 $prefix.transloc_candidate_alignments7_fragments10000 $prefix.1000bp_inversionsclean | sort -k1,1 -k3V  > $prefix.1000bp

##filter out rDNA section##
if [[ $rDNA_filter == "yes" ]]
then
	cat $prefix.1000bp | awk '{if($1=="chrXII" && $3<=445500) {print$0} \
							{if($1=="chrXII" && $4>=485500) {print $0}\
							{if($1!="chrXII") {print $0}}}}' \
								>$prefix.1000bp_r
else
	cp $prefix.1000bp $prefix.1000bp_r
fi

##labelling those that overlap or are very close to the other, seperately defining if an inversion is involved or not##
##then perform it in the opposite sense and once more##

tac $prefix.1000bp_r | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' | \
						awk 'BEGIN{start=0; end=0; sv=""} \
								{if(start>=$3-500 && start<=$3+500 && sv=="inversion" || end>=$4-500 && end<=$4+500 && sv=="inversion" || start<=$3 && end>=$4 && sv=="inversion" || start>=$3 && end<=$4 && sv=="inversion") \
								{print $0"\tdouble_inversion"; start=$3; end=$4} \
								{if(start>=$3-500 && start<=$3+500 && $6=="inversion" || end>=$4-500 && end<=$4+500 && $6=="inversion" || start<=$3 && end>=$4 && $6=="inversion" || start>=$3 && end<=$4 && $6=="inversion") \
								{print $0"\tdouble_inversion"; start=$3; end=$4; sv=$6} \
								{if(start>=$3-500 && start<=$3+500 && sv!="inversion" || end>=$4-500 && end<=$4+500 && sv!="inversion" || start<=$3 && end>=$4 && sv!="inversion" || start>=$3 && end<=$4 && sv!="inversion") \
								{print $0"\tdouble"; start=$3; end=$4; sv=$6} \
								{if(start<=$3-500 || start>=$3+500 || end<=$4-500 || end>=$4+500) {print $0; start=$3; end=$4; sv=$6} \
								{if(start==0 && end==0) {print $0; start=$3; end=$4; sv=$6}}}}}}' \
								> $prefix.1000bp_r_doubles

tac $prefix.1000bp_r_doubles | \
						awk 'BEGIN{start=0; end=0; sv=""} \
								{if(start>=$3-500 && start<=$3+500 && sv=="inversion" || end>=$4-500 && end<=$4+500 && sv=="inversion" || start<=$3 && end>=$4 && sv=="inversion" || start>=$3 && end<=$4 && sv=="inversion") \
								{print $0"\tdouble_inversion"; start=$3; end=$4} \
								{if(start>=$3-500 && start<=$3+500 && $6=="inversion" || end>=$4-500 && end<=$4+500 && $6=="inversion" || start<=$3 && end>=$4 && $6=="inversion" || start>=$3 && end<=$4 && $6=="inversion") \
								{print $0"\tdouble_inversion"; start=$3; end=$4; sv=$6} \
								{if(start>=$3-500 && start<=$3+500 && sv!="inversion" || end>=$4-500 && end<=$4+500 && sv!="inversion" || start<=$3 && end>=$4 && sv!="inversion" || start>=$3 && end<=$4 && sv!="inversion") \
								{print $0"\tdouble"; start=$3; end=$4; sv=$6} \
								{if(start<=$3-500 || start>=$3+500 || end<=$4-500 || end>=$4+500) {print $0; start=$3; end=$4; sv=$6} \
								{if(start==0 && end==0) {print $0; start=$3; end=$4; sv=$6}}}}}}' \
								> $prefix.1000bp_r_doubles2


cat $prefix.1000bp_r_doubles2 | awk '{if($6!="transloc" && $9=="double" && $10=="double_inversion") {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$10} \
									{if($6!="transloc" && $9=="double_inversion" && $10=="double" || $6!="transloc" && $9=="double_inversion" && $10=="" || $9=="double_inversion" && $10=="double_inversion") {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9} \
									{if($6!="transloc" && $9=="double" && $10=="") {print $0}\
									{if($6=="transloc" || $9=="") {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}\
									{if($6!="transloc" && $9=="double" && $10=="double") {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}}}}}}' \
									> $prefix.1000bp_r_doubles3

rm $prefix.gaps_50filter_temp
rm $prefix.gaps_50filtered
rm $prefix.50bp
rm $prefix.50bp_r
rm $prefix.gaps_1000filter_temp
rm $prefix.gaps_1000filtered
rm $prefix.1000bp
rm $prefix.1000bp_r
rm $prefix.1000bp_r_doubles
rm $prefix.1000bp_r_doubles2

echo ""
echo "Removing deletions and insertions called during global alignment due to inversions"
echo ""

##label 'double' as complicated and then remove insertions and/or deletions from 'double_inversion' calls##
##take size and position in reference from a deletion primarily, secondly insertion if no deletion present in mixed calls##

cat $prefix.1000bp_r_doubles3 | awk 'BEGIN{size=0; sv=""; all=""; dup=""; startref=0; endref=0; startq=0; endq=0} \
									{if($9=="double") {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tcomplicated"; dup=""; size=0; sv=""} \
									{if($9=="" && dup=="") {print $0; dup=""; size=0; sv=""} \

									{if($9=="double_inversion" && sv=="" && $6=="insertion") {size=$5; sv=$6; startref=$3; endref=$4; startq=$7; endq=$8} \
									{if($9=="double_inversion" && sv=="" && $6=="deletion") {size=$5; sv=$6; startref=$3; endref=$4; startq=$7; endq=$8} \
									{if($9=="double_inversion" && sv=="" && $6=="inversion") {all=$0; dup=$9; sv=$6} \

									{if($9=="double_inversion" && sv=="deletion" && $6=="insertion") {}} \
									{if($9=="double_inversion" && sv=="deletion" && $6=="inversion") {print $0"\t"startref"\t"endref"\t"size"\t"startq"\t"endq; dup=""; sv=""; size=0} \

									{if($9=="double_inversion" && sv=="insertion" && $6=="inversion") {all=$0; dup=$9; sv=$6} \
									{if($9=="double_inversion" && sv=="insertion" && $6=="deletion") {size=$5; sv=$6; startref=$3; endref=$4; startq=$7; endq=$8} \

									{if($9=="double_inversion" && sv=="inversion" && $6=="deletion") {print all"\t"$3"\t"$4"\t"$5"\t"$7"\t"$8; sv=""; dup=""; all=""; size=0} \
									{if($9=="double_inversion" && sv=="inversion" && $6=="insertion") {size=$5} \

									{if($9=="" && dup=="double_inversion") {print all"\t"size"\n"$0; dup=""; size=0; sv=""} \
									{if($9=="double" && dup=="double_inversion") {print all"\t"size"\n"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\tcomplicated"; dup=""; size=0; sv=""}}}}}}}}}}}}}' |\
									
									awk '{if($10=="") {print $0} \
										{if($10!="") {print $1"\t"$2"\t"$10"\t"$11"\t"$12"\t"$6"\t"$13"\t"$14"\t"$9}}}' \
									> $prefix.1000bp_filteredcalls_temp

cat $prefix.1000bp_filteredcalls_temp |\
	awk '{if($9=="double_inversion") {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8} else {print $0}}' > $prefix.1000bp_filteringgaps1

cat $prefix.1000bp_filteringgaps1 | awk '{if($6=="insertion" || $6=="deletion") print $0}' > $prefix.gaps_1000filteredinv

#cat $prefix.gaps_1000filteredinv ""$prefix"_ref".coords_dups_1000 ""$prefix"_query".coords_dups_1000 ""$prefix"_ref".fragments_more500concat $prefix.1000bp_inversionsclean ""$prefix"_ref".coordsg_matched_largeinv  | sort -k1,1 -k3V  > $prefix.1000bp2
cat $prefix.gaps_1000filteredinv ""$prefix"_ref".coords_dups_1000 ""$prefix"_query".coords_dups_1000 $prefix.transloc_candidate_alignments7_fragments10000 $prefix.1000bp_inversionsclean  | sort -k1,1 -k3V  > $prefix.1000bp2


if [[ $rDNA_filter == "yes" ]]
then
	cat $prefix.1000bp2| awk '{if($1=="chrXII" && $3<=445500) {print$0} \
							{if($1=="chrXII" && $4>=485500) {print $0}\
							{if($1!="chrXII") {print $0}}}}' \
								>$prefix.1000bp_filteredcalls
else
	cp $prefix.1000bp2 $prefix.1000bp_filteredcalls
fi

########REMOVING INSERTIONS THAT ARE ALSO TANDEM DUPLICATIONS##########

##orientate all calls so that it moves from the beginning to the end of each chromosome
##except for translocations to maintain the correct borders with the correct orientation
cat $prefix.1000bp_filteredcalls $prefix.50bp_doubles | sort -k1,1 -k3V | awk '{if($6 != "transloc" && $3 > $4) {print $1"\t"$2"\t"$4"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8} \
						else {print $0}}' \
						| awk '{if($6 != "transloc" && $7 > $8) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$8"\t"$7} \
						else {print $0}}' > $prefix.all_temp

##select only duplication and insertions and filter insertions where they match duplication regions
cat $prefix.all_temp | awk '{if($6=="duplication" || $6=="insertion") print $0}' > $prefix.all_dup_ins
cat $prefix.all_dup_ins | awk  'BEGIN {chrq=""; startr=""; stopq=""} \
									{if($6=="duplication") {print $0; chrq=$2; startr=$3; endq=$8} \
									else {if($6=="insertion" && $2==chrq && $3>startr-50 && $3<startr+50 && $8>endq-50 && $8<endq+50) {} \
									else {print $0}}}' > $prefix.all_dup_ins_filtered

##combine new set with others
cat $prefix.all_temp | awk '{if($6!="duplication" && $6!="insertion") {print $0}}' > $prefix.all_else
cat $prefix.all_else $prefix.all_dup_ins_filtered | awk '!/mitochondrion/{print $0}' | sort -k1,1 -k3V > $prefix.filteredcalls

##select only contraction and deletions and filter insertions where they match contraction regions
cat $prefix.filteredcalls | awk '{if($6=="contraction" || $6=="deletion") print $0}' > $prefix.all_contr_del
cat $prefix.all_contr_del | awk  'BEGIN {chrq=""; startr=""; stopq=""} \
									{if($6=="contraction") {print $0; chrq=$2; startr=$3; endq=$8} \
									else {if($6=="deletion" && $2==chrq && $3>startr-50 && $3<startr+50 && $8>endq-50 && $8<endq+50) {} \
									else {print $0}}}' > $prefix.all_contr_del_filtered

##combine new set with others
cat $prefix.filteredcalls | awk '{if($6!="contraction" && $6!="deletion") {print $0}}' > $prefix.all_else
cat $prefix.all_else $prefix.all_contr_del_filtered | awk '!/mitochondrion/{print $0}' | sort -k1,1 -k3V > $prefix.filteredcalls


###REMVOING ALL INDELS BIGGER THAN 150kb###
cat $prefix.filteredcalls | awk '{if($6 == "deletion" && $5 > 150000 || $6 == "insertion" && $5 > 150000) {} else print $0}' | awk 'NF' > $prefix.filteredcalls2
rm $prefix.filteredcalls
mv $prefix.filteredcalls2 $prefix.filteredcalls


rm $prefix.1000bp_r_doubles3
rm $prefix.1000bp_filteredcalls_temp
rm $prefix.1000bp_filteringgaps1
rm $prefix.gaps_1000filteredinv
rm $prefix.1000bp2
rm $prefix.all_temp
rm $prefix.all_dup_ins
rm $prefix.all_dup_ins_filtered
rm $prefix.all_contr_del
rm $prefix.all_contr_del_filtered
rm $prefix.all_else

#################################
#####BLAST STEP IF SELECTED######
#################################

if [[ $blast_step == "yes" ]]
then
	
	echo ""
	echo "label INDELS events as novel or mobile elements"
	echo ""
	
	##reference assemblies##
	$SAMTOOLS faidx $reference_assembly
	$SAMTOOLS faidx $query_assembly
	$BLASTDB -in $reference_assembly -input_type fasta -parse_seqids -dbtype nucl
	touch $prefix.blast_in
	touch $prefix.insertion_blast
	touch $prefix.deletion_blast
	
	##GET FILE WITH JUST INSERTIONS OR DELETIONS##
	cat $prefix.filteredcalls | awk '{if($6~"insertion") print $0}' > $prefix.insertions
	cat $prefix.filteredcalls | awk '{if($6~"deletion") print $0}' > $prefix.deletions
	
	##read in insertion file, seperate by tab, and blast against blast db made with reference##
	##check to see if inserted fragment matches anything on reference##
	##if matches output same line as in original file and add 'mobile'##
	##add 'novel' if no matches found##

	while IFS=$'\t' read v1 v2 v3 v4 v5 v6 v7 v8 v9
	do
		ref_chr=$v1
		query_chr=$v2
		STARTr=$v3
		STOPr=$v4
		STARTq=$v7
		STOPq=$v8
		size=$v5
		alt=$v9
		echo "$query_chr	$STARTq	$STOPq" >> $prefix.blast_in
		$SAMTOOLS faidx $query_assembly $query_chr:$STARTq-$STOPq > query_segment.fa
		$BLASTN -db $reference_assembly -query query_segment.fa -max_target_seqs 1 -outfmt 6 >> $prefix.blast_in
		$BLASTN -db $reference_assembly -query query_segment.fa -max_target_seqs 1 -outfmt 6 > temp.blast1
		if [ -s temp.blast1 ]; then 
			echo "$ref_chr	$query_chr	$STARTr	$STOPr	$size	insertion_mobile	$STARTq	$STOPq	$alt" >> $prefix.insertion_blast
		else
			echo "$ref_chr	$query_chr	$STARTr	$STOPr	$size	insertion_novel	$STARTq	$STOPq	$alt" >> $prefix.insertion_blast
		fi
	
	done < $prefix.insertions


	$BLASTDB -in $query_assembly -input_type fasta -parse_seqids -dbtype nucl
	touch $prefix.blast_del
	
	while IFS=$'\t' read v1 v2 v3 v4 v5 v6 v7 v8 v9
	do
		ref_chr=$v1
		query_chr=$v2
		STARTr=$v3
		STOPr=$v4
		STARTq=$v7
		STOPq=$v8
		size=$v5
		alt=$v9
		echo "$ref_chr	$STARTr	$STOPr" >> $prefix.blast_del
		$SAMTOOLS faidx $reference_assembly $ref_chr:$STARTr-$STOPr > ref_segment.fa
		$BLASTN -db $query_assembly -query ref_segment.fa -max_target_seqs 1 -outfmt 6 >> $prefix.blast_del
		$BLASTN -db $query_assembly -query ref_segment.fa -max_target_seqs 1 -outfmt 6 > temp.blast2
		if [ -s temp.blast2 ]; then 
			echo "$ref_chr	$query_chr	$STARTr	$STOPr	$size	deletion_mobile	$STARTq	$STOPq	$alt" >> $prefix.deletion_blast
		else
			echo "$ref_chr	$query_chr	$STARTr	$STOPr	$size	deletion_novel	$STARTq	$STOPq	$alt" >> $prefix.deletion_blast
		fi
	
	done < $prefix.deletions

	##combine new INDELS labelled SVs with other calls again##

	cat $prefix.filteredcalls | awk '{if($6!="insertion" && $6!="deletion") print $0}' > $prefix.noINDEL
	echo "ref_chr	query_chr	ref_start	ref_stop	size	SV_type	query_start	query_stop" > $prefix.SVs_all.tsv
	cat $prefix.insertion_blast $prefix.deletion_blast $prefix.noINDEL | sort -k1,1 -k3V >> $prefix.SVs_all.tsv
	
	rm query_segment.fa
	rm ref_segment.fa
	rm temp.blast1
	rm temp.blast2
	rm $prefix.noINDEL

else
	cat $prefix.filteredcalls | sort -k1,1 -k3V  >> $prefix.SVs_all.tsv
fi



echo ""
echo "Extracting DNA involved in SV events"
echo ""

###getting DNA associated with each variant and adding it to output
##get insertion fragments from query assembly
cat ${prefix}.SVs_all.tsv | grep -v "ref_chr" | awk '{if($6 ~ "insertion" || $6 ~ "contraction") print $0}' | while read variant
do
	location=$( echo $variant | awk '{print $2":"$7"-"$8}' )
	fragment=$( samtools faidx ${query_assembly} $location | grep -v ">" | tr -d '\n')
	echo $variant" "$fragment | sed 's/ /\t/g' >> ${prefix}.SVs_all.withfragmentTEMP.tsv
done
# get DNA for all other variants from reference
cat ${prefix}.SVs_all.tsv | grep -v "ref_chr" | awk '{if($6 !~ "insertion" && $6!~"contraction") print $0}' | while read variant
do
	location=$( echo $variant | awk '{if($3 > $4) {print $1":"$4"-"$3} else {print $1":"$3"-"$4}}' )
	fragment=$( samtools faidx ${reference_assembly} $location | grep -v ">" | tr -d '\n')
	echo $variant" "$fragment | sed 's/ /\t/g' >> ${prefix}.SVs_all.withfragmentTEMP.tsv
done
cat ${prefix}.SVs_all.withfragmentTEMP.tsv | sort -k1,1 -k3V > ${prefix}.SVs_all.withfragmentTEMP2.tsv
cat ${prefix}.SVs_all.withfragmentTEMP2.tsv | awk '{if($9 == "complicated" || $9 == "double"){t=$9; $9=$10; $10=t; print} else {print $0}}' | sed 's/ /\t/g' > ${prefix}.SVs_all.withfragment.tsv
rm ${prefix}.SVs_all.withfragmentTEMP*


echo ""
echo "Adding translocation border information to last column"
echo ""

##translocations are currently only labelled as the fragment in involved. It can be used to figure out the borders manually per strain but doesn't allow easy strain to strain comparison
##the last column for translocations now contains the border information showing the edge/edges of the fragment involved in the event and with which other reference chromosome
##added to same column as notes for complicated or doubled calls
##follows similar conventions to vcf format
## edge on outside of hard bracket is where the current fragment is placed, i.e. it is here[*[   or    ]*]here[*[   or    ]*]here    whilst the ref chromosome and position for the adjacent fragment (indicated by the star) is placed within the two brackets


##adding to file without fragment added
##also reorientating the translocation fragments with the opposite sense to the ref after identifying the correct bordering fragment
cat ${prefix}.SVs_all.tsv | awk '{if($6 != "transloc") print}' > ${prefix}.SVs_all.notransloc
cat ${prefix}.SVs_all.tsv | awk '{if($6 == "transloc") print}'  | sort -k2,2 -k8n -k7n | awk '{print $2}' | sort -u | while read chr
do
	cat ${prefix}.SVs_all.tsv |\
		awk -v chr="$chr" '{if($2 == chr && $6 == "transloc" ) print}'  | sort -k2,2 -k8n -k7n | awk '{if(NR==1) {line=$0; prev=$1; pos=$4} else if(NR != "1" && mid != "T") {print line"\t]"$1":"$3"]"; mid="T"; line=$0; prev2=$1; pos2=$4 } else if(NR != "1" && mid == "T") {print line"\t["prev":"pos"[\n"line"\t]"$1":"$3"]" ; prev=prev2; pos=pos2 ;  prev2=$1; pos2=$4; line=$0 }} END{print line"\t["prev":"pos"[" }'
	done > ${prefix}.SVs_all.translocwithborder
echo "ref_chr	query_chr	ref_start	ref_stop	size	SV_type	query_start	query_stop	info	" > $prefix.SVs_all.tsv
cat ${prefix}.SVs_all.notransloc ${prefix}.SVs_all.translocwithborder | sort -k1,1 -k3V |\
awk '{if($3 > $4) {print $1"\t"$2"\t"$4"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9} else {print $0}}'>> ${prefix}.SVs_all.tsv

##adding to file with fragment
cat ${prefix}.SVs_all.withfragment.tsv | awk '{if($6 != "transloc") print}' > ${prefix}.SVs_all.withfragment.notransloc
cat ${prefix}.SVs_all.withfragment.tsv | awk '{if($6 == "transloc") print}'  | sort -k2,2 -k8n -k7n | awk '{print $2}' | sort -u | while read chr
do
	cat ${prefix}.SVs_all.withfragment.tsv |\
		awk -v chr="$chr" '{if($2 == chr && $6 == "transloc" ) print}'  | sort -k2,2 -k8n -k7n | awk '{if(NR==1) {line=$0; prev=$1; pos=$4} else if(NR != "1" && mid != "T") {print line"\t]"$1":"$3"]"; mid="T"; line=$0; prev2=$1; pos2=$4 } else if(NR != "1" && mid == "T") {print line"\t["prev":"pos"[\n"line"\t]"$1":"$3"]" ; prev=prev2; pos=pos2 ;  prev2=$1; pos2=$4; line=$0 }} END{print line"\t["prev":"pos"[" }'
	done > ${prefix}.SVs_all.withfragment.translocwithborder
echo "ref_chr	query_chr	ref_start	ref_stop	size	SV_type	query_start	query_stop	fragments	info" > $prefix.SVs_all.withfragment.tsv
cat ${prefix}.SVs_all.withfragment.notransloc ${prefix}.SVs_all.withfragment.translocwithborder | sort -k1,1 -k3V |\
awk '{if($3 > $4) {print $1"\t"$2"\t"$4"\t"$3"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10} else {print $0}}'>> ${prefix}.SVs_all.withfragment.tsv

rm ${prefix}.SVs_all.notransloc 
rm ${prefix}.SVs_all.translocwithborder
rm ${prefix}.SVs_all.withfragment.notransloc 
rm ${prefix}.SVs_all.withfragment.translocwithborder
 
echo ""
echo "Generating a VCF file"
echo ""

echo "##fileformat=VCFv4.3" >> ${prefix}.SVs_all.vcf
date +##fileDate=%Y%m%d >> ${prefix}.SVs_all.vcf
echo "##source=MUMandCo_v${version}" >> ${prefix}.SVs_all.vcf
echo "##reference=${reference_assembly}" >> ${prefix}.SVs_all.vcf
grep ">" ${reference_assembly} | sed 's/>//g' | while read chr
do
	length=$( cat ${reference_assembly} | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' | grep "\S" | grep -x -A 1 ">$chr" | grep -v '>' | awk '{print length}' )
	echo "##contig=<ID=${chr},length=${length}" >> ${prefix}.SVs_all.vcf
done
echo "##query=${query_assembly}" >> ${prefix}.SVs_all.vcf
grep ">" ${query_assembly} | sed 's/>//g' | while read chr
do
	length=$( cat ${query_assembly} | awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);} END {printf("\n");}' | grep "\S" | grep -x -A 1 ">$chr" | grep -v '>' | awk '{print length}' )
	echo "##query_contig=<ID=${chr},length=${length}" >> ${prefix}.SVs_all.vcf
done
echo "##INFO=<ID=END,Number=1,Type=Integer,Description="End position in the reference genome">" >> ${prefix}.SVs_all.vcf
echo "##INFO=<ID=SVLEN,Number=.,Type=Integer,Description="Difference in length between REF and ALT alleles">" >> ${prefix}.SVs_all.vcf
echo "##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">" >> ${prefix}.SVs_all.vcf
echo "##INFO=<ID=qCHR,Number=1,Type=String,Description="Chromosome in query genome">" >> ${prefix}.SVs_all.vcf
echo "##INFO=<ID=qSTART,Number=1,Type=Integer,Description="Start position in query genome">" >> ${prefix}.SVs_all.vcf
echo "##INFO=<ID=qEND,Number=1,Type=Integer,Description="End position in query genome">" >> ${prefix}.SVs_all.vcf
echo "##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">" >> ${prefix}.SVs_all.vcf
echo "##ALT=<ID=DEL,Description="Deletion">" >> ${prefix}.SVs_all.vcf
echo "##ALT=<ID=DUP:TANDEM,Description="Tandem Duplication">" >> ${prefix}.SVs_all.vcf
echo "##ALT=<ID=CONTR,Description="Contraction">" >> ${prefix}.SVs_all.vcf
echo "##ALT=<ID=INS,Description="Insertion of novel sequence">" >> ${prefix}.SVs_all.vcf
echo "##ALT=<ID=INV,Description="Inversion">" >> ${prefix}.SVs_all.vcf
echo "##ALT=<ID=TRA,Description="Region involved in translocation alternative to a breakend position">" >> ${prefix}.SVs_all.vcf
echo "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	$prefix" >> ${prefix}.SVs_all.vcf
cat ${prefix}.SVs_all.withfragment.tsv | tail -n +2 | while read variant
do
	##split up the line into individual features to fit it into the VCF format
		echo $variant | awk '{if($6 == "deletion") {print $1"\t"$3"\t.\t"$9"\t<DEL>\t.\tPASS\tEND="$4";SVLEN=-"$5";SVTYPE=DEL;qCHR="$2";qSTART="$7";qEND="$8"\tGT\t1/1"}\
								if($6 == "insertion") {print $1"\t"$3"\t.\t<INS>\t"$9"\t.\tPASS\tEND="$4";SVLEN="$5";SVTYPE=INS;qCHR="$2";qSTART="$7";qEND="$8"\tGT\t1/1"}\
								if($6 == "duplication") {print $1"\t"$3"\t.\t"$9"\t<DUP>\t.\tPASS\tEND="$4";SVLEN="$5";SVTYPE=DUP;qCHR="$2";qSTART="$7";qEND="$8"\tGT\t1/1"}\
								if($6 == "contraction") {print $1"\t"$3"\t.\t"$9"\t<CONTR>\t.\tPASS\tEND="$4";SVLEN="$5";SVTYPE=CONTR;qCHR="$2";qSTART="$7";qEND="$8"\tGT\t1/1"}\
								if($6 == "inversion") {print $1"\t"$3"\t.\t"$9"\t<INV>\t.\tPASS\tEND="$4";SVLEN="$5";SVTYPE=INV;qCHR="$2";qSTART="$7";qEND="$8"\tGT\t1/1"}\
								if($6 == "transloc") {print $1"\t"$3"\t.\t"$9"\t"$10"\t.\tPASS\tEND="$4";SVLEN="$5";SVTYPE=TRA;qCHR="$2";qSTART="$7";qEND="$8"\tGT\t1/1"}}' >> ${prefix}.SVs_all.vcf
done



echo ""
echo "Counting number of detected SVs"
echo ""

total=$(cat $prefix.filteredcalls | awk '{print $0}' | wc -l)
total_del=$(cat $prefix.filteredcalls | awk '/deletion/{print $0}' | wc -l)
total_ins=$(cat $prefix.filteredcalls | awk '/insertion/{print $0}' | wc -l)
#total_INDELmean=$(cat $prefix.filteredcalls | awk 'BEGIN{sum=0} !/double/&&/deletion/||/insertion/{sum+=$5} END{print sum/NR}')
total_dup=$(cat $prefix.filteredcalls | awk '/duplication/{print $0}' | wc -l)
total_contra=$(cat $prefix.filteredcalls | awk '/contraction/{print $0}' | wc -l)
total_inv=$(cat $prefix.filteredcalls | awk '/inversion/{print $0}' | wc -l)
total_trans=$(cat $prefix.filteredcalls | awk '!/complicated/&&/transloc/{print $0}' | wc -l)

echo ""
echo $prefix " Total SVs  = "$total
echo $prefix " Deletions  = "$total_del
echo $prefix " Insertions  = "$total_ins
#echo $prefix " INDEL mean size = "$total_INDELmean
echo $prefix " Duplications  = "$total_dup
echo $prefix " Contractions  = "$total_contra
echo $prefix " Inversions  = "$total_inv
echo $prefix " Translocations  = "$total_trans

touch $prefix.summary.txt
cat > $prefix.summary.txt <<EOF
$prefix
Total_SVs	$total
Deletions	$total_del
Insertions	$total_ins
Duplications	$total_dup
Contractions	$total_contra
Inversions	$total_inv
Translocations	$total_trans

EOF

if [[ $cleanup == "yes" ]]
then
	rm ""$prefix"_ref".coordsg_matched
	rm ""$prefix"_ref".coordsg_matched_subtelofilt
	rm ""$prefix"_query".coordsg_matched
	rm ""$prefix"_query".coordsg_matched_subtelofilt
	rm ""$prefix"_ref".coordsg_matched_subtelofilt_gaps
	rm ""$prefix"_query".coordsg_matched_subtelofilt_gaps
	rm $prefix.transloc_candidate_alignments7_fragments10000
	rm ${prefix}_ref.fragments_less10000
	rm ${prefix}_chromosome_pairs.txt
	rm ${prefix}_pre_list.txt
	rm ""$prefix"_ref".gaps_50
	rm ""$prefix"_ref".gaps_1000
	rm ""$prefix"_query".gaps_50
	rm ""$prefix"_query".gaps_1000
	#rm ""$prefix"_ref".fragments_more500concat
	#rm ""$prefix"_ref".coordsg_matched_largeinv
	rm ""$prefix"_ref".coords_matched
	rm ""$prefix"_ref".coords_matched_subtelofilt
	rm ""$prefix"_ref".coords_matched_subtelofilt_sense
	rm ""$prefix"_query".coords_matched
	rm ""$prefix"_query".coords_matched_subtelofilt
	rm ""$prefix"_ref".inv_50
	rm ""$prefix"_ref".inv_1000
	rm ""$prefix"_ref".inv_1000_neighbours_agg
	rm ""$prefix"_ref".coords_dups_50
	rm ""$prefix"_ref".coords_dups_1000
	rm ""$prefix"_query".coords_dups_50
	rm ""$prefix"_query".coords_dups_1000
	rm $prefix.1000bp_inversionsclean
	rm $prefix.50bp_doubles
	rm $prefix.1000bp_filteredcalls
	rm $prefix.filteredcalls
	if [[ $blast_step == "yes" ]]
	then
		rm $prefix.insertions
		rm $prefix.deletions
		rm $prefix.blast_in
		rm $prefix.blast_del
		rm $prefix.deletion_blast
		rm $prefix.insertion_blast
	fi
fi

mkdir $prefix'_output'
mv $alignments_folder $prefix'_output'/
mv $prefix'.'* $prefix'_output'/
