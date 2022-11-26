## This the the way datasets were merged for O'Donnell et al. 2022 in order to create a non-redundant dataset
## In this example most strains contain multiple genomes per strain and therefore in order to verify events, each needed to be found at least twice
## This lead to ~95% of all variants being validated for strains with multiple genomes
## For those strains with only a single publically available genomes, fewer events were validated with the proportion validated depending on the genetic distance of the closest strain


## For the raw input; each genome file and their contigs/scaffolds (in fasta) were renamed to reflect the strain and assembler (in order to generate unique headers)
## e.g. For strain AAB, the canu assembly, the fasta file as renamed 'AAB_canu.fa' and it contained fasta headers such as '>AAB_canu_chrI' and '>AAB_canu_chrXII'
## This therefore allowed it to have a unique name 'AAB_canu' compared to the SMARTdenovo assembly or any other strain
## All genomes were then placed in a single directory, here names 'genomes'

##for the phased genomes there is two different schemes depending on ploidy
##For diploids, there is two additional genomes assembled, named HP1 and HP2 e.g. "AAB_HP1_canu"
##For polyploids they were phased into individual blocks and assembled seperately, additionally they are not combined arbitraraily across blocks and/or chromosomes into a haplotype assembly
##Therefore for these genomes although they are all combined into an assembly labelled HP, these blocks assemblies were split into seperate files e.g. "AAB_block1_canu"


## STEP ONE: MUM&Co loops for raw SV calls per genome
## First step is to therefore generate an output folder for all the MUM&Co runs
mkdir MAC_output
cd MAC_output
## Now just a simple loop to run the MUMandCo script in the intitial directory for each genome
## The strain then takes everything prior to '.fa' in the genome file
## The SGD reference is also placed in the initial directory
ls ../genomes/ | while read genome
do
	strain=$( echo $genome | awk -F ".fa" '{print $1}' )
	echo $strain
	bash ../mumandco_v3.8.sh -r ../SGDref.nuclear_genome.tidy.fa -q ../genomes/${genome} -g 11000000 -o ${strain} -t 20
done
cd ..
## This gives us our raw SV calls from MUMandCo
## The next step is to therefore start comparing them suing some simple rules. This is outline in 02.



## STEP TWO: Overlapping raw calls to generate clusters and eliminating calls that are not validated by at least two genomes
mkdir MAC_overlap
cd MAC_overlap

## Remove any previous results
rm genomes_all.SVs_all_referencebased_strainID*
## Simplify raw calls to only the reference based data and concatenate all calls into a single file (therefore adding the strain tag in order to keep track of where it comes from)
ls ../MAC_output/ | sed 's/_output//g' | while read strain
do 
cat ../MAC_output/${strain}_output/${strain}.SVs_all.withfragment.tsv | tail -n +2 | awk -v strain="$strain" '{print $1"\t"$3"\t"$4"\t"$5"\t"$6"\t"strain}' >> genomes_all.SVs_all_referencebased_strainID.tsv
cat ../MAC_outpu/${strain}_output/${strain}.SVs_all.withfragment.tsv | tail -n +2 | awk -v strain="$strain" '{print $0"\t"strain}' >> genomes_all.SVs_all_referencebased_strainID.withfragment.tsv
done
## Remove any calls that are not translocations (the fragment is the entire translocated fragment) with scaffolds in them, i.e with Ns
cat genomes_all.SVs_all_referencebased_strainID.withfragment.tsv | awk '{if($9 !~ "NNNN" && $6 != "transloc"){print } else if($6 == "transloc") {print}}' > genomes_all.SVs_all_referencebased_strainID.withfragment.noNs.tsv
## Remove the rDNA region in the SGD S288c genome assembly
## TO BE AVOIDED IF USING ANY OTHER REFERENCE ASSEMBLY!
cat genomes_all.SVs_all_referencebased_strainID.withfragment.noNs.tsv | awk '{if($1 == "chrXII" && $2 >= 450000 && $2 <= 490000 || $1 == "chrXII" && $3 >= 450000 && $3 <= 490000){} else {print}}' > genomes_all.SVs_all_referencebased_strainID.withfragment.noNs.norDNA.tsv

## Modify input format for easier comparisons, partcularly for translocations
cat genomes_all.SVs_all_referencebased_strainID.withfragment.noNs.norDNA.tsv | awk '{if($6=="transloc") {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$10} else {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}}' | awk '{if($6 == "transloc") {gsub(/\[/,"",$11) ; gsub(/\]/,"",$11); gsub(/\:/,"\t",$11)  ;print $0} else {print $0}}' | sed 's/ /\t/g' > genomes_all.SVs_all_referencebased_strainID.withfragment.noNs.norDNA.input_format.tsv
## Remove any previous potential run files and create new folders/files to be used
rm -r clusters/
mkdir clusters

## Run a for loop for each window size, here just using the windows selected in the final nonredundant dataset (6500bp for inversions and translocations and 300bp for the rest)
## Can multi-thread this loop in order to speed it up as each individual chromosome or window can run individiually and get combined afterwards
for i in 300 6500
do
mkdir clusters/window${i}
## Split into chromosome and run seperately
cat genomes_all.SVs_all_referencebased_strainID.withfragment.noNs.norDNA.tsv | awk '{print $1}' | sort -u | while read chromosome
do
echo ${i}" "${chromosome}
rm matched.${chromosome}.${i}.tsv
touch matched.${chromosome}.${i}.tsv
rm unmatched.${chromosome}.${i}.tsv
rm nonredundant_calls.${chromosome}.${i}.tsv
## Grab events and read call by call
cat genomes_all.SVs_all_referencebased_strainID.withfragment.noNs.norDNA.input_format.tsv | awk -v chromosome="$chromosome" '{if($1 == chromosome)print}' | while read line
do
#check if line already has been placed in the matches.tsv file to make sure it has not already been placed in a cluster from previous loops
line2=$( echo $line | sed 's/ /\t/g' | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$11"\t"$12}' )
if grep --quiet "${line2}" matched.${chromosome}.${i}.tsv
then
echo "already matched"
else
## Assign necessary variables
chr=$( echo $line | awk '{print $1}' )
start=$( echo $line | awk '{print $3}'  )
end=$( echo $line | awk '{print $4}'  )
size=$( echo $line | awk '{print $5}' )
type=$( echo $line | awk '{print $6}'  )
## Assign additional variables for translocations
chrb=$( echo $line | awk '{print $11}' )
posb=$( echo $line | awk '{print $12}' )
## Search for matches using the window size etc
cat genomes_all.SVs_all_referencebased_strainID.withfragment.noNs.norDNA.input_format.tsv | grep $type | awk -v chromosome="$chromosome" '{if($1 == chromosome)print}' |\
awk -v window="$i" -v chr="$chr" -v start="$start" -v end="$end"  -v size="$size" -v type="$type" -v chrb="$chrb" -v posb="$posb" -F "\t" \
'BEGIN{OFMT="%.f" ; CONVFMT="%.f"} \
{if(type=="insertion" && chr==$1 && $3>=start-window && $3<=start+window && $4>=end-window && $4<=end+window && type==$6 && $5>=size*0.9 && $5<=size*1.1 ||\
type=="transloc" && chr==$1 && $3>=start-window && $3<=start+window && type==$6 && chrb==$11 && $12>=posb-window && $12<=posb+window ||\
type=="transloc" && chr==$1 && $4>=end-window && $4<=end+window && type==$6 && chrb==$11 && $12>=posb-window && $12<=posb+window ||\
type=="inversion" && chr==$1 && $3>=start-window && $3<=start+window && $4>=end-window && $4<=end+window && type==$6 ||\
type!="insertion" && type!="transloc" && type!= "inversion" && chr==$1 && $3>=start-window && $3<=start+window && $4>=end-window && $4<=end+window && type==$6)\
{print $0}}' > current_matches.${chromosome}.${i}.tsv
#see how many matches are found minus the call itself
count=$( cat current_matches.${chromosome}.${i}.tsv | wc -l | awk '{print $1-1}' )
#echo "currently "$count" matches"
## See if there were any matches placed in the newly written current_matches.${chromosome}.${i}.tsv file
## if no matches, skip next step
## else loop through finding matches until no new lines/matches are added to the current_matches.${chromosome}.${i}.tsv file
if [[ $count == "0" ]]
then
#echo "no matches, moving on"
## Save unmatched line split by type
echo $line | sed 's/ /\t/g' >> unmatched.${chromosome}.${i}.tsv
else
## Begin counting the cluster
n=$( echo $n | awk '{print $1+1}'  )
#echo "found matches, cluster " $n
## Reset variables
nloop="0"
wc2=""
wc=$( cat current_matches.${chromosome}.${i}.tsv | wc -l )
until [[ ${wc}  == ${wc2} ]]
do
## Count number of loops performed for new matches being added to the cluster
nloop=$( echo $nloop | awk '{print $0+1}')
#echo "running loop "$nloop
## Get intial count at beginning of loop
wc=$( cat current_matches.${chromosome}.${i}.tsv | sort -u | wc -l  )
#echo "starting with "$wc" paired calls"
## Going to use the max and min values to search for more matches and try combinations of them therefore reducing the number of searches when many matches are found
startmax=$( cat current_matches.${chromosome}.${i}.tsv | awk '{if(start<$3 || start==""){start=$3}} END{print start}' )
startmin=$( cat current_matches.${chromosome}.${i}.tsv | awk '{if(start>$3 || start==""){start=$3}} END{print start}' )
endmax=$( cat current_matches.${chromosome}.${i}.tsv | awk '{if(end<$4 || end==""){end=$4}} END{print end}' )
endmin=$( cat current_matches.${chromosome}.${i}.tsv | awk '{if(end>$4 || end==""){end=$4}} END{print end}' )
sizemax=$( cat current_matches.${chromosome}.${i}.tsv | awk '{if(size<$5 || size==""){size=$5}} END{print size}' )
sizemin=$( cat current_matches.${chromosome}.${i}.tsv | awk '{if(size>$5 || size==""){size=$5}} END{print size}' )
#max and min for translocation joint position
posbmax=$( cat current_matches.${chromosome}.${i}.tsv | awk '{if(posb<$12 || posb==""){posb=$12}} END{print posb}' )
posbmin=$( cat current_matches.${chromosome}.${i}.tsv | awk '{if(posb>$12 || posb==""){posb=$12}} END{print posb}' )
#echo $type $startmin $startmax $endmin $endmax $sizemin $sizemax $chrb $posbmin $posbmax
cat genomes_all.SVs_all_referencebased_strainID.withfragment.noNs.norDNA.input_format.tsv | grep $type | awk -v chromosome="$chromosome" '{if($1 == chromosome)print}' |\
awk -v window="$i" -v chr="$chr" -v startmax="$startmax" -v endmax="$endmax" -v sizemax="$sizemax" -v startmin="$startmin" -v endmin="$endmin"  -v sizemin="$sizemin" -v type="$type" -v chrb="$chrb" -v posbmax="$posbmax" -v posbmin="$posbmin" -F "\t" \
'BEGIN{OFMT="%.f" ; CONVFMT="%.f"} \
{if(type=="insertion" && chr==$1 && $3>=startmin-window && $3<=startmax+window && $4>=endmin-window && $4<=endmax+window && type==$6 && $5>=sizemin*0.9 && $5<=sizemax*1.1 ||\
type=="transloc" && chr==$1 && $3>=startmin-window && $3<=startmax+window && type==$6 && chrb==$11 && $12>=posbmin-window && $12<=posbmax+window ||\
type=="transloc" && chr==$1 && $4>=endmin-window && $4<=endmax+window && type==$6 && chrb==$11 && $12>=posbmin-window && $12<=posbmax+window ||\
type=="inversion" && chr==$1 && $3>=startmin-window && $3<=startmax+window && $4>=endmin-window && $4<=endmax+window && type==$6 ||\
type!="insertion" && type!="transloc" && type!="inversion" && chr==$1 && $3>=startmin-window && $3<=startmax+window && $4>=endmin-window && $4<=endmax+window && type==$6 )\
{print $0}}' >> current_matches.${chromosome}.${i}.tsv
## Remove doubles
cat current_matches.${chromosome}.${i}.tsv | sort -u > temp.${chromosome}.${i}
mv temp.${chromosome}.${i} current_matches.${chromosome}.${i}.tsv
## Get new count of matches to compare to initial count before loops
wc2=$( cat current_matches.${chromosome}.${i}.tsv | wc -l  )
#echo "ending with "$wc2" paired calls"
done
## Can begin averaging out the values to create the merged non-redundant call
startavg=$( cat current_matches.${chromosome}.${i}.tsv | awk 'BEGIN{OFMT="%.f" ; CONVFMT="%.f"} {sum=sum+$3; count=count+1} END{print sum/count}' )
endavg=$( cat current_matches.${chromosome}.${i}.tsv | awk 'BEGIN{OFMT="%.f" ; CONVFMT="%.f"} {sum=sum+$4; count=count+1} END{print sum/count}' )
sizeavg=$( cat current_matches.${chromosome}.${i}.tsv | awk 'BEGIN{OFMT="%.f" ; CONVFMT="%.f"} {sum=sum+$5; count=count+1} END{print sum/count}' )
#genomes=$(  cat current_matches.${chromosome}.${i}.tsv | awk '{print $2}' | tr '\n' ';' )
#startall=$(  cat current_matches.${chromosome}.${i}.tsv | awk '{print $3}' | tr '\n' ';' )
#endall=$(  cat current_matches.${chromosome}.${i}.tsv | awk '{print $4}' | tr '\n' ';' )
#sizeall=$(  cat current_matches.${chromosome}.${i}.tsv | awk '{print $5}' | tr '\n' ';' )
posbavg=$( cat current_matches.${chromosome}.${i}.tsv | awk 'BEGIN{OFMT="%.f" ; CONVFMT="%.f"} {sum=sum+$12; count=count+1} END{print sum/count}' )

#echo $chr $startavg $endavg $sizeavg $type
## Now print out the new non-redundant call with averaged positions and size
cat current_matches.${chromosome}.${i}.tsv | awk -v chr="$chr" -v type="$type" -v startavg="$startavg" -v endavg="$endavg" -v sizeavg="$sizeavg" -F "\t" \
'BEGIN{OFMT="%.f" ; CONVFMT="%.f"} \
{genomes=genomes";"$2 ; regions=regions";"$9 ; startall=startall";"$3 ; endall=endall";"$4 ; sizeall=sizeall";"$5 ; chrbposb=chrbposb";"$10 }\
END{if(type=="transloc") {print chr"\t"genomes"\t"startavg"\t"endavg"\t"sizeavg"\t"type"\t"regions"\t"startall"\t"endall"\t"sizeall"\t"chrbposb} else {print chr"\t"genomes"\t"startavg"\t"endavg"\t"sizeavg"\t"type"\t"regions"\t"startall"\t"endall"\t"sizeall}}' |\
awk -F "]" '{if($0 ~ "transloc") {print $0"\t"((NF-1)/2)} else {print $0}}' |\
awk -F "[" '{if($0 ~ "transloc") {print $0"\t"((NF-1)/2)} else {print $0}}' |\
awk -F "\t" -v chrb="$chrb" -v posbavg="$posbavg" 'BEGIN{OFMT="%f"} {if($6 == "transloc" && $13 > $12) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t["chrb":"posbavg"[\t"$11} \
else if($6 == "transloc" && $13 < $12) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t]"chrb":"posbavg"]\t"$11} \
else if($6 == "transloc" && $13 == $12) {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\tCONFLICT\t"$11} \
else {print $0}}' |\
sed "s/\t\;/\t/g" >> nonredundant_calls.${chromosome}.${i}.tsv


#add the current matches to the file in order to skip already processed lines
cat current_matches.${chromosome}.${i}.tsv | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$11"\t"$12}' >> matched.${chromosome}.${i}.tsv
cat matched.${chromosome}.${i}.tsv | sort -u > temp2.${chromosome}.${i}
mv temp2.${chromosome}.${i} matched.${chromosome}.${i}.tsv
#move cluster file to another position before re-writing
mv current_matches.${chromosome}.${i}.tsv clusters/window${i}/${n}_cluster_matches.${chromosome}.${i}.${type}.tsv

fi
fi
done
done


## Now we take the non-redundant calls paired per chromosome and window size and aggregate them per window then generate the final callset
## We can also calculate some stats for the raw and non-redundant data

for i in 300 6500
do

	cat nonredundant_calls.*.${i}.tsv > nonredundant_calls.$i.tsv
	rm nonredundant_calls.*.${i}.tsv
	cat unmatched.*.${i}.tsv > unmatched.$i.tsv
	rm unmatched.*.${i}.tsv
	cat matched.*.${i}.tsv > matched.$i.tsv
	rm matched.*.${i}.tsv
	rm current_matches.*.${i}.tsv

	echo "window: "$i
	echo " "
	total=$( cat genomes_all.SVs_all_referencebased_strainID.withfragment.noNs.norDNA.tsv | wc -l )
	echo "Total number of calls:" $total
	totalmatch=$( cat unmatched.${i}.tsv | wc -l | awk -v total="$total" '{print total-$0"\t("(total-$0)/total")"}' )
	echo "Total number of calls matched:" $totalmatch
	totalNR=$(cat nonredundant_calls.${i}.tsv |  wc -l)
	echo "Resulting number of non-redundant calls:" $totalNR
	deletion=$(cat nonredundant_calls.${i}.tsv | grep 'deletion' | wc -l)
	echo "non-redundant deletions:" $deletion
	insertion=$(cat nonredundant_calls.${i}.tsv | grep 'insertion' | wc -l)
	echo "non-redundant insertions:" $insertion
	duplication=$(cat nonredundant_calls.${i}.tsv | grep 'duplication' | wc -l)
	echo "non-redundant duplications:" $duplication
	contraction=$(cat nonredundant_calls.${i}.tsv | grep 'contraction' | wc -l)
	echo "non-redundant contractions:" $contraction
	inversion=$(cat nonredundant_calls.${i}.tsv | grep 'inversion' | wc -l)
	echo "non-redundant inversions:" $inversion
	transloc=$(cat nonredundant_calls.${i}.tsv | grep 'transloc' | wc -l)
	echo "non-redundant translocs:" $transloc
	echo "Number of genomes per SV"
	cat nonredundant_calls.${i}.tsv | awk '{print $2}' | awk -F ";" 'NF > m { m = NF } { s += NF } END { printf("Max = %d\nAvg = %g\n", m, s/NR) }'
	echo " "
done


## count the number of events per window,total, validated/mathced and non-redundant
echo "Calculating summary values such as number matched and non-redundant counts per window and per SV type"
echo "SV_type,window,total,matched,proportion,nonredundant" | sed 's/,/\t/g' > summary_counts.tsv

echo "" | awk '{if($6 != "") print $6} END{print "chr"}' - nonredundant_calls.0.tsv | sort -u | while read event
do
	#first of all the events called by MUM&Co
	all=$( cat genomes_all.SVs_all_referencebased_strainID.withfragment.noNs.norDNA.tsv | grep $event | wc -l )

	#then by all the events that found a match and the number of nonredundant
	for i in 300 6500
	do
		unmatched=$( cat unmatched.$i.tsv | grep ${event} | awk '{print $2}' | awk -F ";" '{sum=sum+NF} END{print sum}' )
		nonredundant=$( cat nonredundant_calls.$i.tsv | grep ${event} | wc -l )
	
		echo ${event}","${i}","${all}","${unmatched}","${nonredundant} |\
			sed 's/chr/all/g' |\
			sed 's/,/\t/g' |\
			awk '{print $1"\t"$2"\t"$3"\t"$3-$4"\t"($3-$4)/$3"\t"$5}' >> summary_counts.tsv
	done
done

echo "Combining calls from the 300 window (deletions, insertions, duplications and contractions) and 6500 window (inversions and translocations) into a single file"
## Need to get variants from different window files and generate the final file and sort it again by chr and position
## It is a pain to sort using 'sort' with chr and position and all the other columns so broken up by a chr list and sorted by position for each chromsome seperately
cat nonredundant_calls.300.tsv | awk '{if($6 == "deletion" || $6 == "contraction" || $6 == "insertion" || $6 == "duplication") print $0}' >> SVs.temp.tsv
cat nonredundant_calls.6500.tsv | awk '{if($6 == "transloc" || $6 == "inversion") print $0}' >> SVs.temp.tsv
rm SVs.manualmods.tsv
cat SVs.temp.tsv | awk '{print $1}' | sort -u | while read chr
do
	cat SVs.temp.tsv | awk -v chr="$chr" '{if($1 == chr) print}' | sort -k3n -k4n >> SVs.manualmods.tsv
done
rm SVs.temp.tsv

echo "Moving preliinary files into folders"

##move the raw files into seperate folders
mkdir raw_clustering/
mv matched.*.tsv raw_clustering/
mv unmatched.*.tsv raw_clustering/
mv clusters/ raw_clustering/
mkdir raw_clustering_NR/
mv nonredundant_calls.*.tsv raw_clustering_NR/
mkdir raw_clustering_summaryfiles/
mv calls_within_or_across_strains.* raw_clustering_summaryfiles/
mv genome_matches_proportion.* raw_clustering_summaryfiles/
mv summary_counts.tsv raw_clustering_summaryfiles/

echo "Splitting file into two files. One with transocations with 'conflicts' and the other contains the remaining events"

####need to manually edit the final file a bit for 'conflicting' translocations and large deletions from polyploid blocks
###then can analyse this file
### translocations are in the conflict file, where the sense/orientation of the ajoining fragment to a translocation fragment was evenly waited in both directions
### here they are manually modified, split if containing more than two assemblies, removed if only 2
## after manual mod they will hav a 'DONE' tag in the file name
cat SVs.manualmods.tsv | grep 'CONFLICT' > SVs.manualmods.conflict.tsv
### the other contains the large deletions that are also manually modified, and checked against changed in coverage used for aneuploidy detection
cat SVs.manualmods.tsv | grep -v 'CONFLICT' > SVs.manualmods.noconflict.tsv

##these files then need to be combined back together
echo "MANUAL MODIFICATION IS NOW NEEDED OF THE FILES 'SVs.manualmods.conflict.tsv' AND 'SVs.manualmods.noconflict.tsv' "
echo "'conflict' contains translocations with equal number fo merged variants for each orientation of the adjacent translocation fragment"
echo "after manual assessing, the manually curated and solved CONFLICT translocations need to be put in a file called: 'SVs.manualmods.conflict.DONE.tsv'"
echo "'noconflict' needs to be manually assessed for the validity of large deletions, primarily within polyploid phased blocks. Check 03.after... for definitions. Manual curation is needed paticularly for checking large deletions against coverage analysis."
echo "after manual assessing, the manually curated and solved 'noconflict' events need to be put in a file called: 'SVs.manualmods.noconflict.toremove.tsv'"


##ran these two commands to get the non-conflicting translocations (removing most and just splitting one into two different calls)
##grab only the sixth call, the first two aggregated were used to make another call and the 3rd and 4th for the other call
##re-calculated the averages for the positions etc automatically
##EXCEPT the position for the translocation joined region, was done manually and implicity written into the output
cat SVs.manualmods.conflict.tsv | awk '{if(NR == "6") print}' | sed 's/\t/\n/g' | awk -F ";" '{if($2 != "" ){print $1";"$2} else {print $1}}' | tr '\n' '\t' | sed 's/CONFLICT/\[chrXIII:809444[/g' > temp1.conflict
start=$( cat temp1.conflict | cut -f8 | sed 's/;/\n/g' | awk 'BEGIN{OFMT="%.f" ; CONVFMT="%.f"} {sum=sum+$1} END{print sum/NR}'  )
end=$( cat temp1.conflict | cut -f9 | sed 's/;/\n/g' | awk 'BEGIN{OFMT="%.f" ; CONVFMT="%.f"} {sum=sum+$1} END{print sum/NR}'  )
size=$( cat temp1.conflict | cut -f10 | sed 's/;/\n/g' | awk 'BEGIN{OFMT="%.f" ; CONVFMT="%.f"} {sum=sum+$1} END{print sum/NR}'  )
cat temp1.conflict | awk -v start="$start" -v end="$end" -v size="$size" '{print $1"\t"$2"\t"start"\t"end"\t"size"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' >> SVs.manualmods.conflict.DONE.tsv

cat SVs.manualmods.conflict.tsv | awk '{if(NR == "6") print}' | sed 's/\t/\n/g' | awk -F ";" '{if($2 != "" ){print $3";"$4} else {print $1}}' | tr '\n' '\t' | sed 's/CONFLICT/\]chrXIII:806016]/g' > temp2.conflict
start=$( cat temp2.conflict | cut -f8 | sed 's/;/\n/g' | awk 'BEGIN{OFMT="%.f" ; CONVFMT="%.f"} {sum=sum+$1} END{print sum/NR}'  )
end=$( cat temp2.conflict | cut -f9 | sed 's/;/\n/g' | awk 'BEGIN{OFMT="%.f" ; CONVFMT="%.f"} {sum=sum+$1} END{print sum/NR}'  )
size=$( cat temp2.conflict | cut -f10 | sed 's/;/\n/g' | awk 'BEGIN{OFMT="%.f" ; CONVFMT="%.f"} {sum=sum+$1} END{print sum/NR}'  )
cat temp2.conflict | awk -v start="$start" -v end="$end" -v size="$size" '{print $1"\t"$2"\t"start"\t"end"\t"size"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12}' >> SVs.manualmods.conflict.DONE.tsv
rm temp*.conflict

###first take the list of 'non conflict' events that are to be removed
###these included Deletions and contractions greater than 20kb that were not verified by coverage changes
###the list 'SVs.manualmods.noconflict.toremove.tsv' only needs to kept the first 6 identifying rows of each call, which will make each line unique (from chromosome to SV type)
###the above checking needs to be done manually, the rest will be done automatically, adding the calls to be removed to those identified as above
###the rest being done automatically:
###remove all insertions with a a gap between locations (end-start) bigger than 25kb
##this removed about 210 events (~4%) of the calls based on first and currently only run (value may change but not likely)

cat SVs.manualmods.noconflict.tsv | awk '{if($6 == "insertion" || $6 == "duplication")print}' | awk '{if($4-$3 > 25000) print}' | cut -f1-6 >> SVs.manualmods.noconflict.toremove.tsv 


cp SVs.manualmods.noconflict.tsv SVs.manualmods.noconflict.temp.tsv
cat SVs.manualmods.noconflict.toremove.tsv | while read line
do
	chr=$( echo $line | awk '{print $1}'  )
	start=$( echo $line | awk '{print $3}'  )
	end=$( echo $line | awk '{print $4}'  )
	size=$( echo $line | awk '{print $5}'  )
	type=$( echo $line | awk '{print $6}'  )
	cat SVs.manualmods.noconflict.temp.tsv | awk -v chr="$chr" -v start="$start" -v end="$end" -v size="$size" -v type="$type" '{if($1==chr && $3 == start && $4 == end && $5 == size && $6 == type) {print $0"\tmatched"} else {print $0}}' > SVs.manualmods.noconflict.temp2.tsv
	mv SVs.manualmods.noconflict.temp2.tsv SVs.manualmods.noconflict.temp.tsv
done
grep -v matched SVs.manualmods.noconflict.temp.tsv > SVs.manualmods.noconflict.DONE.tsv
rm SVs.manualmods.noconflict.temp.tsv

###the other file was the 'conflict' file where the orientation of joining events in translocations were not the same in merged events
### these were either split into different events or removed if they had more or less than two genomes for the call respectively
### this file can now be added to the other filtered file to get the final SVs dataset!

#first just merge the manually modified files and sort them by start position
cat SVs.manualmods.noconflict.DONE.tsv SVs.manualmods.conflict.DONE.tsv | sort -k3V > temp.tsv
#then explcitely encode the chromosome order since roman numeral are a pain to order properly
#and grab the SVs associated with each chromosome in this order and add them to the final SV file
echo chrI,chrII,chrIII,chrIV,chrV,chrVI,chrVII,chrVIII,chrIX,chrX,chrXI,chrXII,chrXIII,chrXIV,chrXV,chrXVI | tr ',' '\n' | while read chr
do
	cat temp.tsv | awk -v chr="$chr" '{if($1 == chr) print $0}'
done > SVs.FINAL.tsv
rm temp.tsv


##now can convert the final dataset into a vcf; first containing just a simplified homozygous call for all
rm -r nonredundant_VCF
mkdir nonredundant_VCF

##Can mimic the header from any generic MUM&Co VCF but leave out the query genome parts
##Just steal from AAB_canu output except for the column header line as this needs strains added
date=$( date +%Y%m%d )
grep '##' ../MAC_output/AAB_canu_output/AAB_canu.SVs_all.vcf | grep -v query | awk -v date="$date" '{print} END{print "##source "date"=In-house bash script aggregating calls based on each genomes MUM&Co tsv output then converting to VCF format and aggregating by strain" }' > nonredundant_VCF/generic_header.txt
##so now take the column header and add the strain names alongside it
grep "#" ../MAC_output/AAB_canu_output/AAB_canu.SVs_all.vcf | grep -v "##" | awk -F "\t" '{for(i=1; i <= NF; i++) {if($i !~ "AAB") {print $i}}}' | cat - all_strains.txt | tr '\n' '\t' | awk '{print $0}' > nonredundant_VCF/generic_header2.txt

##Big issue is what becomes the non-ref allele?
##For those taken from the ref it is fine (e.g. everything but insertions and contractions)
##but the insertion allele? a consensus of the insertion sequences for all strains??

###file in nonredundant_vcf/generating_multisample_vcf has the next commands for converting the final tsv file to a multi sample vcf
rm nonredundant_VCF/homoonly.temp.vcf
cat SVs.FINAL.tsv | while read variant
do
##for the reference the fragment is reextracted using the new coordinates
##this is harder for the insertions and contraction as this was taken from the query
##but here we have many query genomes...
##so for now I am just taking the first fragment present in the first genome
type=$( echo $variant | awk '{print $6}' )
if [[ $type == "insertion" || $type == "contraction" ]]
then
fragment=$( echo $variant | awk '{print $7}' | awk -F ";" '{print $1}' )
else
location=$( echo $variant | awk '{if($3 > $4) {print $1":"$4"-"$3} else {print $1":"$3"-"$4}}' )
fragment=$( samtools faidx ../SGDref.nuclear_genome.tidy.fa $location | grep -v ">" | tr -d '\n')
fi

straincount=$( cat all_strains.txt | wc -l )

## create some fill 0/0 place keepers to add to the line
fill=$( printf "%0.s0/0;" $( seq 1 $straincount  ) | sed 's/;/\t/g')

##generate a generic vcf set for the variant and then add the fill above to the end
templine=$( echo $variant","$fragment | sed 's/,/\t/g' | awk -v fill="$fill" '{if($6 == "deletion") {print $1"\t"$3"\t.\t"$11"\t<DEL>\t.\tPASS\tEND="$4";SVLEN=-"$5";SVTYPE=DEL\tGT\t"fill}\
if($6 == "insertion") {print $1"\t"$3"\t.\t<INS>\t"$11"\t.\tPASS\tEND="$4";SVLEN="$5";SVTYPE=INS\tGT\t"fill}\
if($6 == "duplication") {print $1"\t"$3"\t.\t"$11"\t<DUP>\t.\tPASS\tEND="$4";SVLEN="$5";SVTYPE=DUP\tGT\t"fill}\
if($6 == "contraction") {print $1"\t"$3"\t.\t"$11"\t<CONTR>\t.\tPASS\tEND="$4";SVLEN="$5";SVTYPE=CONTR\tGT\t"fill}\
if($6 == "inversion") {print $1"\t"$3"\t.\t"$11"\t<INV>\t.\tPASS\tEND="$4";SVLEN="$5";SVTYPE=INV\tGT\t"fill}\
if($6 == "transloc") {print $1"\t"$3"\t.\t"$13"\t"$11"\t.\tPASS\tEND="$4";SVLEN="$5";SVTYPE=TRA\tGT\t"fill}}' ) 


##here I will change the temporary 0/0 fill with 1/1 for each strains that contain each variant
##to get the right column to modify I just ask take the header which the strain names and find which column the strain would be
##then I only modify that per single strain and save over the previous temp file
## the brackets before the while loop and after the file being saved is to maintain the variable changes made in the while loop until I save it
echo $variant | awk '{print $2}' | awk -F ";" '{for(i=1; i<=NF; i++) print $i}' | awk -F "_canu" '{print $1}' | awk -F "_Sdn" '{print $1}' | awk -F "_mixed" '{print $1}' | awk -F "_chr" '{print $1}' | awk -F "_block" '{print $1}' | awk -F "_HP" '{print $1}' | sort -u | ( while read strain
do
count=$( cat nonredundant_VCF/generic_header2.txt | tr '\t' '\n' | awk -v strain="$strain" '{if($0 == strain) print NR}'  )
templine2=$( echo $templine | awk -v count="$count" '{for(i=1; i<=NF; i++) if(i == count) {print "1/1"} else {print $i}}' | tr '\n' '\t'  )

templine=$( echo $templine2 | sed 's/ /\t/g' )
done 

echo $templine | sed 's/ /\t/g' >> nonredundant_VCF/homoonly.temp.vcf )

done 

cat nonredundant_VCF/generic_header.txt nonredundant_VCF/generic_header2.txt nonredundant_VCF/homoonly.temp.vcf | sed 's/\t$//g' > nonredundant_VCF/homoonly.vcf
rm homoonly.temp.vcf

###one last step is to just sort the file to make sure it will work with other downstream tools
bcftools sort nonredundant_VCF/homoonly.vcf > nonredundant_VCF/homoonly2.vcf
mv nonredundant_VCF/homoonly2.vcf nonredundant_VCF/homoonly.vcf

###that gives a VCF file for all the strains with homozygous labels for every call

###In order to make the next step to add heterozygosity to the vcf we need to calculate whether variants are heterowygous or not!
##because this is using genome comparisons only phased genomes can possibly have heterozygous variants called
##in this dataset we have two formats, diploids and polyploids

## calculate those calls from heterozygote diploids that are present in only a single haplotype (hetero) or both haplotypes (homo) whether within one or multiple assemblies for each
# current example is just using the 0 window output raw file
rm -r diploid_zygosity
mkdir diploid_zygosity
##also adding the HP involved in the hetero calls too, added to the end of lines with hetero calls
cat genomes_all.SVs_all_referencebased_strainID.withfragment.noNs.norDNA.tsv  | grep 'HP'  | awk '{print $2}' | sort -u | awk -F "_" '{print $1}' | sort -u | while read strain
do
cat SVs.FINAL.tsv | awk -v strain="$strain" '{if($2 ~ strain"_HP" ) print}' | while read call
do
zygosity=$( echo $call |  awk '{print $2}' | awk  -v strain="$strain" -F ";" '{for(i=1; i <= NF; i++ ) if($i ~ "HP" && $i ~ strain) print $i}' | awk -F "_" '{print $2}' | sort -u  | awk 'END{if(NR > 1) {print "homo"} else if( $0 ~ "HP1") {print "hetero\tHP1"} else {print "hetero\tHP2"}}' )
echo ${call} ${zygosity} | sed 's/ /\t/g' >> diploid_zygosity/${strain}.calls_zygosity.FINAL.tsv
done
done
##get the type of call and proportion homo and hetero and for the hetero calls the proportion from HP1 and HP2
echo "strain,SV_type,total,homozygous,homozygous_proportion,heterozygous,heterozygous_proportion,HP1_proportion,HP2_proportion" | sed 's/,/\t/g' > diploid_zygosity/calls_by_type_homo_hetero_count.tsv
ls diploid_zygosity/  | grep ".calls_zygosity" | while read file
do
strain=$( echo $file | awk -F "." '{print $1}'  )
cat diploid_zygosity/$file | awk '{print $6} END{print "chr"}' | sort -u | while read event
do
homo=$( grep $event diploid_zygosity/$file | awk '{if($NF == "homo") print}' | wc -l  )
hetero=$( grep $event diploid_zygosity/$file | awk '{step=NF-1; if($step == "hetero") print}' | wc -l)
HP1=$( grep $event diploid_zygosity/$file | awk '{step=NF-1; if($step == "hetero") print $NF}' | grep HP1 | wc -l  )
HP2=$( grep $event diploid_zygosity/$file | awk '{step=NF-1; if($step == "hetero") print $NF}' | grep HP2 | wc -l  )
grep $event diploid_zygosity/$file | wc -l | awk -v strain="$strain" -v event="$event" -v homo="$homo" -v hetero="$hetero" -v HP1="$HP1" -v HP2="$HP2" '{if(event == "chr") {event = "all" ; print strain"\t"event"\t"$0"\t"homo"\t"homo/$0"\t"hetero"\t"hetero/$0"\t"HP1/hetero"\t"HP2/hetero} else {print strain"\t"event"\t"$0"\t"homo"\t"homo/$0"\t"hetero"\t"hetero/$0"\t"HP1/hetero"\t"HP2/hetero}}'
done
done >> diploid_zygosity/calls_by_type_homo_hetero_count.tsv

### calculating whether events are heterozygous for polyploids
### more difficult as for each region in the reference genome I don't know how many haplotypes have been assembled
### therefore it is hard to say one is missing for an event and therefore it is heterozygous etc
### however I have aligned all the 'best' assembly blocks, in a single file, against the ref and calculated the coverage
### This should tell me how many blocks there are and therefore how many should be present for a homozygous call
### will use this info to detect heterozygous/homozygous calls
##get the strains
rm -r polyploid_zygosity
mkdir polyploid_zygosity
##first we need to actually calculate the coverage for each region which is an indicator of the number of haplotypes for each region
##here only the number of unique haplotypes have been assembled do to using nphase phasing and assembling only the blocks
##this is not the same method as used for diplods as for them the unphased reads were added to each haplotype set in order to fix assembly issues that would be created in the face of LOH
mkdir polyploid_zygosity/polyploid_block_coverage/
ls ../genomes | grep "_block" | awk -F "_block" '{print $1}'  | while read strain
do
	cat ../genomes/${strain}_block* > polyploid_zygosity/polyploid_block_coverage/${strain}_HP.fa
	strain=$( echo $genome | awk -F "." '{print $1}'  )
	minimap2 -ax asm5 --secondary=no ../SGDref.nuclear_genome.tidy.fa polyploid_zygosity/polyploid_block_coverage/${strain}_HP.fa | samtools sort -o polyploid_zygosity/polyploid_block_coverage/${strain}.ref_minimap.bam -
	bedtools genomecov -d -ibam polyploid_zygosity/polyploid_block_coverage/${strain}.ref_minimap.bam | gzip > polyploid_zygosity/polyploid_block_coverage/${strain}.ref_minimap.cov.tsv.gz
	zcat polyploid_zygosity/polyploid_block_coverage/${strain}.ref_minimap.cov.tsv.gz | awk '{if($3 > 0) print $3}' | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }' | awk -v strain="$strain" '{print strain"\t"$0}' >> all_strains.cov_median.tsv
done

##now we can determine if each call is homo or heterowygote considering we know the number of haplotypes
ls ../genomes | grep "_block" | awk -F "_block" '{print $1}' | while read strain
do
	cat SVs.FINAL.tsv | awk -v strain="$strain" '{if($2 ~ strain"_block" ) print}' | while read call
	do
		chr=$( echo $call | awk '{print $1}' )
		start=$( echo $call | awk '{print $3}' )
		end=$( echo $call | awk '{print $4}' )
		type=$( echo $call | awk '{print $6}' )
		if [[ $type == "deletion" || $type == "contraction" ]]
		then
		## get 20kb up and downstream of the event start and stop since that region will have lower coverage then take the median coverage
		coverage=$( zcat polyploid_zygosity/polyploid_block_coverage/${strain}_HP.ref_minimap.cov.tsv.gz | awk -v chr="$chr" -v start="$start" -v end="$end" '{if($1 == chr && $2 >= start-20000 && $2 >= start || $1 == chr && $2 >= end && $2 <= end+20000) {print $0}}' | awk '{if($3 > 0) print $3}' | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }' | awk '{print $0}' )
		else
		## get 20kb up and downstream PLUS the region involved in the event then take the median coverage
		coverage=$( zcat polyploid_zygosity/polyploid_block_coverage/${strain}_HP.ref_minimap.cov.tsv.gz | awk -v chr="$chr" -v start="$start" -v end="$end" '{if($1 == chr && $2 >= start-20000 && $2 <= end+20000) {print $0}}' | awk '{if($3 > 0) print $3}' | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }' | awk '{print $0}' )
		fi
		zygosity=$( echo $call | awk '{print $2}' | awk  -v strain="$strain" -F ";" '{for(i=1; i <= NF; i++ ) if($i ~ "block" && $i ~ strain) print $i}' | awk -F "_" '{print $2}' | sort -u  | awk -v coverage="$coverage" 'END{if(NR < coverage) {print "hetero"} else {print "homo"}}' )
		echo $call $zygosity | sed 's/ /\t/g' >> polyploid_zygosity/${strain}.calls_zygosity.FINAL.tsv
	done
done
##gets stats for the zygosity of polyploid strain calls
##get the type of call and proportion homo and hetero
echo "strain,SV_type,total,homozygous,homozygous_proportion,heterozygous,heterozygous_proportion" | sed 's/,/\t/g' > polyploid_zygosity/calls_by_type_homo_hetero_count.tsv
ls polyploid_zygosity/  | grep ".calls_zygosity" | while read file
do
strain=$( echo $file | awk -F "." '{print $1}'  )
cat polyploid_zygosity/$file | awk '{print $6} END{print "chr"}' | sort -u | while read event
do
homo=$( grep $event polyploid_zygosity/$file | awk '{if($NF == "homo") print}' | wc -l  )
hetero=$( grep $event polyploid_zygosity/$file | awk '{if($NF == "hetero") print}' | wc -l)
grep $event polyploid_zygosity/$file | wc -l | awk -v strain="$strain" -v event="$event" -v homo="$homo" -v hetero="$hetero" '{if(event == "chr") {event = "all"; print strain"\t"event"\t"$0"\t"homo"\t"homo/$0"\t"hetero"\t"hetero/$0} else {print strain"\t"event"\t"$0"\t"homo"\t"homo/$0"\t"hetero"\t"hetero/$0}}'
done
done >> polyploid_zygosity/calls_by_type_homo_hetero_count.tsv


##combine diploid and polyploid heterozygosity estimation
##also add the SNPs amount per strains and the count and proportion of thos ethat are heterozygous
cat diploid_zygosity/calls_by_type_homo_hetero_count.tsv | cut -f1-7 | tail -n+2 > diploid.temp
cat polyploid_zygosity/calls_by_type_homo_hetero_count.tsv | tail -n+2 > polyploid.temp
cat polyploid_zygosity/calls_by_type_homo_hetero_count.tsv | cut -f1-7 | head -n1 | awk '{print $0"\tploidy\tSNPs_total\tSNPs_hetero\tSNPs_hetero_prop"}' > header.temp
cat diploid.temp polyploid.temp > diploid_polyploid.temp
cat diploid_polyploid.temp | while read line
do
strain=$( echo $line | awk '{print $1}' )
illumina=$(  grep $strain ../../illumina_variants.phenovar.stats.trim.tsv | awk '{print $5+$6"\t"$6"\t"$6/($5+$6)}' )
cat ../Table_nuclear_genome_stats.tsv | grep -v SGDref | awk -F "\t" -v strain="$strain" -v line="$line" '{if($3 == strain) {print line"\t"$5}}' | sort -u | awk -v illumina="$illumina" '{print $0"\t"illumina}'
done > diploid_polyploid.plus.temp
cat header.temp diploid_polyploid.plus.temp > SVs.calls_by_type_homo_hetero_count.tsv
rm *.temp

###now we can modify the vcf in order to incorperate the information on heterowygous calls
###just need to to diplay 0/1 instead of 1/1
cp nonredundant_VCF/homoonly.vcf nonredundant_VCF/homo_and_hetero.temp.vcf
ls *_zygosity/* | grep '.calls_zygosity' | while read file
do
strain=$( echo $file | awk -F "/" '{print $2}' | awk -F "." '{print $1}' )
cat $file | while read call
do
##see if the call is homozygous or heterozygous and make the variable the VCF notation associated with it
zygo=$( echo $call | awk '{ if($(NF-1) == "homo" || $NF == "homo") {print "1/1"} else if($(NF-1) == "hetero" || $NF == "hetero") {print "0/1"} }' )
##this is the trick to change the right call, it counts which column the strain occupies using the temp header containing column names
col=$( cat nonredundant_VCF/generic_header2.txt | awk -v strain="$strain" '{for(i=1; i <= NF; i++){if($i == strain) {print i}}}' )
##set of SV variables
chr=$( echo $call | awk '{print $1}' )
start=$( echo $call | awk '{print $3}' )
end=$( echo $call | awk '{print $4}' )
type=$( echo $call | awk '{print $6}'  | awk '{if($1 == "deletion") {print "DEL"} else if($1 == "insertion") {print "INS"} else if($1 == "duplication") {print "DUP"} else if($1 == "contraction") {print "CONT"} else if($1 == "inversion") {print "INV"} else if($1 == "transloc") {print "TRA"}}' )
size=$( echo $call | awk '{print $5}'  )
##just some read out
echo $strain $col $chr $start $end $size $type $zygo
## the col variable will be used to check that this cell contains currently a homozygous call, then if so, and all the other SV variables match, then it will change it according to the zygosity variable
cat nonredundant_VCF/homo_and_hetero.temp.vcf | grep -v "#" | awk -v col="$col" -v strain="$strain" -v chr="$chr" -v start="$start" -v end="$end" -v size="$size" -v type="$type" -v zygo="$zygo" '{if($col == "1/1" && $1 == chr && $2 == start && $8 ~ type && $8 ~ end";" && $8 ~ "SVLEN="size";" || $col == "1/1" && $1 == chr && $2 == start && $0 ~ type && $8 ~ end";" && $8 ~ "SVLEN=-"size";") {$col = zygo; print $0} else {print $0}}' | sed 's/ /\t/g' > nonredundant_VCF/homo_and_hetero.temp2.vcf
##move temp save to file to be continuously updated
mv nonredundant_VCF/homo_and_hetero.temp2.vcf nonredundant_VCF/homo_and_hetero.temp.vcf
done
done

###add the headers back
cat nonredundant_VCF/generic_header.txt nonredundant_VCF/generic_header2.txt nonredundant_VCF/homo_and_hetero.temp.vcf | sed 's/\t$//g' > nonredundant_VCF/homo_and_hetero.vcf
rm nonredundant_VCF/homo_and_hetero.temp.vcf

###Now I have the final dataset, they can be split by strain for a strain specific SV VCF
##can just use a loop with bcftools

###first need to compress with bgzip and create an index file with tabix
bgzip nonredundant_VCF/homo_and_hetero.vcf
tabix nonredundant_VCF/homo_and_hetero.vcf.gz


