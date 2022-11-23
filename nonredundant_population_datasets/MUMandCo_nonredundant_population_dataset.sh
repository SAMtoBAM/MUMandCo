## This the the way datasets were merged for O'Donnell et al. 2022 in order to create a non-redundant dataset
## In this example most strains contain multiple genomes per strain and therefore in order to verify events, each needed to be found at least twice
## This lead to ~95% of all variants being validated for strains with multiple genomes
## For those strains with only a single publically available genomes, fewer events were validated with the proportion validated depending on the genetic distance of the closest strain


## For the raw input; each genome file and their contigs/scaffolds (in fasta) were renamed to reflect the strain and assembler (in order to generate unique headers)
## e.g. For strain AAB, the canu assembly, the fasta file as renamed 'AAB_canu.fa' and it contained fasta headers such as '>AAB_canu_chrI' and '>AAB_canu_chrXII'
## This therefore allowed it to have a unique name 'AAB_canu' compared to the SMARTdenovo assembly or any other strain
## All genomes were then placed in a single directory, here names 'genomes'


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




