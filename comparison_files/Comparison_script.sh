#!/bin/bash


prefix="show-diff_human_TRA_nofilt"
tools_vcf="show-diff/show-diff_human_TRA5"
species="human"
SV_type="TRA"
###TRANSOLCATIONS ACTUALLY DONE BY EYE###

##pick MUMandCO, paftools, assemblytics,show-diff,svmu
dataset="show-diff"

mkdir COMPARISONS
touch COMPARISONS/$prefix.summary.txt

if [[ $SV_type == "DEL" ]]
then
	cat $species/$SV_type/vcf/*.vcf |\
	tail -n +23 |\
	awk -F "\t" '{print $1"\t"$2"\t"$8}' | awk -F "=" '{print $1"\t"$4"\t"$2}' |\
	sed 's/;EVENT//g' |\
	awk '{print $1"\t"$2"\t"$4"\t"$5"\tblank"}' |\
	sed 's/DEL/deletion/g' > COMPARISONS/$species.$SV_type.filtered
else
	if [[ $SV_type == "DUPd"  ]]
	then
		cat $species/$SV_type/vcf/*.vcf |\
		tail -n +23 |\
		awk -F "\t" 'BEGIN{count="1"} {if(count=="1") {start=$5; count="2"} else if(count=="2") {print $1"\t"$2"\t"$2-1"\tinsertion\tblank"; count="1"}}' \
		> COMPARISONS/$species.$SV_type.filtered
	else
		if [[ $SV_type == "DUPt"  ]]
		then
			cat $species/$SV_type/vcf/*.vcf |\
			tail -n +23 |\
			awk -F "\t" 'BEGIN{count="1"} {if(count=="1") {start=$5; count="2"} else if(count=="2") {print $1"\t"start"\t"$2"\tinsertion\tblank"; count="1"}}' |\
			sed 's/:/\t/g' | awk '{print $1"\t"$3"\t"$4"\t"$5"\t"$6}' | sed 's/\[//g' \
			> COMPARISONS/$species.$SV_type.filtered
		else
			if [[ $SV_type == "INV"  ]]
			then
				cat $species/$SV_type/vcf/*.vcf |\
				tail -n +16 |\
				awk -F "\t" '{print $1"\t"$2"\t"$8}' | awk -F "=" '{print $1"\t"$2"\t"$6}' | awk '{print $1"\t"$2"\t"$5"\tinversion\tblank"}' \
				> COMPARISONS/$species.$SV_type.filtered
			else 
				if [[ $SV_type == "TRA"  ]]
				then
					cat $species/$SV_type/vcf/*.vcf |\
					tail -n +15 |\
					awk -F "\t" 'BEGIN{count="1"} {if(count=="1") {start=$5; count="2"} else if(count=="2") {print $1"\t"$2"\t"$2-1"\ttransloc\tblank"; count="1"}}' \
					> COMPARISONS/$species.$SV_type.filtered
				fi
			fi
		fi
	fi
fi

truth_vcf="COMPARISONS/$species.$SV_type.filtered"

if [[ $dataset == "MUMandCO" ]]
then
	cat $tools_vcf | cut -f 1-6 | sed 's/duplication/insertion/g' > COMPARISONS/$prefix.filtered
	tools_vcf_final=COMPARISONS/$prefix.filtered
else
	if [[ $dataset == "paftools" ]]
	then
		cat $tools_vcf | awk '{if($1=="V")print $2"\t"$9"\t"$3"\t"$4}' > $prefix.temp_final
		cat $tools_vcf |\
		sed 's/-/blank/g' |\
			awk '{if($1=="V") print $0}' |\
				awk '{if($7=="blank") {printf "\t"length($8)"\n"} else if($8=="blank") {printf length($7)"\t""\n"} else {printf length($7)"\t"length($8)"\n"}}' |\
					awk -F '\t' '{if($1>$2) {print $1-$2"\tdeletion"} else if($2>$1) {print $2-$1"\tinsertion"} else print $1"\tignore"}' > $prefix.temp_final2
		paste $prefix.temp_final $prefix.temp_final2 >> $prefix.temp_final3
		cat $prefix.temp_final3 | awk '{if($5 >=50) print $0}' > COMPARISONS/$prefix.filtered
		tools_vcf_final=COMPARISONS/$prefix.filtered
		rm $prefix.temp*
	else
		if [[ $dataset == "assemblytics" ]]
		then
			cat $tools_vcf | awk '{print $1"\t"$10"\t"$2"\t"$3"\t"$5"\t"$7}' |\
				sed 's/:.*+//g' | sed 's/:.*-//g' |\
				sed 's/Tandem_contraction/deletion/g' |\
				sed 's/Repeat_contraction/deletion/g' |\
				sed 's/Tandem_expansion/insertion/g' |\
				sed 's/Repeat_expansion/insertion/g' |\
				sed 's/Insertion/insertion/g' | sed 's/Deletion/deletion/g' > COMPARISONS/$prefix.filtered
			tools_vcf_final=COMPARISONS/$prefix.filtered
		else
			if [[ $dataset == "show-diff" ]]
			then
				cat $tools_vcf | tail -n +5 |\
					awk '{if($2== "GAP") {print $1"\t"$2"\t"$3"\t"$4"\t"$7} \
					else if( $2 == "INV" || $2 == "BRK" || $2 == "DUP" || $2 == "SEQ") {print $1"\t"$2"\t"$3"\t"$4"\t"$5} \
					else if( $2 == "JMP" || $5 >=50) {print $1"\tinsertion\t"$3"\t"$4"\t"$5} \
					else if( $2 == "JMP" || $5 <0) {print $1"\tinsertion\t"$4"\t"$3"\t"$5}}' |\
					awk '{if($2=="GAP" && $5>0) {print $1"\tblank\t"$3"\t"$4"\t"$5"\tdeletion"}\
					else if($2=="GAP" && $5<0) {print $1"\tblank\t"$3"\t"$4"\t"$5"\tinsertion"}\
					else if($2=="GAP" && $5==0) {print $1"\tblank\t"$3"\t"$4"\t"$5"\tNULL"}\
					else print $1"\tblank\t"$3"\t"$4"\t"$5"\t"$2}' |\
					sed 's/BRK/deletion/g' | sed 's/DUP/deletion/g' | sed 's/INV/inversion/g' | sed 's/SEQ/translocation/g' |\
					sed 's/-//g' |	awk '{if($5 >= 50 || $6 == "translocation" || $6 == "NULL" || $6 == "inversion") print $0}' > COMPARISONS/$prefix.filtered
				tools_vcf_final=COMPARISONS/$prefix.filtered
			else 
				if [[ $dataset == "svmu" ]]
				then
					cat $tools_vcf |\
					sed 's/DEL/deletion/g' |\
					sed 's/INS/insertion/g'|\
					sed 's/nCNV-Q/insertion/g'|\
					sed 's/CNV-Q/insertion/g'|\
					sed 's/nCNV-R/deletion/g'|\
					sed 's/CNV-R/deletion/g' |\
					sed 's/INV/inversion/g' |\
					awk '{if($4=="deletion" || $4=="insertion" || $4=="inversion" ) {print $1"\t"$5"\t"$2"\t"$3"\t"$9"\t"$4}}' > COMPARISONS/$prefix.filtered
					tools_vcf_final=COMPARISONS/$prefix.filtered
				fi
			fi
		fi
	fi
fi

mkdir COMPARISONS/${prefix}
if [[ $SV_type == "DUPt" ]]
then
	for i in 50 100 250 500 1000 2500 5000
	do
		size_of_window=$i
		buff=$(( size_of_window/2 ))
		
		cat $truth_vcf |\
			awk -v buff="$buff" '{print $1"\t"$2-buff"\t"$2+buff"\t"$3-buff"\t"$3+buff"\t"$4"\t"$5}' > ${truth_vcf}_buffer${size_of_window}.csv
		
		mkdir COMPARISONS/${prefix}_window_${size_of_window}
		cat ${truth_vcf}_buffer${size_of_window}.csv | while read line
		do
			count=$((count + 1))
			cat $tools_vcf_final |\
				awk -v line="$line" '{print $0"\t"line}' |\
					awk '{if($1==$7 && $3>=$8 && $3<=$9 && $4>=$10 && $4<=$11 && $6==$12 || $1==$7 && $4>=$8 && $4<=$9 && $3>=$10 && $3<=$11 && $6==$12) print $0}' > COMPARISONS/${prefix}_window_${size_of_window}/$count.txt
		done
		
		mv ${truth_vcf}_buffer* COMPARISONS/${prefix}
	
		cat COMPARISONS/${prefix}_window_${size_of_window}/*.txt | cut -f 1-6 | sort -k1,1 -k3,3V > COMPARISONS/${prefix}_${size_of_window}.csv
		
		total=$(wc -l COMPARISONS/${prefix}_${size_of_window}.csv)
		deletion=$(grep 'deletion' COMPARISONS/${prefix}_${size_of_window}.csv | wc -l)
		insertion=$(grep 'insertion' COMPARISONS/${prefix}_${size_of_window}.csv | wc -l)
		inversion=$(grep 'inversion' COMPARISONS/${prefix}_${size_of_window}.csv | wc -l)
		translocation=$(grep 'translocation' COMPARISONS/${prefix}_${size_of_window}.csv | wc -l)
			
		echo ${prefix}_SVs	$size_of_window
		echo deletion	$deletion
		echo insertion	$insertion
		echo inversion	$inversion
		echo translocation	$translocation
		echo ${prefix}_SVs	$size_of_window>>COMPARISONS/$prefix.summary.txt
		echo deletion	$deletion>>COMPARISONS/$prefix.summary.txt
		echo insertion	$insertion>>COMPARISONS/$prefix.summary.txt
		echo inversion	$inversion>>COMPARISONS/$prefix.summary.txt
		echo translocation	$translocation>>COMPARISONS/$prefix.summary.txt
	
		cat COMPARISONS/${prefix}_window_${size_of_window}/*.txt | cut -f 7-13 | sort -k1,1 -k3,3V |\
			awk -v buff="$buff" '{print $1"\t"$2+buff"\t"$4+buff"\t"$6"\t"$7}' > COMPARISONS/${prefix}_${size_of_window}_TRUTH.csv
		
		total_truth=$(wc -l COMPARISONS/${prefix}_${size_of_window}_TRUTH.csv)
		deletion_truth=$(grep 'deletion' COMPARISONS/${prefix}_${size_of_window}_TRUTH.csv | wc -l)
		insertion_truth=$(grep 'insertion' COMPARISONS/${prefix}_${size_of_window}_TRUTH.csv | wc -l)
		inversion_truth=$(grep 'inversion' COMPARISONS/${prefix}_${size_of_window}_TRUTH.csv | wc -l)
		translocation_truth=$(grep 'translocation' COMPARISONS/${prefix}_${size_of_window}_TRUTH.csv | wc -l)
		
		echo ${prefix}_Truth	$size_of_window
		echo deletion	$deletion_truth
		echo insertion	$insertion_truth
		echo inversion	$inversion_truth
		echo translocation	$translocation_truth
		echo ${prefix}_Truth	$size_of_window>>COMPARISONS/$prefix.summary.txt
		echo deletion	$deletion_truth>>COMPARISONS/$prefix.summary.txt
		echo insertion	$insertion_truth>>COMPARISONS/$prefix.summary.txt
		echo inversion	$inversion_truth>>COMPARISONS/$prefix.summary.txt
		echo translocation	$translocation_truth>>COMPARISONS/$prefix.summary.txt

		mv COMPARISONS/${prefix}_* COMPARISONS/$prefix/
	done
else
	for i in 50 100 250 500 1000 2500 5000
	do
		size_of_window=$i
		buff=$(( size_of_window/2 ))
		
		cat $truth_vcf |\
			awk -v buff="$buff" '{print $1"\t"$2-buff"\t"$2+buff"\t"$3-buff"\t"$3+buff"\t"$4"\t"$5}' > ${truth_vcf}_buffer${size_of_window}.csv
		
		mkdir COMPARISONS/${prefix}_window_${size_of_window}
		cat ${truth_vcf}_buffer${size_of_window}.csv | while read line
		do
			count=$((count + 1))
			cat $tools_vcf_final |\
				awk -v line="$line" '{print $0"\t"line}' |\
					awk '{if($1==$7 && $3>=$8 && $3<=$9 && $4>=$10 && $4<=$11 && $6==$12 || $1==$7 && $4>=$8 && $4<=$9 && $3>=$10 && $3<=$11 && $6==$12) print $0}' > COMPARISONS/${prefix}_window_${size_of_window}/$count.txt
		done
		
		mv ${truth_vcf}_buffer* COMPARISONS/${prefix}
	
		cat COMPARISONS/${prefix}_window_${size_of_window}/*.txt | cut -f 1-6 | sort| uniq | sort -k1,1 -k3,3V > COMPARISONS/${prefix}_${size_of_window}.csv
		
		total=$(wc -l COMPARISONS/${prefix}_${size_of_window}.csv)
		deletion=$(grep 'deletion' COMPARISONS/${prefix}_${size_of_window}.csv | wc -l)
		insertion=$(grep 'insertion' COMPARISONS/${prefix}_${size_of_window}.csv | wc -l)
		inversion=$(grep 'inversion' COMPARISONS/${prefix}_${size_of_window}.csv | wc -l)
		translocation=$(grep 'translocation' COMPARISONS/${prefix}_${size_of_window}.csv | wc -l)
			
		echo ${prefix}_SVs	$size_of_window
		echo deletion	$deletion
		echo insertion	$insertion
		echo inversion	$inversion
		echo translocation	$translocation
		echo ${prefix}_SVs	$size_of_window>>COMPARISONS/$prefix.summary.txt
		echo deletion	$deletion>>COMPARISONS/$prefix.summary.txt
		echo insertion	$insertion>>COMPARISONS/$prefix.summary.txt
		echo inversion	$inversion>>COMPARISONS/$prefix.summary.txt
		echo translocation	$translocation>>COMPARISONS/$prefix.summary.txt
	
		cat COMPARISONS/${prefix}_window_${size_of_window}/*.txt | cut -f 7-13 | sort| uniq | sort -k1,1 -k3,3V |\
			awk -v buff="$buff" '{print $1"\t"$2+buff"\t"$4+buff"\t"$6"\t"$7}' > COMPARISONS/${prefix}_${size_of_window}_TRUTH.csv
		
		total_truth=$(wc -l COMPARISONS/${prefix}_${size_of_window}_TRUTH.csv)
		deletion_truth=$(grep 'deletion' COMPARISONS/${prefix}_${size_of_window}_TRUTH.csv | wc -l)
		insertion_truth=$(grep 'insertion' COMPARISONS/${prefix}_${size_of_window}_TRUTH.csv | wc -l)
		inversion_truth=$(grep 'inversion' COMPARISONS/${prefix}_${size_of_window}_TRUTH.csv | wc -l)
		translocation_truth=$(grep 'translocation' COMPARISONS/${prefix}_${size_of_window}_TRUTH.csv | wc -l)
		
		echo ${prefix}_Truth	$size_of_window
		echo deletion	$deletion_truth
		echo insertion	$insertion_truth
		echo inversion	$inversion_truth
		echo translocation	$translocation_truth
		echo ${prefix}_Truth	$size_of_window>>COMPARISONS/$prefix.summary.txt
		echo deletion	$deletion_truth>>COMPARISONS/$prefix.summary.txt
		echo insertion	$insertion_truth>>COMPARISONS/$prefix.summary.txt
		echo inversion	$inversion_truth>>COMPARISONS/$prefix.summary.txt
		echo translocation	$translocation_truth>>COMPARISONS/$prefix.summary.txt

		mv COMPARISONS/${prefix}_* COMPARISONS/$prefix/
	done
fi

mv COMPARISONS/$species.$SV_type.filtered COMPARISONS/$prefix/
mv COMPARISONS/$prefix.summary.txt COMPARISONS/$prefix/
mv COMPARISONS/$prefix.filtered COMPARISONS/$prefix/