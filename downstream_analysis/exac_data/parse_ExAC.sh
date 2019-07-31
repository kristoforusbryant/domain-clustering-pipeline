#!/bin/bash
# Filter raw ExAC vcf and remove unnecessary columns
# NOTE: !!! this may take a while, we recommend running it in parallel !!!

if [ -z $1 ]; then echo "Error: no input"; exit; fi 
if [ -z $2 ]; then echo "Error: no output"; exit; fi

input=$1
name=$(echo $input | sed 's/.*\///' | sed 's/.vcf//')
output=$2 

vcftools --gzvcf $input --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --chr X --chr Y --remove-filtered-all --min-meanDP 10.0 --recode --recode-INFO CSQ --recode-INFO AF --out ${name}.filtered &&

# Parse CSQ 
while read line 
do
	# keep all headers from original vcf file 
	regex="^\#.*"
	if [[ "$line" =~ $regex ]]; then echo $line; continue; fi 	
	
	CSQ=$(echo $line | grep -o -e ';CSQ=[^;]*')
	for i in $(echo $CSQ | sed 's/^;CSQ=//' | sed 's/,/ /g')
	do
		echo $(echo $line | sed "s/;CSQ=[^;]*//")\;CSQ=$i # keep only CSQ INFO 
	done

# if mutation span a range of residues, summarise it as mutation on the first residue position 
done < ${name}.filtered.recode.vcf | grep protein_coding |
	awk '{sub(/AF=/, "", $8); sub(/;CSQ=/,"|", $8); print $8}' |
	awk -F"|" 'OFS=","{ 
		$1 = gensub(/,.*$/,"","g",$1);
		residue_position = gensub(/-.*$/, "","g", $16); 
		if(length(residue_position) > 0) {print $30"-"residue_position" "$6,$1,$3,$34,$35}  
	}' > $output &&

rm ${name}.filtered
