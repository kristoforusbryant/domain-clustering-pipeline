#!/bin/bash 
# Remove sequence alignments which have gaps above a certain cutoff 

cutoff=0.3

cat | paste -d"@" - - | awk -F "@" -v cutoff=$cutoff 'BEGIN {OFS="@"}{
	len= length ($2)
	gap_len = gsub(/-/, "-", $2)
	gap_perc = gap_len/len
	if (gap_perc < cutoff) print $0;
}' | sed 's/@/\n/' 
