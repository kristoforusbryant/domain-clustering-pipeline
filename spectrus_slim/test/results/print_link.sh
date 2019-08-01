#!/bin/bash 
# print the diagonal with the two adjacet diagonals beside it 

count=1
while IFS= read -r line 
do 
	if [ ! $count -eq 1 ]; then echo $line | awk -v count=$count -F, '{ prev = count-1; nex = count+1; print $prev, $count, $nex }';fi 
	count=$(( $count + 1 ))

done < $1

