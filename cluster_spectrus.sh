#!/bin/bash 
# Cluster one DCA using Spectral Clustering 

DCA_path=$1
sample=$2
output=$3

cd spectrus_slim_$sample
# extract only the weigths
cat $DCA_path | awk -F, '{print $3}' > ${sample}.temp
# get the number of residues 
n_residues=$(cat $DCA_path | tail -1 | awk -F, '{ print $2}')
echo $n_residues
# run spectrus 
./spectrus ${sample}.temp $n_residues 
 
# move results we are interested in to clust_bins
best_cluster=$(cat results_${sample}.temp/quality_score.dat | sort -k2 -g -r | head -1 | awk '{print $1}')
mv results_${sample}.temp/final_labels_kmed-${best_cluster}.dat ${output}

# remove the results directory 
rm -r results_${sample}.temp
rm ${sample}.temp

cd .. && rm -r spectrus_slim_$sample
