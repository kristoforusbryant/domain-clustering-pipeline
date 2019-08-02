#!/bin/bash 
# Cluster one DCA using Spectral Clustering 

# make spectrus if spectrus don't already exist
if ! [ -e spectrus ]; then
	echo "spectrus not found, initializing spectrus"
	cd third-party-code/spectrus 
	if [ -e spectrus ]; then rm spectrus;fi 
	if make; then
		cp ./spectrus ./INPUT_PARAMS.DAT ../..
		cd ../..
	else
		echo "Unable to make spectrus, check if the packages in env/spectrus.yaml are properly installed; exiting.." 
		exit 1
	fi
fi  

DCA_path=$1
sample=$2
output=$3

# extract only the weigths
cat $DCA_path | awk -F, '{print $3}' > ${sample}.temp
# get the number of residues 
n_residues=$(cat $DCA_path | tail -1 | awk -F, '{ print $2}')
echo $n_residues
echo ${sample}.temp
pwd

# run spectrus 
./spectrus ${sample}.temp $n_residues 
 
# move results we are interested in to clust_bins
best_cluster=$(cat results_${sample}.temp/quality_score.dat | sort -k2 -g -r | head -1 | awk '{print $1}')
mv results_${sample}.temp/final_labels_kmed-${best_cluster}.dat ${output}

# move results directory to results directory 
rm -rf results/clustering_stats/results_${sample}
mv -f results_${sample}.temp/ results/clustering_stats/results_${sample}
rm ${sample}.temp
