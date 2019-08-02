import json
with open('sample_list.dat', 'r') as filename:
	samples=[row[:-1] for row in filename]
with open('params.json', 'r') as filename: 
	params=json.loads(filename.read())


#configfile: "config.yaml"
#expand("results/plots/{sample}.SC.clustering.plot", sample=fasta_names["samples"])

rule targets: 
	input:
		expand("results/clustering/{sample}.SC", sample=samples), 
		expand("results/plots/{sample}.SC.clustering.png", sample=samples), 
		expand("results/plots/{sample}.SC.quality.png", sample=samples) 
		
rule generateMSA:
	input:
		fa="fasta_input/{sample}.fa"
	output:
		"results/msa/{sample}.msa.a3m"
	params: 
		cpu_per_tasks=params["hhblits"]["cpu_per_task"],
		n_of_iterations=params["hhblits"]["n_of_iterations"],
		e_value=params["hhblits"]["e_value"],
		database=params["hhblits"]["database"]
	log:
		err="results/logs/generateMSA/{sample}.error.log", 
		summary="results/logs/generateMSA/{sample}.summary.log"
	threads:params["hhblits"]["cpu_per_task"]
	conda: 
		"env/hhsuite.yaml"
	shell: 
		"hhblits -cpu {params.cpu_per_tasks} -n {params.n_of_iterations} -e {params.e_value} -i {input.fa} -oa3m {output} -d {params.database} 1> {log.summary} 2> {log.err}"

rule filterMSA_by_gap: 
	input: 
		msa="results/msa/{sample}.msa.a3m"
	output:
		"results/msa/{sample}.filtered.gap.a3m"
	shell: 
		"cat {input.msa} | ./filter_by_gap.sh > {output}"

rule reformatMSA: 
	input:
		"results/msa/{sample}.filtered.gap.a3m"
	output:
		"results/msa/{sample}.filtered.reformatted.a3m"
	shell:
		"cat {input} | ./reformat_a3m_fa.py |"
		"sed 's/-/./g' > {output}"

rule filterMSA_by_count:
	input:
		"results/msa/{sample}.filtered.reformatted.a3m"
	output:
		"results/msa/{sample}.filtered.count.a3m"
	shell:
		"""
		if [ $(cat {input} | grep -c "^>") -lt 100 ]; then 
			echo "Stopping job due to not enough sequence in the MSA (< 100) "
			sample=$(echo {input} | sed 's/.*\///' | sed 's/.filtered.reformatted.a3m//'); echo $sample
			sed -i '/$sample/d' sample_list.dat && echo "$sample FAILED: TOO FEW MSA SEQUENCES" >> fasta_done.dat
			exit 1 
		else cp {input} {output}
		fi
		"""
"""
rule conserve_scoring:
	input:
		"results/msa/{sample}.filtered.count.a3m"
	output:
		"results/covscore/{sample}.covscore"
	shell: 
	... 
"""

rule generateMSA_stats:
	input: 
		fa="fasta_input/{sample}.fa",
		count="results/msa/{sample}.filtered.gap.a3m"

	output: 
		"results/msa_stats/{sample}.msa.stats"
	shell: 
		"""
		cat {input.fa} | sed -n '1,2p' >> results/msa_stats/{sample}.msa.stats
		echo 'Protein Sequence Length:' >> results/msa_stats/{sample}.msa.stats 
		cat {input.fa} | sed -n '2p' | wc -c >> results/msa_stats/{sample}.msa.stats
		echo 'MSA count' >> results/msa_stats/{sample}.msa.stats
		wc -l {input.count} >> results/msa_stats/{sample}.msa.stats
		"""
rule gplmDCA: 
	input: 
		msa="results/msa/{sample}.filtered.count.a3m"
	output: 
		"results/dca/{sample}.gplmDCA"
	threads:params["gplmDCA"]["nr_of_cores"]
	params: 
		lambda_h=params["gplmDCA"]["lambda_h"],
		lambda_J=params["gplmDCA"]["lambda_J"],
		lambda_chi=params["gplmDCA"]["lambda_chi"],
		reweighting_threshold=params["gplmDCA"]["reweighting_threshold"],
		nr_of_cores=params["gplmDCA"]["nr_of_cores"], 
		M=params["gplmDCA"]["M"]
	log: 
		"results/gplmDCA/{sample}.log"
	shell: 
		"matlab -nodisplay -nosplash -r \"gplmDCA_asymmetric('{input.msa}', '{output}', {params.lambda_h}, {params.lambda_J}, {params.lambda_chi}, {params.reweighting_threshold}, {params.nr_of_cores}, {params.M})\" 2> {log}"

rule spectral_clustering: 
	input: 
		dca="results/dca/{sample}.gplmDCA"
	output:
		clust="results/clustering/{sample}.SC",
		stats="results/clustering_stats/{sample}.SCstats",
		toplot="results/clustering_stats/results_{sample}"
	conda:
		"env/spectrus.yaml"
	shell:
		"""
		sample=$(echo {input.dca} | sed 's/.*\///' | sed 's/.gplmDCA//')
		./cluster_spectrus.sh {input.dca} $sample {output.clust} results/clustering_stats/ > {output.stats} &&
		sed -i '/$sample/d' sample_list.dat && echo "$sample OK" >> fasta_done.dat
		rm -rf results/clustering_stats/results_$sample
		mv -f results_"$sample".temp/ {output.toplot}
		""" 

rule output_graph: 
	input: 
		SC="results/clustering_stats/results_{sample}"
	output: 
		clust="results/plots/{sample}.SC.clustering.png", 
		qual="results/plots/{sample}.SC.quality.png"
	conda:
		"env/Rplot.yaml"
	shell: """
		mkdir -p results/plots;
		Rscript {input.SC} {output.clust} {output.qual}
		"""
