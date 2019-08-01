configfile: "fasta_names.yaml"

rule targets: 
	input:
		expand("results/plots/{sample}.SC.clustering.plot", sample=fasta_names["samples"])

rule generateMSA:
	input:
		fa="fasta_input/{sample}.fa"
		cpu_per_tasks=1
		n_of_iterations=5
		e_value=0.001
		database="db/uniclust30_2017_10/uniclust30_2017_10"
	output:
		"results/msa/{sample}.msa.a3m"
	conda: 
		"envs/hhsuite.yaml"
	shell: 
		""" hhblits 
				-cpu cpu_per_tasks
				-n n_of_iterations 
            	-e e_value
            	-i {input.fa}  
            	-oa3m {sample}.msa.oa3m 
            	-d database """

rule filterMSA_by_gap: 
	input: 
		msa="results/msa/{sample}.msa.a3m"
	output:
		"results/msa/{sample}.filtered.gap.a3m"
	shell: 
		"cat {input.msa} | ./filterMSA_by_gap.sh >" 
		"results/msa/{sample}.filtered.gap.a3m"

rule reformatMSA: 
	input:
		"results/msa/{sample}.filtered.gap.a3m"
	output:
		"results/msa/{sample}.filtered.reformated.a3m"
	shell:
		"cat {input} | ./reformat_a3m_fa.py |"
		"sed 's/-/./g' > results/msa/{sample}.filtered.reformated.a3m"

rule filterMSA_by_count:
	input:
		"results/msa/{sample}.filtered.reformatted.a3m"
	output:
		"results/msa/{sample}.filtered.count.a3m"
	shell:
		"if [ $(cat {input} | grep -c "^>") -lt 100 ]; then"
		"exit 1; else cp {input} results/msa/{sample}.filtered.reformated.a3m;fi"

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
		"cat {input.fa} | sed -n '1,2p' >> results/msa_stats/{sample}.msa.stats;"
		"echo 'Protein Sequence Length:' "
		"cat {input.fa} | sed -n '2p' | wc -c >> results/msa_stats/{sample}.msa.stats;"
		"echo 'MSA count' >> results/msa_stats/{sample}.msa.stats"
		"wc -l {input.count} >> results/msa_stats/{sample}.msa.stats"

rule gplmDCA: 
	input: 
		msa="results/msa/{sample}.filtered.count.a3m",
		lambda_h="0.01",
		lambda_J="0.01",
		lambda_chi="0.001",
		reweighting_threshold="0.1",
		nr_of_cores="10", # may need changing depending on number of CPUs available
		M="-1"
	output: 
		"results/dca/{sample}.gplmDCA"
	shell: 
		"matlab -nodisplay -nosplash -r "
		"\"gplmDCA_asymmetric('{input.msa}',"
		"'{output}', {input.lambda_h}, {input.lambda_J},"
		" {input.lambda_chi}, {input.reweighting_threshold},"
		" {input.nr_of_cores}, {n.M})\""

rule spectral_clustering: 
	input: 
		"results/dca/{sample}.gplmDCA"
	output:
		clust="results/clustering/{sample}.SC",
		stats="results/clustering_stats/{sample}.SCstats"
	conda: 
		"env/slepc.yaml"
	shell:
		"cp -r ./spectrus_slim/ ./spectrus_slim_{sample}/;"
		"./cluster_spectrus.sh {input} {sample} {output.clust} 2> {output.stats}" 
			# this command (and cluster_spectrus.sh) is a little hack-y

rule output_graph: 
	input: 
		SC="results/clustering/{sample}.SC"
	output: 
		clust="results/plots/{sample}.SC.clustering.png", 
		qual="results/plots/{sample}.SC.quality.png"
	conda:
		"env/Rplot.yaml"
	shell: 
		"Rscript {input.SC} {output.clust} {output.qual}"

