# domain-clustering-pipeline

## Introduction
This pipeline is designed to automate the prediction of protein domains given only the protein's sequence, based on the method described in (Granata _et al._, 2017). The pipeline is written in Snakemake framework (Koster and Rahmann, 2012) to allow easy execution both in compute clusters and in local machines and to ensure reproducibility of the generated data. We also strive to make installation as convenient as possible, by having most of our dependencies available on conda. 

The pipeline consists of three main processes: 
- Generation of multiple sequence alignment using an iterative protein sequence search tool called hhblits   
- Prediction of direct coupling between every pair of protein residues using gplmDCA 
- Generation of protein domains clusters using spectral clustering 

This pipeline can be run end-to-end -- i.e. from FASTA sequence input toa table of every residue position and their corresponding domain membership -- or can be started in from intermediate files. 

This repository also include downstream analyses which aims to inspect the possible relationship between the number and identity of ExAC mutations (Lek _et al._, 2016) and pathogenic variants.  

## Dependencies
Python, Conda, MATLAB

## Setup
Before using the pipeline, you will need to install snakemake on top of the dependencies listed above
```
conda install -c bioconda -c conda-forge snakemake
```

Also, hhblits requires a database to be downloaded for which they recommend [UniClust30](https://uniclust.mmseqs.com/). Download this database into the `db/` directory.    

## Running the pipeline on local machine
On the local machine you could simply run the pipeline by firstly populating the `fasta_list.dat` file with linebreak-delimieted name of your fasta files (without the .fa suffix). You would also need to modify the parameters in `params.json`.

Once those are done, you can run 
```
snakemake -k --use-conda 
```

## Running the pipeline on a compute cluster
For parallelise implementation on a compute cluster, you would need to create a configuration file called `cluster.json`, for example
```
{
  "__default__":
    {
      "partition": "general",
      "mem": "1G",
      "time": "00:30:00"
    }
}          
```
then, simply run
```
snakemake -k -j 100 --use-conda --cluster-config cluster.json --cluster "sbatch -p {cluster.partition} --mem={cluster.mem} -t {cluster.time} -c {threads}"
```
where `-j 100` means that you only allow a maxiumum of 100 jobs to be run simultaneously on the cluster. 

## Results 
The result along with all intermediate and log files can be found in the `results/` directory, where the final result of every residue's domain membership is found in `results/clustering`. 

## Analysis 



## Bibliography 
  Granata, D., Ponzoni, L., Micheletti, C., & Carnevale, V. (2017). Patterns of coevolving amino acids unveil structural and dynamical domains. Proceedings of the National Academy of Sciences, 114(50), E10612-E10621.
  Köster, J., & Rahmann, S. (2012). Snakemake—a scalable bioinformatics workflow engine. Bioinformatics, 28(19), 2520-2522.
  Lek, M., Karczewski, K. J., Minikel, E. V., Samocha, K. E., Banks, E., Fennell, T., ... & Tukiainen, T. (2016). Analysis of protein-coding genetic variation in 60,706 humans. Nature, 536(7616), 285.
