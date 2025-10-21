# CoverageStats
Calculate coverage statistics after mapping short-read sequenced isolates to a reference sequence

# Installation of package
git clone https://github.com/kylegontjes/CoverageStats.git

# Before running the pipeline
## 1. Set up the reference genome 
### Approach #1: Helper script 
bash index_reference_genome.sh [reference_genome.fasta]

### Approach #2: Manual curation
### A. Create index of reference genome
module load Bioinformatics 

module load bwa 

bwa index [reference_genome.fasta]

### B. Create fai index file
module load samtools 

samtools faidx [reference_genome.fasta]

### C. Create sequence dictionary
module load gatk 

gatk CreateSequenceDictionary -R [reference_genome.fasta]

## 2. Generate isolate list
path="/nfs/esnitkin/Project_Penn_KPC/Sequence_data/fastq/Penn/SRA_submission/"

sample_id="sample_id" 

sample_names=$(ls -1 $path | grep _R1 | cut -d. -f1 | sed 's/_R1//')

echo -e\n $sample_id $sample_names | tr ' ' '\n' > config/sample.tsv

## 3. Edit the config/cluster.json file with your account, email address (email), and desired walltime (walltime)

## 4. Update the config/config.yaml or config/config_pretrimmed.yaml with sample list (samples), results folder name (prefix), input reads' path (input_reads), and the reference genome's path (reference_genome)

# Running snakemake workflow
module load singularity
module load snakemake

# Dry run
snakemake -s CoverageStats.smk --dryrun -p

# Preferred method: Running sbat script 
## On non-trimmed files
sbatch CoverageStats.sbat 

## On pre-trimmed files
sbatch CoverageStats_pretrimmed.sbat


# Alternative method: Run command in terminal as interactive job
## On non-trimmed files
snakemake -s CoverageStats.smk --use-singularity -j 999 --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes} -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem} --output=slurm_out/slurm-%j.out" --cluster-config config/cluster.json --configfile config/config.yaml --latency-wait 30 --keep-going 
## On pre-trimmed files
snakemake -s CoverageStats_pretrimmed.smk --use-singularity -j 999 --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes} -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem} --output=slurm_out/slurm-%j.out" --cluster-config config/cluster.json --configfile config/config_pretrimmed.yaml --latency-wait 30 --keep-going 