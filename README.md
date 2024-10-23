# CoverageStats
Calculate coverage statistics after mapping short-read sequenced isolates to a reference sequence

# Installation of package
git clone https://github.com/kylegontjes/CoverageStats.git

# Before running the pipeline
## 1. Set up the reference genome  
### Create index of reference genome
module load Bioinformatics 

module load bwa 

bwa index [reference genome] 
### Create fai index file
module load samtools 

samtools faidx [reference genome] 

### Create sequence dictionary
module load gatk 

gatk CreateSequenceDictionary -R [reference genome] 

## 2. Generate isolate list
path="/nfs/esnitkin/Project_Penn_KPC/Sequence_data/fastq/Penn/SRA_submission/"

sample_id="sample_id" 

sample_names=$(ls -1 $path | grep _R1 | cut -d. -f1  | sed 's\_R1\\' | sed 's\_R2\\' | sort | uniq)

echo -e\n $sample_id $sample_names | tr ' ' '\n' > config/sample.tsv

# Running snakemake workflow
module load singularity

module load snakemake

snakemake -s CoverageStats.smk --use-conda --use-singularity -j 999 --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes} -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem} --output=slurm_out/slurm-%j.out" --conda-frontend conda --cluster-config config/cluster.json --configfile config/config.yaml --latency-wait 30 --keep-going 

# To run many isolates at the same time (and possibly close the computer, try running this command)
srun --account=esnitkin1 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=5GB --cpus-per-task=1 --time=12:00:00 --pty /bin/bash
