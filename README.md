# CoverageStats
Calculate coverage statistics after mapping short-read sequenced isolates to a reference sequence

# Generate isolate list
path="/nfs/esnitkin/Project_Penn_KPC/Sequence_data/fastq/Penn/SRA_submission/"

sample_id="sample_id" 

sample_names=$(ls -1 $path | grep _R1 | cut -d. -f1  | sed 's\_R1\\' | sed 's\_R2\\' | sort | uniq | head -n3)

echo -e\n $sample_id $sample_names | tr ' ' '\n' > config/sample.tsv

# Run command
snakemake -s CoverageStats.smk --use-conda --use-singularity -j 999 --cluster "sbatch -A {cluster.account} -p {cluster.partition} -N {cluster.nodes} -t {cluster.walltime} -c {cluster.procs} --mem-per-cpu {cluster.pmem}" --conda-frontend conda --cluster-config config/cluster.json --configfile config/config.yaml --latency-wait 60