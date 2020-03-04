function mytest {
    "$@"
    local status=$?
    if [ $status -ne 0 ]; then
        echo "error with $@" >> errors.txt
    fi
    return $status
}
echo "Pipeline issues:" > errors.txt

mytest snakemake -n output/sim-divergedPair_ar-10_seed-1/trimmomatic/bowtie2/pol/acme/sr-0/sorted.fasta
mytest snakemake -n output/sim-divergedPair_ar-10_seed-1/trimmomatic/bowtie2/pol/sorted.bam
mytest snakemake -n output/sim-divergedPair_ar-10_seed-1/trimmomatic/bowtie2/pol/acme/superreads.json
mytest snakemake -n output/truth/sim-divergedPair/pol_gene.json
mytest snakemake -nf output/sim-divergedPair_ar-10_seed-1/reads.fastq

errors=$(grep -c error errors.txt)
if [ $errors -ne 0 ]; then
  echo "ERRORS IN PIPELINE LOGIC!"
  exit 1
fi
exit 0
