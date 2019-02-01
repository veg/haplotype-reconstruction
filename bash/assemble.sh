qfilt -Q data/$1.fastq -q 15 -l 50 -P - -R 8 -j >> data/$1_qfilt.fna 2>> data/$1_qfilt.json
bealign -r $2.fas -e 0.5 -m HIV_BETWEEN_F -D data/$1_$2_discards.fna -R data/$1_qfilt.fna data/dataset-$1_reference-$2_aligned.bam
