python extract_sequences.py
art_illumina -ss HS25 -i data/related_1.fasta -l 120 -s 50 -c 150 -o data/related_1_reads
art_illumina -ss HS25 -i data/related_2.fasta -l 120 -s 50 -c 150 -o data/related_2_reads
art_illumina -ss HS25 -i data/diverged_1.fasta -l 120 -s 50 -c 150 -o data/diverged_1_reads
art_illumina -ss HS25 -i data/diverged_2.fasta -l 120 -s 50 -c 150 -o data/diverged_2_reads
cat data/related_1_reads.fq data/related_2_reads.fq > data/related.fastq
cat data/diverged_1_reads.fq data/diverged_2_reads.fq > data/diverged.fastq
