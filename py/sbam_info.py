import argparse
import json

import pysam


parser = argparse.ArgumentParser(
    description='Count insertions and deletions that occur in a SAM/BAM file.'
)

parser.add_argument(
    '-i', '--input',
    help='input SAM/BAM file'
)

args = parser.parse_args()
input_filename = args.input
open_mode = 'r' if input_filename.split('.')[-1] == 'sam' else 'rb'
pysam_alignment = pysam.AlignmentFile(input_filename, open_mode)

information = {
    'insertions': [],
    'deletions': []
}

total_reads = 0
reads_with_insertion = 0
reads_with_deletion = 0
reads_with_nonzero_query_start = 0
max_read_length = 0
min_read_length = 100000
for read in pysam_alignment.fetch():
    total_reads += 1
    saw_insertion = False
    saw_deletion = False
    for cigar_tuple in read.cigartuples:
        action = cigar_tuple[0]
        if action == 1:
            saw_insertion = True
            information['insertions'].append(read.query_name)
        if action == 2:
            saw_deletion = True
            information['deletions'].append(read.query_name)
    if saw_insertion:
        reads_with_insertion += 1
    if saw_deletion:
        reads_with_deletion += 1
    if read.query_alignment_start != 0:
        reads_with_nonzero_query_start += 1
    max_read_length = max(max_read_length, read.infer_query_length())
    min_read_length = min(min_read_length, read.infer_query_length())

information['numberOfInsertions'] = len(information['insertions'])
information['numberOfDeletions'] = len(information['deletions'])
information['totalNumberOfReads'] = total_reads
information['numberOfReadsWithInsertions'] = reads_with_insertion
information['numberOfReadsWithDeletions'] = reads_with_deletion
information['numberOfReadsWithNonzeroQueryStart'] = reads_with_nonzero_query_start
information['minimumReadLength'] = min_read_length
information['maximumReadLength'] = max_read_length
print(json.dumps(information, indent=2))
