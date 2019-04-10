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

information['numberOfInsertions'] = len(information['insertions'])
information['numberOfDeletions'] = len(information['deletions'])
information['totalNumberOfReads'] = total_reads
information['numberOfReadsWithInsertions'] = reads_with_insertion
information['numberOfReadsWithDeletions'] = reads_with_deletion
print(json.dumps(information, indent=2))
