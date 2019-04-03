import json
import os
import subprocess


datasets = []
for x in os.walk('input/compartmentalization/'):
    if len(x[2]) > 0:
        directory = x[0]
        reads_filename = x[2][0] if x[2][0].split('.')[-1] == 'fna' else x[2][1]
        qual_filename = x[2][0] if x[2][0].split('.')[-1] == 'qual' else x[2][1]
        reads = directory + "/" + reads_filename
        quality = directory + "/" + qual_filename
        new_reads = directory + "/reads.fasta"
        new_quality = directory + "/scores.qual"
        subprocess.run(["mv", reads, new_reads])
        subprocess.run(["mv", quality, new_quality])
        datasets.append('-'.join(directory.replace('CSF-PELLET','CSF_PELLET').split('/')[2:]))

with open('compartmentalization.json', 'w') as json_file:
    json.dump(datasets, json_file, indent=2)

