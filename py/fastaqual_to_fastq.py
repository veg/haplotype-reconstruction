# Converts a FASTA and QUAL file to FASTQ
# using Biopython (Bio). Runs on Python 2.

# The libraries we need #
import sys, os, Bio
from Bio import SeqIO
from Bio.SeqIO.QualityIO import PairedFastaQualIterator
# Get the shell arguments #
fa_path = sys.argv[1]
qa_path = sys.argv[2]
fq_path = sys.argv[3]
# Check that the paths are valid #
if not os.path.exists(fa_path): raise Exception("No file at %s." % fa_path)
if not os.path.exists(qa_path): raise Exception("No file at %s." % qa_path)
# Do it #
with open(fq_path, "w") as handle:
    records = PairedFastaQualIterator(open(fa_path), open(qa_path))
    count = SeqIO.write(records, handle, "fastq")
# Report success #
print "Converted %i records" % count
