import pysam

from .error_correction import ErrorCorrection


def error_correction_io(input_bam, output_bam):
    alignment = pysam.AlignmentFile(input_bam, 'rb')
    error_correction = ErrorCorrection(alignment)
    error_correction.write_corrected_reads(output_bam)

