import numpy as np


class AlignedSegment:

    def __init__(self, pysam_aligned_segment):
        self.pysam_aligned_segment = pysam_aligned_segment
        self.current_cigar_tuple_index = 0
        self.position_in_current_cigar_tuple = 0
        self.position_along_reference = pysam_aligned_segment.reference_start
        self.active = None

    @property
    def cigartuples(self):
        return self.pysam_aligned_segment.cigartuples

    @property
    def current_cigar_tuple(self):
        return self.cigartuples[self.current_cigar_tuple_index]

    def advance_cigar(self):
        current_action = self.current_cigar_tuple[0]
        current_stride = self.current_cigar_tuple[1]
        if self.position_in_current_cigar_tuple == current_stride:
            self.current_cigar_tuple_index += 1
            self.position_in_current_cigar_tuple = 0
        else:
            self.position_in_current_cigar_tuple += 1
        if current_action == 0 or current_action == 2:
            self.position_along_reference += 1

    def initiate_left_extended_segment(self, window_start):
        while self.position_along_reference < window_start:
            self.advance_cigar()
        self.active = True

    def initiate_exactly_placed_segment(self):
        self.active = True

    def initiate_right_extended_segment(self):
        self.active = False

    def initiate_conversion(self, window_start):
        if self.position_along_reference < window_start:
            self.initiate_left_extended_segment(window_start)
        elif self.position_along_reference == window_start:
            self.initiate_exactly_placed_segment()
        else:
            self.initiate_right_extended_segment()


class SAMFASTAConverter:

    def __init__(self):
        self.nucleotides = list('ACGTRYSWKMBDHVN-~')
        self.n_nucleotides = len(self.nucleotides)
        self.embedding = np.array([
            [1,   0,   0,    0,   0, 0], # A
            [0,   1,   0,    0,   0, 0], # C
            [0,   0,   1,    0,   0, 0], # G
            [0,   0,   0,    1,   0, 0], # T
            [.5,  0,  .5,    0,   0, 0], # R
            [0,   .5,  0,    .5,  0, 0], # Y
            [0,   .5,  .5,   0,   0, 0], # S
            [.5,  0,   0,    .5,  0, 0], # W
            [0,   0,   .5,   .5,  0, 0], # K
            [.5,  .5,  0,    0,   0, 0], # M
            [0,   1/3, 1/3,  1/3, 0, 0], # B
            [1/3, 0,   1/3,  1/3, 0, 0], # D
            [1/3, 1/3, 0,    0/3, 0, 0], # H
            [1/3, 1/3, 1/3,  0,   0, 0], # V
            [.25, .25, .25, .25,  0, 0], # N
            [0,   0,   0,    0,   1, 0], # -
            [0,   0,   0,    0,   0, 1], # ~
        ])
        self.nuc2ind = {
            nuc: i for i, nuc in enumerate(self.nucleotides)
        }

    def get_numeric_representation(self, record):
        return np.array(
            [self.nuc2ind[char] for char in record.seq],
            dtype=np.int
        )

    def single_aligned_segment_to_fasta(self, aligned_segment):
        read_sequence = aligned_segment.query_alignment_sequence
        position = 0
        alignment = ''
        cigar_tuples = aligned_segment.cigartuples
        number_of_cigar_tuples = len(cigar_tuples)
        for i, cigar_tuple in enumerate(cigar_tuples):
            action = cigar_tuple[0]
            stride = cigar_tuple[1]
            match = action == 0
            insertion = action == 1
            deletion = action == 2
            if match or insertion:
                alignment += read_sequence[position: position + stride]
                position += stride
            elif deletion:
                if len(alignment) > 0 and i < number_of_cigar_tuples:
                    alignment += stride * '-'
        full_seq_np = np.array(
            list(alignment),
            dtype='<U1'
        )
        start = np.argmax(full_seq_np != '-')
        stop = len(full_seq_np) - np.argmax(full_seq_np[::-1] != '-')
        return full_seq_np[start:stop]

    def multiple_aligned_segments_in_window_to_fasta(self, aligned_segments,
            reference_length, window_start, window_end, outer_gap_char='~'):
        number_of_rows = len(aligned_segments)
        number_of_columns = 2*reference_length
        fasta = np.array(
            [number_of_columns*[outer_gap_char] for _ in range(number_of_rows)],
            dtype='<U1'
        )

        current_cigar_tuple = np.zeros(number_of_rows, dtype=np.int)
        place_in_cigar_tuple = np.zeros(number_of_rows, dtype=np.int)
        place_in_pairwise_aligned_segment = np.zeros(number_of_rows, dtype=np.int)
        for i, segment in enumerate(aligned_segments):
            pass 


    def embed_numeric_fasta(numeric_fasta):
        number_of_rows = numeric_fasta.shape[0]
        number_of_molecules = numeric_fasta.shape[1]
        size_of_embedding = self.embedding.shape[1]
        number_of_columns = number_of_molecules*size_of_embedding
        embedded_fasta = np.reshape(
            embedding[numeric_fasta, :],
            newshape=(number_of_rows, number_of_columns) 
        )

