#!/usr/bin/env python3
from Bio import SeqIO

FASTA = "data/AC122_unique_1a_48.fas"


class SymMatrixWithoutDiagonal(object):
# it is upper-packed matrix (see http://www.ibm.com/support/knowledgecenter/SSFHY8_5.3.0/com.ibm.cluster.essl.v5r3.essl100.doc/am5gr_upsm.htm#am5gr_upsm).
    def __init__(self, vector):
        self.vector = vector

    def __getitem__(self, index):
        if len(index) != 2:
            return self.vector[index]
        if index[0] == index[1]:
            return 0
        if index[0] > index[1]:
            j, i = index[0]-1, index[1]
        else:
            i, j = index[0], index[1]-1
        return self.vector[i+int(j*(j+1)/2)]


def hamming_distance(s1, s2):
    #Return the Hamming distance between equal-length sequences
    if len(s1) != len(s2):
        raise ValueError("Undefined for sequences of unequal length")
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))


def filter_repeated_sequences(sequences):
    if not sequences:
        return sequences
    return list(sequences[i] for i in
                [0] + list(
                    filter(
                        lambda i: not(True in map(lambda s2: sequences[i] == s2, sequences[:i])),
                        list(range(1, len(sequences))))))


def get_sequences_distance_matrix(sequences):
    return list(map(
        lambda i: list(map(
            lambda s2: hamming_distance(sequences[i], s2),
            sequences[:i])),
        list(range(1, len(sequences)))))


class Graph(object):
    def __init__(self, fasta):
        self.sequences = self.parse_fasta(fasta)
        self.vertices = self.infer_vertices(self.sequences)

    @staticmethod
    def parse_fasta(fasta):
        seqs = []
        for seq_record in SeqIO.parse(fasta, "fasta"):
            seqs.append(str(seq_record.seq))
        return seqs

    def infer_vertices(self, sequences):
        filtered_sequences = filter_repeated_sequences(sequences)
        distance_matrix = get_sequences_distance_matrix(filtered_sequences)
        for item in distance_matrix:
            print(', '.join(map(str, item[:])))
        return 0


def main(fasta):
#    graph = Graph(fasta)
    matrix = SymMatrixWithoutDiagonal([21, 31, 32, 41, 42, 43, 51, 52, 53, 54])
    print(matrix[0, 1])
    pass

if __name__ == "__main__":
    main(FASTA)
