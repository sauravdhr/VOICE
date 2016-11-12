import unittest
import network_creator
SEQS_0 = ["AAA"]
SEQS_1 = ["AAA", "AAA"]
SEQS_2 = ["AAA", "AAB"]
SEQS_3 = ["AAA", "AAB", "AAA"]
SEQS_4 = ["AAA", "AAB", "AAA", "AAB"]
SEQS_5 = ["AA", "AB", "AC"]
SEQS_6 = ["AAA", "ABA", "ABB"]
SEQS_7 = ["AAA", "AAB", "AAA", "AAB", "AAC"]
SEQS_8 = ["BAB", "AAB", "BAC", "ACB", "AAC"]

MATRIX = [21, 31, 32, 41, 42, 43, 51, 52, 53, 54]
MATRIX_FULL = [[0, 21, 31, 41, 51],
               [21, 0, 32, 42, 52],
               [31, 32, 0, 43, 53],
               [41, 42, 43, 0, 54],
               [51, 52, 53, 54, 0]]


class Test(unittest.TestCase):
    def test(self):
        self.assertEqual(network_creator.hamming_distance("AAA", "AAB"), 1)
        self.assertEqual(network_creator.hamming_distance("AAA", "AAA"), 0)

        self.assertEqual(network_creator.filter_repeated_sequences(None), None)
        self.assertEqual(network_creator.filter_repeated_sequences([]), [])
        self.assertEqual(network_creator.filter_repeated_sequences(SEQS_0), SEQS_0)
        self.assertEqual(network_creator.filter_repeated_sequences(SEQS_1), SEQS_0)
        self.assertEqual(network_creator.filter_repeated_sequences(SEQS_2), SEQS_2)
        self.assertEqual(network_creator.filter_repeated_sequences(SEQS_3), SEQS_2)
        self.assertEqual(network_creator.filter_repeated_sequences(SEQS_4), SEQS_2)

        self.assertEqual(network_creator.get_sequence_sets_difference(SEQS_0, SEQS_0), [])
        self.assertEqual(network_creator.get_sequence_sets_difference(SEQS_8, SEQS_8), [])
        self.assertEqual(network_creator.get_sequence_sets_difference(SEQS_7, SEQS_8), ["AAA", "AAA"])
        self.assertEqual(network_creator.get_sequence_sets_difference(SEQS_8, SEQS_7), ["BAB", "BAC", "ACB"])


class TriangleMatrixTest(unittest.TestCase):
    def test(self):
        self.assertEqual(network_creator.SymMatrixWithoutDiagonal(MATRIX)[0, 0], 0)
        self.assertEqual(network_creator.SymMatrixWithoutDiagonal(MATRIX)[0, 1], 21)
        self.assertEqual(network_creator.SymMatrixWithoutDiagonal(MATRIX)[1, 0], 21)
        self.assertEqual(network_creator.SymMatrixWithoutDiagonal(MATRIX)[1, 1], 0)
        self.assertEqual(network_creator.SymMatrixWithoutDiagonal(MATRIX)[0, 2], 31)
        self.assertEqual(network_creator.SymMatrixWithoutDiagonal(MATRIX)[2, 0], 31)
        self.assertEqual(network_creator.SymMatrixWithoutDiagonal(MATRIX)[2, 1], 32)
        self.assertEqual(network_creator.SymMatrixWithoutDiagonal(MATRIX)[1, 2], 32)
        self.assertEqual(network_creator.SymMatrixWithoutDiagonal(MATRIX)[2, 2], 0)
        self.assertEqual([[e for e in r] for r in network_creator.SymMatrixWithoutDiagonal(MATRIX)], MATRIX_FULL)


class MediansCreationTest(unittest.TestCase):
    def test(self):
        self.assertEqual(network_creator.find_all_medians(SEQS_5), [])
        self.assertEqual(network_creator.find_all_medians(SEQS_6), ["ABA"])
        self.assertEqual(network_creator.find_all_medians(SEQS_4), ["AAA", "AAB"])
        self.assertEqual(network_creator.find_all_medians(SEQS_7), ["AAA", "AAB"])
        self.assertEqual(network_creator.find_all_medians(SEQS_8), ["BAB", "AAB", "BAC", "AAC"])
