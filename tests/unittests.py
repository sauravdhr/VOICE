import unittest
import main
SEQS_0 = ["AAA"]
SEQS_1 = ["AAA", "AAA"]
SEQS_2 = ["AAA", "AAB"]
SEQS_3 = ["AAA", "AAB", "AAA"]
SEQS_4 = ["AAA", "AAB", "AAA", "AAB"]
MATRIX = [21, 31, 32, 41, 42, 43, 51, 52, 53, 54]
MATRIX_FULL = [[0, 21, 31, 41, 51],
               [21, 0, 32, 42, 52],
               [31, 32, 0, 43, 53],
               [41, 42, 43, 0, 54],
               [51, 52, 53, 54, 0]]


class Test(unittest.TestCase):
    def test(self):
        self.assertEqual(main.hamming_distance("AAA", "AAB"), 1)
        self.assertEqual(main.hamming_distance("AAA", "AAA"), 0)

        self.assertEqual(main.filter_repeated_sequences(None), None)
        self.assertEqual(main.filter_repeated_sequences([]), [])
        self.assertEqual(main.filter_repeated_sequences(SEQS_0), SEQS_0)
        self.assertEqual(main.filter_repeated_sequences(SEQS_1), SEQS_0)
        self.assertEqual(main.filter_repeated_sequences(SEQS_2), SEQS_2)
        self.assertEqual(main.filter_repeated_sequences(SEQS_3), SEQS_2)
        self.assertEqual(main.filter_repeated_sequences(SEQS_4), SEQS_2)


class TriangleMatrixTest(unittest.TestCase):
    def test(self):
        self.assertEqual(main.SymMatrixWithoutDiagonal(MATRIX)[0, 0], 0)
        self.assertEqual(main.SymMatrixWithoutDiagonal(MATRIX)[0, 1], 21)
        self.assertEqual(main.SymMatrixWithoutDiagonal(MATRIX)[1, 0], 21)
        self.assertEqual(main.SymMatrixWithoutDiagonal(MATRIX)[1, 1], 0)
        self.assertEqual(main.SymMatrixWithoutDiagonal(MATRIX)[0, 2], 31)
        self.assertEqual(main.SymMatrixWithoutDiagonal(MATRIX)[2, 0], 31)
        self.assertEqual(main.SymMatrixWithoutDiagonal(MATRIX)[2, 1], 32)
        self.assertEqual(main.SymMatrixWithoutDiagonal(MATRIX)[1, 2], 32)
        self.assertEqual(main.SymMatrixWithoutDiagonal(MATRIX)[2, 2], 0)
        self.assertEqual([[e for e in r] for r in main.SymMatrixWithoutDiagonal(MATRIX)], MATRIX_FULL)
