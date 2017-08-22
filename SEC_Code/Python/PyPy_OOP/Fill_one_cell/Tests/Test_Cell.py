import unittest
from parameterized import parameterized
from mock import patch, Mock

from Cell import Cell
from Population import Population


class TestCell(unittest.TestCase):

    # [[1, 1, 0],
    #  [0, 0, 1]]

    def setUp(self):
        self.num_rows = 2
        self.num_cols = 3

        self.Matrix = Cell.initialize_matrix(self.num_cols, self.num_rows)

        # N = 10, pA = 0.5, pB = 0.5
        mock_pop_1 = Mock(spec=Population, size=10, locus_A=(['A'] * 5) + (['a'] * 5), locus_B=(['B'] * 5) + (['b'] * 5))
        # N = 17, pA = 0.7059, pB = 0.1176
        mock_pop_2 = Mock(spec=Population, size=17, locus_A=(['A'] * 12) + (['a'] * 5), locus_B=(['B'] * 2) + (['b'] * 15))
        # N = 32, pA = 0, pB =0.6875
        mock_pop_3 = Mock(spec=Population, size=32, locus_A=(['A'] * 0) + (['a'] * 32), locus_B=(['B'] * 22) + (['b'] * 10))

        self.Matrix[0][0].pop = True
        self.Matrix[0][1].pop = True
        self.Matrix[1][2].pop = True

        self.Matrix[0][0].population = mock_pop_1
        self.Matrix[0][1].population = mock_pop_2
        self.Matrix[1][2].population = mock_pop_3

        # print Population.allele_freq(self.Matrix[1][2].population.locus_B)

    def tearDown(self):
        del self.Matrix

    def test_initialize_matrix(self):

        num_rows = 2
        num_cols = 5
        Matrix = Cell.initialize_matrix(num_cols, num_rows)

        for i in range(num_rows):
            for j in range(num_cols):
                assert isinstance(Matrix[i][j], Cell)

    def test_init(self):

        for i in range(self.num_rows):
            for j in range(self.num_cols):
                self.assertEqual(self.Matrix[i][j].i, i)
                self.assertEqual(self.Matrix[i][j].j, j)

        self.assertTrue(self.Matrix[0][0].pop)
        self.assertTrue(self.Matrix[0][1].pop)
        self.assertTrue(self.Matrix[1][2].pop)

        self.assertFalse(self.Matrix[0][2].pop)
        self.assertFalse(self.Matrix[1][0].pop)
        self.assertFalse(self.Matrix[1][1].pop)

        assert isinstance(self.Matrix[0][0].population, Population)
        self.assertEqual(self.Matrix[0][0].population.size, 10)

if __name__ == '__main__':
    unittest.main()


# def test_init():
#     assert Matrix[0][0].i == 0


# def test_prob_create():
#     assert Matrix[0][0].prob_create(1000, 0) == 0

#     assert Matrix[0][0].prob_create(1000, 0.5) == 0.005
#     assert Matrix[0][0].prob_create(10, 0.5) == 0.5

#     assert Matrix[0][0].prob_create(1000, 1) == 0.01


# def test_create_population():
#     Matrix[0][0].create_population(rows, cols, Matrix, 1.0, 1.0, 10, 1.5)
#     assert Matrix[0][1].pop
#     assert Matrix[0][1].population.size == 10
