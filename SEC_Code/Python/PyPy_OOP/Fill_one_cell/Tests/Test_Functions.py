import unittest
from parameterized import parameterized
from mock import patch

import Functions


class TestPopulation(unittest.TestCase):

    @patch('Cell.Cell')
    def test_matrix_full(self, mock_cell):

        Matrix = [[0, 0], [0, 0]]

        for i in range(len(Matrix)):
            for j in range(len(Matrix[0])):
                Matrix[i][j] = mock_cell
                Matrix[i][j].pop = False

        self.assertEqual(Functions.matrix_full(Matrix), 0)

        for i in range(len(Matrix)):
            for j in range(len(Matrix[0])):
                Matrix[i][j] = mock_cell
                Matrix[i][j].pop = True

        self.assertEqual(Functions.matrix_full(Matrix), 1)

        Matrix[0][0].pop = True
        Matrix[0][1].pop = True
        Matrix[1][0].pop = True
        Matrix[1][1].pop = False

        self.assertEqual(Functions.matrix_full(Matrix), 0)

    @parameterized.expand([([0.0, 1.0], [0.0, 1.0]), ([0.5, 0.5], [0.5, 1.0]), ([0.2, 0.8], [0.2, 1.0])])
    def test_cdf(self, weights, expected):

        self.assertEqual(Functions.cdf(weights), expected)

    @patch('random.random')
    @patch('Functions.cdf')
    def test_choice(self, mock_cdf, mock_random):

        possibilities = 'Aa'

        mock_cdf.return_value = [0.0, 1.0]
        mock_random.return_value = 0.5
        weights = [0.0, 1.0]
        expected = 'a'

        self.assertEqual(Functions.choice(possibilities, weights), expected)

        mock_cdf.return_value = [0.5, 1.0]
        mock_random.return_value = 0.4
        weights = [0.5, 0.5]
        expected = 'A'

        self.assertEqual(Functions.choice(possibilities, weights), expected)

        mock_cdf.return_value = [0.5, 1.0]
        mock_random.return_value = 0.6
        weights = [0.5, 0.5]
        expected = 'a'

        self.assertEqual(Functions.choice(possibilities, weights), expected)

        mock_cdf.return_value = [0.2, 1.0]
        mock_random.return_value = 0.1
        weights = [0.2, 0.8]
        expected = 'A'

        self.assertEqual(Functions.choice(possibilities, weights), expected)

        mock_cdf.return_value = [0.2, 1.0]
        mock_random.return_value = 0.3
        weights = [0.2, 0.8]
        expected = 'a'

        self.assertEqual(Functions.choice(possibilities, weights), expected)

        mock_cdf.return_value = [0.2, 1.0]
        mock_random.return_value = 0.201
        weights = [0.2, 0.8]
        expected = 'A'

        self.assertNotEqual(Functions.choice(possibilities, weights), expected)
