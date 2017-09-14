import unittest
from parameterized import parameterized
from mock import patch

from simulations import functions


class TestFunctions(unittest.TestCase):
    """Unit tests for Functions.py"""

    @patch('simulations.cell.Cell')
    def test_matrix_full(self, mock_cell):
        """Tests matrix_full function

        Args:
            mock_cell (:obj:'class instance'): Mick instance of Cell class from Cell.py
        """
        Matrix = [[0, 0], [0, 0]]

        # Correctly returns 0 if Matrix is completely empty
        for i in range(len(Matrix)):
            for j in range(len(Matrix[0])):
                Matrix[i][j] = mock_cell
                Matrix[i][j].pop = False

        self.assertEqual(functions.matrix_full(Matrix), 0)

        # Correctly returns 1 if Matrix is completely full
        for i in range(len(Matrix)):
            for j in range(len(Matrix[0])):
                Matrix[i][j] = mock_cell
                Matrix[i][j].pop = True

        self.assertEqual(functions.matrix_full(Matrix), 1)

        # Correctly returns 0 if Matrix is only partially full
        Matrix[0][0].pop = True
        Matrix[0][1].pop = True
        Matrix[1][0].pop = True
        Matrix[1][1].pop = False

        self.assertEqual(functions.matrix_full(Matrix), 0)

    @parameterized.expand([([0.0, 1.0], [0.0, 1.0]), ([0.5, 0.5], [0.5, 1.0]), ([0.2, 0.8], [0.2, 1.0])])
    def test_cdf(self, weights, expected):
        """Tests cdf function

        Args:
            weights (:obj:'list' of :obj:'float'): list containing weighting factors. Must sum to 1.
            expected: Expected output from input list of weights.
        """

        # Returns correct cdf for various input lists specified in parameterized.expand.
        self.assertEqual(functions.cdf(weights), expected)

    @patch('random.random')
    @patch('simulations.functions.cdf')
    def test_choice(self, mock_cdf, mock_random):
        """Tests choice function

        Args:
            mock_cdf (:obj:'mocked function'): Mocked version of cdf function from Functions.py
            mock_random (:obj:'mocked method'): Mocked random method from random module
        """
        possibilities = 'Aa'

        # Correctly returns 'a' when probability of sampling 'a' is 1.0
        mock_cdf.return_value = [0.0, 1.0]
        mock_random.return_value = 0.5
        weights = [0.0, 1.0]
        expected = 'a'

        self.assertEqual(functions.choice(possibilities, weights), expected)

        # Correctly returns 'A' when probability of sampling 'A' is 0.5 and random value is less than 0.5
        mock_cdf.return_value = [0.5, 1.0]
        mock_random.return_value = 0.4
        weights = [0.5, 0.5]
        expected = 'A'

        self.assertEqual(functions.choice(possibilities, weights), expected)

        # Correctly returns 'a' when probability of sampling 'A' is 0.5 and random value is greater than 0.5
        mock_cdf.return_value = [0.5, 1.0]
        mock_random.return_value = 0.6
        weights = [0.5, 0.5]
        expected = 'a'

        self.assertEqual(functions.choice(possibilities, weights), expected)

        # Correctly returns 'A' when probability of sampling 'A' is 0.8 and random value is less than 0.2
        mock_cdf.return_value = [0.2, 1.0]
        mock_random.return_value = 0.199
        weights = [0.2, 0.8]
        expected = 'A'

        self.assertEqual(functions.choice(possibilities, weights), expected)

        # Correctly returns 'a' when probability of sampling 'a' is 0.2 and random value is greater than 0.2
        mock_cdf.return_value = [0.2, 1.0]
        mock_random.return_value = 0.201
        weights = [0.2, 0.8]
        expected = 'a'

        self.assertEqual(functions.choice(possibilities, weights), expected)

if __name__ == '__main__':
    unittest.main()
