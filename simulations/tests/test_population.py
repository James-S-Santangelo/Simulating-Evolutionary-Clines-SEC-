import unittest
from parameterized import parameterized
from mock import patch
import math

from simulations import population


class TestPopulation(unittest.TestCase):
    """Tests Population.py"""

    def setUp(self):
        """Setup to be run before each test

        Sets up a single instance of the population class to be used for testing.
        """
        N = 10
        pA = 0.5
        pB = 0.5
        qA = 1 - pA
        qB = 1 - pB
        locus_A = (['A'] * int(N * pA)) + (['a'] * int(round(N * qA)))
        locus_B = (['B'] * int(N * pB)) + (['b'] * int(round(N * qB)))
        self.pop = population.Population(N, locus_A, locus_B)

    def tearDown(self):
        """tearDown function that deletes landscape matrix after each test"""
        del self.pop

    def test_init(self):
        """Tests contructor method"""

        # Tests that Population instance has size and locus attributes consistent with those assigned in setUp
        self.assertEqual(self.pop.size, 10)
        self.assertListEqual(self.pop.locus_A, ['A', 'A', 'A', 'A', 'A', 'a', 'a', 'a', 'a', 'a'])
        self.assertListEqual(self.pop.locus_B, ['B', 'B', 'B', 'B', 'B', 'b', 'b', 'b', 'b', 'b'])

    @parameterized.expand([("100", 15), ("10", 10), ("5", 7)])
    def test_pop_growth(self, K, expected):
        """Tests pop_growth method

        Args:
            K (int): Maximum population size (i.e. carrying capacity). Specified by string in parameterized.expand.
            expected (int): Expected population size returned from pop_growth method. Specified by int in parameterized.expand.
        """

        # Test that 3 inputs from parameterized.expand provide expected output.
        self.assertEqual(self.pop.pop_growth(K), expected)

    def test_allele_freq(self):
        """Tests allele_freq method"""

        # Correctly calculated frequency of 'A' allele from Population setup in setUp method
        self.assertEqual(self.pop.allele_freq(self.pop.locus_A), 0.5)

        # Correctly calculates frequency of 'B' allele as 0 and identifies return as float object
        locus = (['B'] * 0) + (['b'] * 10)
        freq = self.pop.allele_freq(locus)
        self.assertEqual(freq, 0)
        assert isinstance(freq, float)

        # Correctly reaises error when given empty list as input
        locus = []
        self.assertRaises(ZeroDivisionError, self.pop.allele_freq, locus)

    def test_phenotype(self):
        """Tests phenotype method"""

        # Correctly calculated phenotype frequency of Population setup in setUp method
        self.assertEqual(self.pop.phenotype(self.pop.locus_A, self.pop.locus_B), 0.5625)

        # Correctly calculates frequency of phenotype when dominant allele is 0
        self.pop.locus_B = (['B'] * 0) + (['b'] * 10)
        self.assertEqual(self.pop.phenotype(self.pop.locus_A, self.pop.locus_B), 0)

        # Correctly calculates frequency of phenotype when recessive allele is 0
        self.pop.locus_B = (['B'] * 10) + (['b'] * 0)
        self.assertEqual(self.pop.phenotype(self.pop.locus_A, self.pop.locus_B), 0.75)

        # Correctly raises error when allele list is empty
        self.pop.locus_B = []
        self.assertRaises(ZeroDivisionError, self.pop.phenotype, self.pop.locus_A, self.pop.locus_B)

    def test_sample_alleles(self):
        """Tests sample_alleles method"""

        # Correctly return all recessive or dominant alleles and given probabilities of 0 or 1, respectively.
        self.assertListEqual(self.pop.sample_alleles(0, 'Aa'), (['A'] * 0) + (['a'] * 10))
        self.assertListEqual(self.pop.sample_alleles(1, 'Aa'), (['A'] * 10) + (['a'] * 0))

    @patch('random.choice')
    def test_sample_population(self, mock_rand_choice):
        """Tests sample_population method

        Args:
            mock_rand_choice (:obj:'mocked function'): Mocked version of the choice function from the random module.
        """

        # Correctly returns list with half 'A' and half 'a' and asserts that mock function is successfully calleld with allele list.
        expected = ['A', 'A', 'A', 'a', 'a', 'a']
        mock_rand_choice.side_effect = ['A', 'A', 'A', 'a', 'a', 'a']

        N = 6
        actual = self.pop.sample_population(self.pop.locus_A, N)
        self.assertEqual(actual, expected)
        mock_rand_choice.assert_called_with(self.pop.locus_A)

    @patch('simulations.functions.choice')
    def test_sample_alleles_mocked(self, mock_choice):
        """Tests sample_allele with mocked function

        Args:
            mock_choice (:obj:'mocked method'): Mocked version of the choice function from the Functions module.
        """

        # Correctly returns list with half 'A' and half 'a' and asserts that mock function is successfully calleld with allele list.
        expected = ['A', 'A', 'A', 'A', 'A', 'a', 'a', 'a', 'a', 'a']
        mock_choice.side_effect = ['A', 'A', 'A', 'A', 'A', 'a', 'a', 'a', 'a', 'a']
        actual = self.pop.sample_alleles(0.5, 'Aa')

        self.assertEqual(actual, expected)
        mock_choice.assert_called_with('Aa', [0.5, 0.5])


if __name__ == '__main__':
    unittest.main()
