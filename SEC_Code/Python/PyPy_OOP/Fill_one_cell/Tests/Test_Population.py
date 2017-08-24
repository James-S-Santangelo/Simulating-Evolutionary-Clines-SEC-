import unittest
from parameterized import parameterized
from mock import patch
import math

import Population


class TestPopulation(unittest.TestCase):

    def setUp(self):
        N = 10
        pA = 0.5
        pB = 0.5
        qA = 1 - pA
        qB = 1 - pB
        locus_A = (['A'] * int(N * pA)) + (['a'] * int(round(N * qA)))
        locus_B = (['B'] * int(N * pB)) + (['b'] * int(round(N * qB)))
        self.pop = Population.Population(N, locus_A, locus_B)

    def tearDown(self):
        del self.pop

    def test_init(self):
        self.assertEqual(self.pop.size, 10)
        self.assertListEqual(self.pop.locus_A, ['A', 'A', 'A', 'A', 'A', 'a', 'a', 'a', 'a', 'a'])
        self.assertListEqual(self.pop.locus_B, ['B', 'B', 'B', 'B', 'B', 'b', 'b', 'b', 'b', 'b'])

    @parameterized.expand([("100", 15), ("10", 10), ("5", 7)])
    def test_pop_growth(self, K, expected):
        self.assertEqual(self.pop.pop_growth(self.pop.size, K, math.log(1.5)), expected)

    def test_allele_freq(self):
        self.assertEqual(self.pop.allele_freq(self.pop.locus_A), 0.5)

        locus = (['B'] * 0) + (['b'] * 10)
        freq = self.pop.allele_freq(locus)
        self.assertEqual(freq, 0)
        assert isinstance(freq, float)

        locus = []
        self.assertRaises(ZeroDivisionError, self.pop.allele_freq, locus)

    def test_phenotype(self):
        self.assertEqual(self.pop.phenotype(self.pop.locus_A, self.pop.locus_B), 0.5625)

        self.pop.locus_B = (['B'] * 0) + (['b'] * 10)
        self.assertEqual(self.pop.phenotype(self.pop.locus_A, self.pop.locus_B), 0)

        self.pop.locus_B = (['B'] * 10) + (['b'] * 0)
        self.assertEqual(self.pop.phenotype(self.pop.locus_A, self.pop.locus_B), 0.75)

        self.pop.locus_B = []
        self.assertRaises(ZeroDivisionError, self.pop.phenotype, self.pop.locus_A, self.pop.locus_B)

    def test_sample_alleles(self):
        self.assertListEqual(self.pop.sample_alleles(0, 'Aa'), (['A'] * 0) + (['a'] * 10))
        self.assertListEqual(self.pop.sample_alleles(1, 'Aa'), (['A'] * 10) + (['a'] * 0))

    @patch('random.choice')
    def test_sample_population(self, mock_rand_choice):

        expected = ['A', 'A', 'A', 'a', 'a', 'a']
        mock_rand_choice.side_effect = ['A', 'A', 'A', 'a', 'a', 'a']

        N = 6
        actual = self.pop.sample_population(self.pop.locus_A, N)
        self.assertEqual(actual, expected)
        mock_rand_choice.assert_called_with(self.pop.locus_A)

    @patch('Functions.choice')
    def test_sample_alleles_mocked(self, mock_choice):

        expected = ['A', 'A', 'A', 'A', 'A', 'a', 'a', 'a', 'a', 'a']
        mock_choice.side_effect = ['A', 'A', 'A', 'A', 'A', 'a', 'a', 'a', 'a', 'a']
        actual = self.pop.sample_alleles(0.5, 'Aa')

        self.assertEqual(actual, expected)
        mock_choice.assert_called_with('Aa', [0.5, 0.5])


if __name__ == '__main__':
    unittest.main()
