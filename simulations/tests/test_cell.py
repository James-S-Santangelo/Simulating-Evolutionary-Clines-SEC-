import unittest
import math
from mock import patch, Mock

from simulations.cell import Cell
from simulations.population import Population


class TestCell(unittest.TestCase):
    """Tests Cell.py"""

    def setUp(self):
        """setUp to be run before each test function

        Creates a landscape with 2 rows and 3 columns. Initializes
        mocked populations of varying sizes and allele frequencies into
        3 of the 6 landscape matrix cells, according to the following pattern

        [[1, 2, 0],
         [0, 0, 3]]

        where cells that are not 0 contain populations. Numbers correspond to
        the mocked populations created in function.
        """
        Cell.num_rows = 2
        Cell.num_cols = 3

        self.Matrix = Cell.initialize_matrix()

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

    def tearDown(self):
        """tearDown function that deletes landscape matrix after each test"""
        del self.Matrix

    def test_initialize_matrix(self):
        """Tests initialize_matrix method"""
        Cell.num_rows = 2
        Cell.num_cols = 5
        Matrix = Cell.initialize_matrix()

        # Checks that initialize_matrix successfully instantiates Cell objects in cells.
        for i in range(Cell.num_rows):
            for j in range(Cell.num_cols):
                assert isinstance(Matrix[i][j], Cell)

    def test_init(self):
        """Tests contructor method"""

        # Tests that row and column indices are successfully assigned.
        for i in range(Cell.num_rows):
            for j in range(Cell.num_cols):
                self.assertEqual(self.Matrix[i][j].i, i)
                self.assertEqual(self.Matrix[i][j].j, j)

        # Tests that cells initialized with population have pop set to True
        self.assertIs(self.Matrix[0][0].pop, True)
        self.assertIs(self.Matrix[0][1].pop, True)
        self.assertIs(self.Matrix[1][2].pop, True)

        # Tests that cells NOT initialized with population have pop set to False
        self.assertIs(self.Matrix[0][2].pop, False)
        self.assertIs(self.Matrix[1][0].pop, False)
        self.assertIs(self.Matrix[1][1].pop, False)

        # Tests that one of the cells with a population is correctly assigned an instance of Population class
        assert isinstance(self.Matrix[0][0].population, Population)
        # Tests that instance of Population class has correctly assigned size attribute.
        self.assertEqual(self.Matrix[0][0].population.size, 10)

    def test_weight_create_prob(self):
        """Tests weight_create_prob method"""

        # Zero returned when max_prob_create is 0
        Cell.max_create_prob = 0
        self.assertEqual(self.Matrix[0][0].weight_create_prob(1000), 0)
        # Correct computation of creation probability
        Cell.max_create_prob = 0.5
        self.assertEqual(self.Matrix[0][0].weight_create_prob(100), 0.05)
        # Computes current probability as max_prob when current size is equal to max_size
        Cell.max_create_prob = 0.5
        self.assertEqual(self.Matrix[0][0].weight_create_prob(10), 0.5)
        # Correct computation when max_size is 1
        Cell.max_create_prob = 1
        self.assertEqual(self.Matrix[0][0].weight_create_prob(1000), 0.01)

        # Raises ValueError if creation_prob is greater than 1
        Cell.max_create_prob = 1.5
        self.assertRaises(ValueError, self.Matrix[0][0].weight_create_prob, 1000)
        # Raises ValueError if creation_prob is negative
        Cell.max_create_prob = -0.5
        self.assertRaises(ValueError, self.Matrix[0][0].weight_create_prob, 1000)
        # Raises attribute error if function called on cell with no population
        Cell.max_create_prob = 1
        self.assertRaises(AttributeError, self.Matrix[0][2].weight_create_prob, 1000)

    @patch('simulations.functions.choice')
    def test_will_create(self, mock_choice):
        """Tests will_create method

        args:
            mock_choice (:obj:'mocked function'): Mocked version of choice functon from Functions.py
        """

        # True when weighted_create_prob is 1
        mock_choice.return_value = '1'
        self.assertIs(self.Matrix[0][1].will_create(1), True)

        # False when weighted_create_prob is 0
        mock_choice.return_value = '0'
        self.assertIs(self.Matrix[0][1].will_create(0), False)

        # True when weighted_create_prob is 0.5 and choice returns 1
        mock_choice.return_value = '1'
        self.assertIs(self.Matrix[0][1].will_create(0.5), True)

    def test_empty_neighbors(self):
        """Tests empty_neighbors method"""

        # Correctly returns empty adjacent cells. Tested on all mocked populations
        expected_list = [(1, 0), (1, 1)]
        self.assertListEqual(self.Matrix[0][0].empty_neighbors(self.Matrix), expected_list)
        expected_list = [(0, 2), (1, 0), (1, 1)]
        self.assertListEqual(self.Matrix[0][1].empty_neighbors(self.Matrix), expected_list)
        expected_list = [(0, 2), (1, 1)]
        self.assertListEqual(self.Matrix[1][2].empty_neighbors(self.Matrix), expected_list)

        # Fails since no population in cell calling neighbors method
        self.assertRaises(Exception, self.Matrix[0][2].empty_neighbors, self.Matrix)

        # Return empty list if all adjacent cells are occupied
        mock_pop_4 = Mock(spec=Population, size=10, locus_A=(['A'] * 5) + (['a'] * 5), locus_B=(['B'] * 5) + (['b'] * 5))
        mock_pop_5 = Mock(spec=Population, size=10, locus_A=(['A'] * 5) + (['a'] * 5), locus_B=(['B'] * 5) + (['b'] * 5))

        self.Matrix[1][0] = mock_pop_4
        self.Matrix[1][1] = mock_pop_5

        self.Matrix[1][0].pop = True
        self.Matrix[1][1].pop = True

        expected_list = []
        self.assertListEqual(self.Matrix[0][0].empty_neighbors(self.Matrix), expected_list)

    @patch('simulations.population.Population')
    @patch('random.randint')
    @patch('simulations.cell.Cell.empty_neighbors')
    @patch('simulations.cell.Cell.will_create')
    @patch('simulations.cell.Cell.weight_create_prob')
    def test_create_population(self, mock_weighted, mock_will_create, mock_neighbors, mock_random, mock_new_pop):
        """Tests create_population method

        Args:
            mock_weighted (:obj:'mocked method'): Mocked weight_create_prob method
            mock_will_create (:obj:'mocked method'): Mocked will_create method
            mock_neighbors (:obj:'mocked method'): Mocked empty_neighbors method
            mock_random (:obj:'mocked method'): Mocked randint method from random module
            mock_new_pop (:obj:'mocked class'): Mocked instance of Population class
        """

        max_pop_size = 1000

        # Testing that None is returned if 'mock_will_create' is False
        mock_weighted.return_value = 0
        mock_will_create.return_value = False
        Cell.max_create_prob = 0  # Not relevant since mocked
        Cell.max_pop_size = 100  # Not relevant since mocked
        Cell.bot_prop = 1.0

        self.assertIs(self.Matrix[0][0].create_population(self.Matrix, max_pop_size), None)

        # Testing that None is returned if 'mock_will_create' is True but all adjacent cells contain populations (i.e. 'empty_neighbors' returns empty list)
        mock_weighted.return_value = 1.0
        mock_will_create.return_value = True
        mock_neighbors.return_value = []
        Cell.max_create_prob = 1.0  # Not relevant since mocked
        Cell.max_pop_size = 100  # Not relevant since mocked
        Cell.bot_prop = 1.0

        self.assertIs(self.Matrix[0][0].create_population(self.Matrix, max_pop_size), None)

        # Testing that Exception is raised if the bottleneck proportion is 0
        mock_weighted.return_value = 1.0
        mock_will_create.return_value = True
        mock_neighbors.return_value = [(0, 2), (1, 1)]
        mock_random.return_value = 0
        Cell.max_create_prob = 1.0  # Not relevant since mocked
        Cell.max_pop_size = 100  # Not relevant since mocked
        Cell.bot_prop = 0

        self.assertRaises(Exception, self.Matrix[0][0].create_population, self.Matrix, max_pop_size)

        # Testing that a new population with correct size and loci is created in only one adjacent cell. Also tests conversion of string to float.
        mock_weighted.return_value = 1.0
        mock_will_create.return_value = True
        mock_neighbors.return_value = [(0, 2), (1, 1)]
        mock_random.return_value = 0
        Cell.max_create_prob = 1.0  # Not relevant since mocked
        Cell.max_pop_size = 100  # Not relevant since mocked
        Cell.bot_prop = 0.2

        self.Matrix[1][2].population.sample_population.side_effect = [(['A'] * 3) + (['a'] * 4), (['B'] * 4) + (['b'] * 3)]

        new_size = int(math.ceil(Cell.bot_prop * self.Matrix[1][2].population.size))
        new_locus_A = (['A'] * 3) + (['a'] * 4)
        new_locus_B = (['B'] * 4) + (['b'] * 3)

        mock_new_pop.return_value = Population(new_size, new_locus_A, new_locus_B)

        Cell.bot_prop = "0.2"
        self.Matrix[1][2].create_population(self.Matrix, max_pop_size)

        self.assertIs(self.Matrix[0][2].pop, True)
        self.assertIs(self.Matrix[1][1].pop, False)
        self.assertIsInstance(self.Matrix[0][2].population, Population)
        mock_new_pop.assert_called_once_with(new_size, new_locus_A, new_locus_B)
        self.assertEqual(self.Matrix[0][2].population.size, 7)
        self.assertListEqual(self.Matrix[0][2].population.locus_A, (['A'] * 3) + (['a'] * 4))
        self.assertListEqual(self.Matrix[0][2].population.locus_B, (['B'] * 4) + (['b'] * 3))
        assert isinstance(Cell.bot_prop, float)

    def test_real_migration_rate(self):
        """Tests real_migration_rate method"""

        # Correctly returns migration rate from adjacent cell
        Cell.max_mig_rate = 0.05
        source_y = 0
        source_x = 0
        self.assertEqual(self.Matrix[0][1].real_migration_rate(source_y, source_x), 0.0276)

        # Correactly returns 0 when maximum migration rate is 0
        Cell.max_mig_rate = 0
        self.assertEqual(self.Matrix[0][0].real_migration_rate(source_y, source_x), 0)

        # Correctly returns migration rate from non-adjacent cell
        Cell.max_mig_rate = 0.05
        source_y = 0
        source_x = 2
        self.assertEqual(self.Matrix[0][0].real_migration_rate(source_y, source_x), 0.0053)

        # Correctly returns migration rate of 0 when cell is at maximum dtstance
        Cell.max_mig_rate = 0.05
        source_y = 1
        source_x = 2
        self.assertEqual(self.Matrix[0][0].real_migration_rate(source_y, source_x), 0)

        # Correctly converts values supplied as strings to float
        Cell.max_mig_rate = "0.05"
        source_y = 1
        source_x = 2
        self.assertEqual(self.Matrix[0][0].real_migration_rate(source_y, source_x), 0)
        assert isinstance(Cell.max_mig_rate, float)

    def test_real_K(self):
        """Tests real_K function"""

        # Correctly return carrying capacity for all cells with populations defined in setUp
        Cell.min_K = 10
        Cell.max_K = 1000
        self.assertEqual(self.Matrix[0][0].real_K(), 1000)
        self.assertEqual(self.Matrix[1][2].real_K(), 10)
        self.assertEqual(self.Matrix[0][1].real_K(), 558)

        # Correctly returns carrying capacity when K is not set to vary
        Cell.min_K = 1000
        Cell.max_K = 1000
        self.assertEqual(self.Matrix[0][0].real_K(), 1000)
        self.assertEqual(self.Matrix[1][2].real_K(), 1000)
        self.assertEqual(self.Matrix[0][1].real_K(), 1000)


    def test_real_s(self):
        """Tests real_s function"""

        # Correctly return selection coefficient for all cells with populations defined in setUp
        Cell.max_s = 0.001
        self.assertEqual(self.Matrix[0][0].real_s(), 0)
        self.assertEqual(self.Matrix[1][2].real_s(), 0.001)
        self.assertEqual(self.Matrix[0][1].real_s(), 0.00045)


    @patch('simulations.cell.Cell.real_migration_rate')
    def test_source_population_info(self, mock_real):
        """Tests source_population_info method

        Args:
            mock_real (:onj:'Mocked method'): Mocked real_migration_rate method
        """

        pop_list = [(i, j) for j in range(Cell.num_cols) for i in range(Cell.num_rows) if self.Matrix[i][j].pop]

        Cell.max_mig_rate = 0.05

        mock_real.side_effect = [0.0276, 0.0]
        self.Matrix[0][1].population.allele_freq.side_effect = [0.7059, 0.1176]
        self.Matrix[1][2].population.allele_freq.side_effect = [0, 0.6875]

        # Correctly returns all 4 lists with weighted values
        self.assertTupleEqual(self.Matrix[0][0].source_population_info(pop_list, self.Matrix), ([0.0276, 0.0], [0.7059, 0], [0.1176, 0.6875], [17, 32]))

    def test_weighting(self):
        """Tests weighting method"""

        # Correctly returns original value when values and factors don't vary
        values_list = [30, 30]
        factors_list = [10, 10]
        self.assertEqual(self.Matrix[0][0].weighting(values_list, factors_list), 30)

        # Correctly return midpoint of values when one is 2x the other and factors doesn't vary
        values_list = [30, 60]
        factors_list = [10, 10]
        self.assertEqual(self.Matrix[0][0].weighting(values_list, factors_list), 45)

        # Correctly returns weighted values when initial values and factors more complex
        values_list = [0.1176, 0.6875]
        factors_list = [17, 32]
        self.assertEqual(self.Matrix[0][0].weighting(values_list, factors_list), 0.4898)

        # Correctly returns weighted values when one value is 0
        values_list = [0.0276, 0.0]
        factors_list = [17, 32]
        self.assertEqual(self.Matrix[0][0].weighting(values_list, factors_list), 0.0096)

    @patch('simulations.cell.Cell.real_s')
    def test_freq_after_selection(self, mock_s):
        """Tests freq_after_selection method"""

        # Correctly return allale freuencies following weak selection
        self.Matrix[0][1].population.allele_freq.side_effect = [0.7059, 0.1176]
        Cell.max_s = 0.001

        mock_s.return_value = 0.00045

        self.assertTupleEqual(self.Matrix[0][1].freq_after_selection(), (0.705894, 0.117562, 0))

        # Correctly returns same allele frequencies when no selection
        self.Matrix[0][1].population.allele_freq.side_effect = [0.7059, 0.1176]
        Cell.max_s = 0.001

        mock_s.return_value = 0.0

        self.assertTupleEqual(self.Matrix[0][1].freq_after_selection(), (0.7059, 0.1176, 0))

        # Correclty return allele frequencies when selection is strong
        self.Matrix[0][1].population.allele_freq.side_effect = [0.7059, 0.1176]
        Cell.max_s = 0.01

        mock_s.return_value = 0.01

        self.assertTupleEqual(self.Matrix[0][1].freq_after_selection(), (0.705765, 0.116762, 0))

    def test_freq_after_migration(self):
        """Tests freq_after_migration method"""

        # Correctly returns new allele frequency of mock_pop_2 based on migration from other populations
        self.Matrix[0][1].population.allele_freq.return_value = 0.1176
        migration_weighted = 0.0206
        allele_weighted = 0.6429
        locus = self.Matrix[0][1].population.locus_B

        self.assertEqual(self.Matrix[0][1].freq_after_migration(migration_weighted, allele_weighted, locus), 0.1284)

        # Correctly returns new allele frequency of mock_pop_1 as current frequency if there is no migration.
        self.Matrix[0][0].population.allele_freq.return_value = 0.5
        migration_weighted = 0
        allele_weighted = 0.4898
        locus = self.Matrix[0][0].population.locus_A

        self.assertEqual(self.Matrix[0][0].freq_after_migration(migration_weighted, allele_weighted, locus), 0.5)


if __name__ == '__main__':
    unittest.main()


