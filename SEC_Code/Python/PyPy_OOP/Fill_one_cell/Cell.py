import math
import random

import Population
import Functions


class Cell(object):
    """Controls cells in the landscape matrix"""

    def __init__(self, i, j, pop):
        """Class constructor function

        Args:
            i (int): Row number of cell in landscape matrix
            j (int): Column nujmber of cell in landscape matrix
            Pop (bool): True if cell contains population. False if no population
        """
        self.i = i  # Rows
        self.j = j  # Columns
        self.pop = False

    @classmethod
    def initialize_matrix(cls, num_cols, num_rows):
        """Creates the landscape matrix and intantiates each cell with a Cell object

        Args:
            num_cols (int): Number of columns in landscape
            num_rows (int): Number of rows in landscape

        Returns:
            Matrix (:obj:'list' of :obj:'int'): 2D array storing instance of Cell at every position
        """
        Matrix = [[0] * num_cols for _ in range(num_rows)]
        for i in range(num_rows):
            for j in range(num_cols):
                Matrix[i][j] = cls(i, j, pop=False)
        return Matrix

    def weight_create_prob(self, max_pop_size, max_create_prob):
        """Weights the probability of creating a population

        Based on maximum desired creation probability and a populations
        current size, determine the weighted probability that a given
        population will create a new one

        Args:
            max_pop_size (int): Maximum population size (i.e. carrying capacity)
            max_create_prob (float): Maximum probability of creating a new population

        Returns:
            weighted_create_prob (float): weighted probability of creating a population

        Raises:
            ValueError: If max_create_prop is negative or greater than 1
            """
        if 0 <= max_create_prob <= 1:
            max_create_prob = float(max_create_prob)
            slope = (max_create_prob - 0) / (max_pop_size - 0)
            size = self.population.size
            weighted_create_prob = size * slope
            return round(float(weighted_create_prob), 4)
        else:
            raise ValueError('Creation probability out of range!')

    @staticmethod
    def will_create(weighted_create_prob):
        """Determine if current population will create a new one

        From the weighted creation probability, determine if the
        current population successfully creates a new one.

        Args:
            weighted_create_prob (float): Weighted probability of creating a new population

        Returns:
            bool: True if new population will be created. False if not.
        """
        possibilities = '10'
        weighted_create_prob = float(weighted_create_prob)
        weights = [weighted_create_prob, (1 - weighted_create_prob)]
        will_create = Functions.choice(possibilities, weights)
        if will_create == '1':
            return True
        else:
            return False

    def empty_neighbors(self, num_rows, num_cols, Matrix):
        """Checks all adjacent cells in landscape matrix for populations

        Check whether adjacent cells in the landscape matrix contains populations,
        ignoring boundaries.

        Args:
            num_rows (int): Number of rows in landscape matrix
            num_cols (int): Number of columns in landscape matrix
            Matrix (:obj:'list' of :obj:'int'): 2D array storing instance of Cell at every position.

        Returns:
            empty_neighbors (:obj:'list' of :obj:'tuple'): Tuples are positions (i, j) of adjacent
            cells that do not contain populations

        Raises:
            Exception: If focal cell contain no population
        """
        if self.pop:
            neighbors = [(x2, y2) for x2 in range(self.i - 1, self.i + 2)
                         for y2 in range(self.j - 1, self.j + 2)
                         if (-1 < self.i <= num_rows - 1
                             and -1 < self.j <= num_cols - 1
                             and (self.i != x2
                                  or self.j != y2)
                             and (0 <= x2 <= num_rows - 1)
                             and (0 <= y2 <= num_cols - 1))]

            empty_neighbors = []
            for item in neighbors:
                i, j = item[0], item[1]
                if not Matrix[i][j].pop:
                    empty_neighbors.append(item)
                else:
                    continue

            return empty_neighbors

        else:
            raise Exception('Cell without population cannot create new populations (obviously')

    def create_population(self, num_rows, num_cols, Matrix, max_create_prob, bot_prop, max_pop_size):
        """Creates a population in randomly chosen, empty, adjacent cell

        Creates a new population in a randomly selected adjacent cell that
        contains no population. New populations are formed by sampling alleles
        from the focal population creating the new one. The number of alleles
        sampled depends on the bottleneck proportion.

        Args:
            num_rows (int): Number of rows in the landscape matrix
            num_cols (int): Number of columns in the landscape matrix
            Matrix (:obj:'list' of :obj:'int'): 2D array storing instance of Cell at every position.
            max_create_prob (float): Maximum probability of creating a new population.
            bot_prop (float): Bottleneck proportion. Proportion of alleles sampled upon creation of
            new population.
            max_pop_size (int): Maximum population size (i.e. carrying capacity).

            Returns:
                None: If empty_neighbors returns empty list (i.e. all adjacent cells contain populations)
                or if will_create returns False (i.e. new population not being created).

            Raises:
                Exception: if bottleneck proportion is 0.
        """
        weighted_create_prob = self.weight_create_prob(max_pop_size, max_create_prob)
        will_create = self.will_create(weighted_create_prob)

        if will_create:
            empty_neighbors = self.empty_neighbors(num_rows, num_cols, Matrix)

            if empty_neighbors:
                if bot_prop == 0:
                    raise Exception('You cannot create new empty populations!')
                else:
                    index = random.randint(0, len(empty_neighbors) - 1)
                    i, j = empty_neighbors[index][0], empty_neighbors[index][1]
                    Matrix[i][j].pop = True
                    new_size = int(math.ceil(bot_prop * self.population.size))
                    new_locus_A = self.population.sample_population(self.population.locus_A, new_size)
                    new_locus_B = self.population.sample_population(self.population.locus_B, new_size)
                    Matrix[i][j].population = Population.Population(new_size, new_locus_A, new_locus_B)
            else:
                return None
        else:
            return None

    def real_migration_rate(self, source_y, source_x, num_rows, num_cols, max_mig_rate):
        """Calculates migration rate between populations

        Calculates the realized migration rate between two populations, which
        declines linearly with increasing distance between them.

        Args:
            source_y (int): Row index of source population (i.e. from which alleles are arriving)
            source_x (int): Column index of source population (i.e. from which alleles are arriving)
            num_rows (int): Number of rows in landscape matrix.
            num_cols (int): Number of columns in landscape matrix
            max_mig_rate (float): Maximum possible migration between two populations

        Returns:
            mig_prop (float): Proportion of alleles migrating between focal and source population.
        """
        dist = (((source_y - self.i) ** 2) + ((source_x - self.j) ** 2)) ** (0.5)
        max_i = num_rows - 1
        max_j = num_cols - 1
        start_y = 0
        start_x = 0
        max_dist = (((max_j - start_x) ** 2) + ((max_i - start_y) ** 2)) ** (0.5)
        # max_dist = max_dist
        slope_mig = (0 - max_mig_rate) / (max_dist - 0)
        mig_prop = slope_mig * dist + max_mig_rate
        return round(mig_prop, 4)

    def real_K(self, num_rows, num_cols, min_K, max_K):
        """Calculates carrying capacity of cell in the matrix

        Calculates the relized carrying capacity of cell in the matrix, which
        declines linearly with increasing distance from the first
        initialized landscape cell.

        Args:
            num_rows (int): Number of rows in the landscape matrix
            num_cols (int): Number of columns in the landscape matrix
            min_K (int): Minimum carrying capacity at cell furthest away from start
            max_K (int): Maximum possible carrying capacity (i.e. at start).

            returns:
                K (float): Realized carrying capacity of landscape matrix cell.
        """
        start_y = 0
        start_x = 0
        dist = (((start_y - self.i) ** 2) + ((start_x - self.j) ** 2)) ** (0.5)
        max_i = num_rows - 1
        max_j = num_cols - 1
        max_dist = (((max_j - start_x) ** 2) + ((max_i - start_y) ** 2)) ** (0.5)
        slope_K = (min_K - max_K) / (max_dist - 0)
        K = max_K + slope_K * dist
        return int(math.ceil(K))

    def source_population_info(self, pop_list, Matrix, num_rows, num_cols, max_mig_rate):
        """Collects migration, allele frequency, and population size info from all source populations

        Args:
            pop_list (:obj:'list' of :obj:'tuple'): Tuples with row and columns indices of all landscape
            cells with populations.
            Matrix (:obj:'list' of :obj:'int'): 2D array storing instance of Cell at every position
            num_rows (int): Number of rows in landscape matrix.
            num_cols (int): Number of columns in landscape matrix
            max_mig_rate (float): Maximum possible migration between two populations

        Returns:
            migration_weighted (:obj:`list` of :obj:`float`): Realized migration rates from all source
            populations
            allele_weighted_A (:obj:`list` of :obj:`float`): Frequency of 'A' allele in all source populations
            allele_weighted_B (:obj:`list` of :obj:`float`): Frequency of 'B' allele in all source populations
            pop_size (:obj:`list` of :obj:`int`): Population sizes of all source populations
        """
        migration_weighted = []
        allele_weighted_A = []
        allele_weighted_B = []
        pop_size = []

        for pop in pop_list:
            source_y, source_x = pop[0], pop[1]

            if source_y == self.i and source_x == self.j:
                pass
            else:

                source_pop = Matrix[source_y][source_x]

                mig_rate = self.real_migration_rate(source_y, source_x, num_rows, num_cols, max_mig_rate)
                allele_freq_A = source_pop.population.allele_freq(source_pop.population.locus_A)
                allele_freq_B = source_pop.population.allele_freq(source_pop.population.locus_B)
                population_size = source_pop.population.size

                migration_weighted.append(mig_rate)
                allele_weighted_A.append(allele_freq_A)
                allele_weighted_B.append(allele_freq_B)
                pop_size.append(population_size)

        return migration_weighted, allele_weighted_A, allele_weighted_B, pop_size

    @staticmethod
    def weighting(values_list, factor_list):
        """Weights a list of focal values based on list of weighting factors

        Args:
            values_list (:obj:'list' of :obj:'float'): List of elements to be weighted
            factor_list (:obj:'list' of :obj:'int'): List of elements to be used as weighting factors

        Returns:
            value_weighted (float): Weighted value
        """
        value_weighted = sum(values_list[i] * factor_list[i] / sum(factor_list) for i in range(len(values_list)))
        return round(value_weighted, 4)

    def freq_after_migration(self, migration_weighted, allele_weighted, locus):
        """Calculates the frequency of an allele following migration

        Args:
            migration_weighted (float): Weighted migration rate across all source populations
            allele_weighted (float): Weighted frequency of allele across all source populations
            locus (:obj:'list' of :obj:'str'): List containing alleles at either 'A' or 'B' locus

        Returns:
            p_1 (float): Frequency of allele following migration
        """
        p_0 = self.population.allele_freq(locus)
        p_1 = ((1 - migration_weighted) * p_0) + (migration_weighted * allele_weighted)

        if p_1 < 0:
            p_1 = 0
        elif p_1 > 1.0:
            p_1 = 1.0

        return round(p_1, 4)
