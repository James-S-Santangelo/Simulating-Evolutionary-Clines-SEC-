import math
import random
import sys

from simulations import population
from simulations import functions


class Cell(object):
    """Controls cells in the landscape matrix"""

    max_mig_rate = sys.argv[2]
    """float: Maximum migration rate between any two populations"""

    max_K = 1000
    min_K = sys.argv[1]
    """int: Maximum and minimum carying capacity of cells across the matrix"""

    num_rows = 1
    num_cols = 15
    """int: Number of rows and number of columns in landscape matrix"""

    bot_prop = 0
    """float: Bottleneck proportion

    Proportion of alleles sampled when new population is created
    """

    max_create_prob = 1.0
    """float: Maximum probability of creating a new population"""

    max_s = 0
    """float: Maximum strength of selection acting against the cyanogenic phenotype"""

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
    def initialize_matrix(cls):
        """Creates the landscape matrix and intantiates each cell with a Cell object

        Args:
            None

        Returns:
            Matrix (:obj:'list' of :obj:'int'): 2D array storing instance of Cell at every position
        """
        Matrix = [[0] * Cell.num_cols for _ in range(Cell.num_rows)]
        for i in range(Cell.num_rows):
            for j in range(Cell.num_cols):
                Matrix[i][j] = cls(i, j, pop=False)
        return Matrix

    def weight_create_prob(self, max_pop_size):
        """Weights the probability of creating a population

        Based on maximum desired creation probability and a populations
        current size, determine the weighted probability that a given
        population will create a new one

        Args:
            max_pop_size (int): Maximum population size (i.e. carrying capacity)

        Returns:
            weighted_create_prob (float): weighted probability of creating a population

        Raises:
            ValueError: If max_create_prop is negative or greater than 1
            """
        if 0 <= Cell.max_create_prob <= 1:
            max_create_prob = float(Cell.max_create_prob)
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
        will_create = functions.choice(possibilities, weights)
        if will_create == '1':
            return True
        else:
            return False

    def empty_neighbors(self, Matrix):
        """Checks all adjacent cells in landscape matrix for populations

        Check whether adjacent cells in the landscape matrix contains populations,
        ignoring boundaries.

        Args:
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
                         if (-1 < self.i <= Cell.num_rows - 1
                             and -1 < self.j <= Cell.num_cols - 1
                             and (self.i != x2
                                  or self.j != y2)
                             and (0 <= x2 <= Cell.num_rows - 1)
                             and (0 <= y2 <= Cell.num_cols - 1))]

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

    def create_population(self, Matrix, max_pop_size):
        """Creates a population in randomly chosen, empty, adjacent cell

        Creates a new population in a randomly selected adjacent cell that
        contains no population. New populations are formed by sampling alleles
        from the focal population creating the new one. The number of alleles
        sampled depends on the bottleneck proportion.

        Args:
            Matrix (:obj:'list' of :obj:'int'): 2D array storing instance of Cell at every position.
            max_pop_size (int): Maximum population size (i.e. carrying capacity).

            Returns:
                None: If empty_neighbors returns empty list (i.e. all adjacent cells contain populations)
                or if will_create returns False (i.e. new population not being created).

            Raises:
                Exception: if bottleneck proportion is 0.
        """
        weighted_create_prob = self.weight_create_prob(max_pop_size)
        will_create = self.will_create(weighted_create_prob)
        Cell.bot_prop = float(Cell.bot_prop)

        if will_create:
            empty_neighbors = self.empty_neighbors(Matrix)

            if empty_neighbors:
                if Cell.bot_prop == 0:
                    raise Exception('You cannot create new empty populations!')
                else:
                    index = random.randint(0, len(empty_neighbors) - 1)
                    i, j = empty_neighbors[index][0], empty_neighbors[index][1]
                    Matrix[i][j].pop = True
                    new_size = int(math.ceil(Cell.bot_prop * self.population.size))
                    new_locus_A = self.population.sample_population(self.population.locus_A, new_size)
                    new_locus_B = self.population.sample_population(self.population.locus_B, new_size)
                    Matrix[i][j].population = population.Population(new_size, new_locus_A, new_locus_B)
            else:
                return None
        else:
            return None

    def real_migration_rate(self, source_y, source_x):
        """Calculates migration rate between populations

        Calculates the realized migration rate between two populations, which
        declines linearly with increasing distance between them.

        Args:
            source_y (int): Row index of source population (i.e. from which alleles are arriving)
            source_x (int): Column index of source population (i.e. from which alleles are arriving)

        Returns:
            mig_prop (float): Proportion of alleles migrating between focal and source population.
        """
        Cell.max_mig_rate = float(Cell.max_mig_rate)
        dist = (((source_y - self.i) ** 2) + ((source_x - self.j) ** 2)) ** (0.5)
        max_i = Cell.num_rows - 1
        max_j = Cell.num_cols - 1
        start_y = 0
        start_x = 0
        max_dist = (((max_j - start_x) ** 2) + ((max_i - start_y) ** 2)) ** (0.5)
        # max_dist = max_dist
        slope_mig = (0 - Cell.max_mig_rate) / (max_dist - 0)
        mig_prop = slope_mig * dist + Cell.max_mig_rate
        return round(mig_prop, 4)

    def real_K(self):
        """Calculates carrying capacity of cell in the matrix

        Calculates the relized carrying capacity of cell in the matrix, which
        declines linearly with increasing distance from the first
        initialized landscape cell.

        Args:
            None

        returns:
            K (float): Realized carrying capacity of landscape matrix cell.
        """
        Cell.min_K = int(Cell.min_K)
        Cell.max_K = int(Cell.max_K)
        start_y = 0
        start_x = 0
        dist = (((start_y - self.i) ** 2) + ((start_x - self.j) ** 2)) ** (0.5)
        max_i = Cell.num_rows - 1
        max_j = Cell.num_cols - 1
        max_dist = (((max_j - start_x) ** 2) + ((max_i - start_y) ** 2)) ** (0.5)

        # Slope and equation if gradient is from Rural to Urban
        slope_K = (Cell.min_K - Cell.max_K) / (max_dist - 0)
        K = Cell.max_K + slope_K * dist

        # Slope if gradient is from Urban to Rural
        # slope_K = (Cell.max_K - Cell.min_K) / (max_dist - 0)
        # K = Cell.min_K + slope_K * dist


        return int(math.ceil(K))

    def real_s(self):
        """Calculates selection coefficient of cells in the matrix

        Calculates the relized selection coefficient of cells in the matrix. Selection favours HCN+ genotypes in the rural-most population and declines linearly until the middle of the matrix. Selection then switches and increases linearly, increasingly favouring HCN- towards the urban-most population.

        Args:
            None

        returns:
            K (float): Realized selection coefficient of landscape matrix cell.
        """
        start_y = 0
        start_x = 0
        max_i = Cell.num_rows - 1
        max_j = Cell.num_cols - 1

        max_dist = (((max_j - start_x) ** 2) + ((max_i - start_y) ** 2)) ** (0.5)
        mid_dist = round(max_dist / 2)

        dist = (((start_y - self.i) ** 2) + ((start_x - self.j) ** 2)) ** (0.5)

        if dist == mid_dist:
            s = 0
            return round(s, 0)
        elif dist < mid_dist:
            dist_from_mid = mid_dist - dist
            max_s = float(Cell.max_s)
            slope_s = (0 - (-max_s)) / (mid_dist - 0)
            s = slope_s * dist - max_s
            return round(s, 5)
        elif dist > mid_dist:
            dist_from_mid = dist - mid_dist
            max_s = float(Cell.max_s)
            slope_s = (max_s - 0) / (max_dist - mid_dist)
            s = slope_s * dist_from_mid
            return round(s, 5)

    def source_population_info(self, pop_list, Matrix):
        """Collects migration, allele frequency, and population size info from all source populations

        Args:
            pop_list (:obj:'list' of :obj:'tuple'): Tuples with row and columns indices of all landscape
            cells with populations.
            Matrix (:obj:'list' of :obj:'int'): 2D array storing instance of Cell at every position

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

                mig_rate = self.real_migration_rate(source_y, source_x)
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

    def freq_after_selection(self):
        """Calculates frequency of alleles following selection

        Uses two-locus selection model to calculate the frequency of both dominant alleles the following generation. Also calculates linkage disequilibrium.

        Args:
            None

        Returns:
            pA1 (float): Frequency of 'A' allele the next generation
            pB1 (float): Frequency of 'B' allele the next generation
            D (float): Linkage disequilibrium between 'A' and 'B' loci.
        """

        pA = self.population.allele_freq(self.population.locus_A)
        pB = self.population.allele_freq(self.population.locus_B)

        qA = 1 - pA
        qB = 1 - pB

        s = self.real_s()
        rec = 0.5
        # h = 0

        Gametes = {'AB': pA * pB,
                   'Ab': pA * qB,
                   'aB': qA * pB,
                   'ab': qA * qB}

        D = (Gametes['AB'] * Gametes['ab']) - (Gametes['Ab'] * Gametes['aB'])

        Genotypes = {'AB_AB': [Gametes['AB'] ** 2, 1 - s],

                     'AB_Ab': [Gametes['AB'] * Gametes['Ab'], 1 - s],
                     'Ab_AB': [Gametes['AB'] * Gametes['Ab'], 1 - s],

                     'AB_aB': [Gametes['AB'] * Gametes['aB'], 1 - s],
                     'aB_AB': [Gametes['AB'] * Gametes['aB'], 1 - s],

                     'AB_ab': [Gametes['AB'] * Gametes['ab'], 1 - s],
                     'ab_AB': [Gametes['AB'] * Gametes['ab'], 1 - s],

                     'Ab_Ab': [Gametes['Ab'] ** 2, 1],

                     'Ab_aB': [Gametes['Ab'] * Gametes['aB'], 1 - s],
                     'aB_Ab': [Gametes['Ab'] * Gametes['aB'], 1 - s],

                     'Ab_ab': [Gametes['Ab'] * Gametes['ab'], 1],
                     'ab_Ab': [Gametes['Ab'] * Gametes['ab'], 1],

                     'aB_aB': [Gametes['aB'] ** 2, 1],

                     'aB_ab': [Gametes['aB'] * Gametes['ab'], 1],
                     'ab_aB': [Gametes['aB'] * Gametes['ab'], 1],

                     'ab_ab': [Gametes['ab'] ** 2, 1]}

        marg_fit = {'{0}'.format(key1): [] for key1 in Gametes.keys()}
        for i in Gametes.keys():
            for j in Genotypes.keys():
                if i == j.split("_")[0]:
                    marg_fit[i].append(Genotypes[j][1] * Gametes[j.split("_")[1]])
        for i, v in marg_fit.items():
            marg_fit[i] = sum(v)

        mean_fit = []
        for i1, v1 in Gametes.items():
            for i2, v2 in Gametes.items():
                con = '{0}_{1}'.format(i1, i2)
                mean_fit.append(v1 * v2 * Genotypes[con][1])
        mean_fit = sum(mean_fit)

        AB1 = (Gametes['AB'] * (marg_fit['AB'] / mean_fit)) - (rec * D * (Genotypes['AB_ab'][1] / mean_fit))
        Ab1 = (Gametes['Ab'] * (marg_fit['Ab'] / mean_fit)) + (rec * D * (Genotypes['AB_ab'][1] / mean_fit))
        aB1 = (Gametes['aB'] * (marg_fit['aB'] / mean_fit)) + (rec * D * (Genotypes['AB_ab'][1] / mean_fit))
        ab1 = (Gametes['ab'] * (marg_fit['ab'] / mean_fit)) - (rec * D * (Genotypes['AB_ab'][1] / mean_fit))

        pA1 = AB1 + Ab1
        pB1 = AB1 + aB1

        return round(pA1, 6), round(pB1, 6), round(D, 6)

    def freq_after_migration(self, migration_weighted, allele_weighted, allele_freq):
        """Calculates the frequency of an allele following migration

        Args:
            migration_weighted (float): Weighted migration rate across all source populations
            allele_weighted (float): Weighted frequency of allele across all source populations
            locus (:obj:'list' of :obj:'str'): List containing alleles at either 'A' or 'B' locus

        Returns:
            p_1 (float): Frequency of allele following migration
        """
        p_0 = allele_freq
        p_1 = ((1 - migration_weighted) * p_0) + (migration_weighted * allele_weighted)

        if p_1 < 0:
            p_1 = 0
        elif p_1 > 1.0:
            p_1 = 1.0

        return round(p_1, 4)

    def alleles_next_gen(self, pop_list, Matrix):
        """Calculates frequency of allele in next generation following migration and selection

        Args:
            pop_list (:obj:'list' of :obj:'tuple'): Tuples with row and columns indices of all landscape
            cells with populations.
            Matrix (:obj:'list' of :obj:'int'): 2D array storing instance of Cell at every position

        Returns:
            pA1 and pB1 (float): Frequency of A and B allele following migration and selection
        """
        weighting_factor_list = self.source_population_info(pop_list, Matrix)[3]
        migration_weighting_list = self.source_population_info(pop_list, Matrix)[0]
        allele_weighted_A_list = self.source_population_info(pop_list, Matrix)[1]
        allele_weighted_B_list = self.source_population_info(pop_list, Matrix)[2]

        migration_weighted = self.weighting(migration_weighting_list, weighting_factor_list)
        allele_weighted_A = self.weighting(allele_weighted_A_list, weighting_factor_list)
        allele_weighted_B = self.weighting(allele_weighted_B_list, weighting_factor_list)

        pA_after_sel = self.freq_after_selection()[0]
        pB_after_sel = self.freq_after_selection()[1]

        pA1 = self.freq_after_migration(migration_weighted, allele_weighted_A, pA_after_sel)
        pB1 = self.freq_after_migration(migration_weighted, allele_weighted_B, pB_after_sel)

        return pA1, pB1
