import math
import random

import Population
import Functions


class Cell(object):
    def __init__(self, i, j, pop):
        self.i = i  # Rows
        self.j = j  # Columns

        self.pop = False
        # self.population = Population.Population(0, [], [])

    @classmethod
    def initialize_matrix(cls, num_cols, num_rows):
        Matrix = [[0] * num_cols for _ in range(num_rows)]
        for i in range(num_rows):
            for j in range(num_cols):
                Matrix[i][j] = cls(i, j, pop=False)
        return Matrix

    def weight_create_prob(self, max_pop_size, max_create_prob):
        if 0 <= max_create_prob <= 1:
            max_create_prob = float(max_create_prob)
            slope = (max_create_prob - 0) / (max_pop_size - 0)
            size = self.population.size
            weighted_create_prob = size * slope
            return round(float(weighted_create_prob), 4)
        else:
            raise ValueError('Creation probability out of range!')

    @staticmethod
    def will_create(weighted_create_prob, possibilities):
        weighted_create_prob = float(weighted_create_prob)
        weights = [weighted_create_prob, (1 - weighted_create_prob)]
        will_create = Functions.choice(possibilities, weights)
        if will_create == '1':
            return True
        else:
            return False

    def empty_neighbors(self, num_rows, num_cols, Matrix):
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
        weighted_create_prob = self.weight_create_prob(max_pop_size, max_create_prob)
        will_create = self.will_create(weighted_create_prob, possibilities='10')

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
            value_weighted = sum(values_list[i] * factor_list[i] / sum(factor_list) for i in range(len(values_list)))
            return round(value_weighted, 4)

    def freq_after_migration(self, migration_weighted, allele_weighted, locus):
        p_0 = self.population.allele_freq(locus)
        p_1 = ((1 - migration_weighted) * p_0) + (migration_weighted * allele_weighted)

        if p_1 < 0:
            p_1 = 0
        elif p_1 > 1.0:
            p_1 = 1.0

        return round(p_1, 4)
