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
        max_create_prob = float(max_create_prob)
        slope = (max_create_prob - 0) / (max_pop_size - 0)
        size = self.population.size
        weighted_create_prob = size * slope
        return float(weighted_create_prob)

    @staticmethod
    def real_create_prob(weighted_create_prob, possibilities, num_of_items):
        weights = [weighted_create_prob, (1 - weighted_create_prob)]
        real_create_prob = [Functions.choice(possibilities, weights) for i in range(num_of_items)]
        if real_create_prob == '1':
            return True
        else:
            return False

    def empty_neighbors(self, rows, cols, Matrix):
        neighbors = [(x2, y2) for x2 in range(self.i - 1, self.i + 2)
                     for y2 in range(self.j - 1, self.j + 2)
                     if (-1 < self.i <= rows - 1
                         and -1 < self.j <= cols - 1
                         and (self.i != x2
                              or self.j != y2)
                         and (0 <= x2 <= rows - 1)
                         and (0 <= y2 <= cols - 1))]

        empty_neighbors = []
        for item in neighbors:
            i, j = item[0], item[1]
            if not Matrix[i][j].pop:
                empty_neighbors.append(item)

        return empty_neighbors

    def create_population(self, rows, cols, Matrix, max_p_create, bot, K, r):
        weighted_create_prob = self.weight_create_prob(max_pop_size=K, max_create_prob=max_p_create)
        real_create_prob = self.real_create_prob(weighted_create_prob, possibilities='10', number_of_items=1)

        if real_create_prob:
            empty_neighbors = self.empty_neighbors(self, rows, cols, Matrix)

            if empty_neighbors:
                Nsam = random.randint(0, len(empty_neighbors) - 1)
                i, j = empty_neighbors[Nsam][0], empty_neighbors[Nsam][1]
                Matrix[i][j].pop = True
                new_size = int(math.ceil(bot * self.population.size))
                new_locus_A = self.population.sample_population(self.population.locus_A, new_size)
                new_locus_B = self.population.sample_population(self.population.locus_B, new_size)
                Matrix[i][j].population = Population.Population(new_size, new_locus_A, new_locus_B)
            else:
                pass
        else:
            pass

    def real_migration_rate(self, source_y, source_x, rows, cols, max_mig_rate):
        dist = (((source_y - self.i) ** 2) + ((source_x - self.j) ** 2)) ** (0.5)
        max_i = rows - 1
        max_j = cols - 1
        start_y = 0
        start_x = 0
        max_dist = (((max_j - start_x) ** 2) + ((max_i - start_y) ** 2)) ** (0.5)
        slope_mig = (max_mig_rate - 0) / (max_dist - 0)
        mig_prop = max_mig_rate - slope_mig * dist
        return round(mig_prop, 4)

    def real_K(self, rows, cols, min_K, max_K):
        dist = (((0 - self.i) ** 2) + ((0 - self.j) ** 2)) ** (0.5)
        max_i = rows - 1
        max_j = cols - 1
        start_y = 0
        start_x = 0
        max_dist = (((max_j - start_x) ** 2) + ((max_i - start_y) ** 2)) ** (0.5)
        slope_K = (min_K - max_K) / (max_dist - 0)
        K = max_K + slope_K * dist
        return int(math.ceil(K))

    def source_population_info(self, pop_list, Matrix, rows, cols, max_mig_rate):

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

                mig_rate = self.real_migration_rate(source_y, source_x, rows, cols, max_mig_rate)
                allele_freq_A = source_pop.population.allele_freq(source_pop.population.locus_A)
                allele_freq_B = source_pop.population.allele_freq(source_pop.population.locus_B)
                population_size = source_pop.population.size

                migration_weighted.append(mig_rate)
                allele_weighted_A.append(allele_freq_A)
                allele_weighted_B.append(allele_freq_B)
                pop_size.append(population_size)

        return migration_weighted, allele_weighted_A, allele_weighted_B, pop_size

    @staticmethod
    def weighting(values, factor):
            weighted_value = sum(values[i] * factor[i] / sum(factor) for i in range(len(values)))
            return weighted_value

    def freq_after_migration(self, weighted_migration, weighted_allele, locus):
        p_0 = self.population.allele_freq(locus)
        p_1 = ((1 - weighted_migration) * p_0) + (weighted_migration * weighted_allele)

        if p_1 < 0:
            p_1 = 0
        elif p_1 > 1.0:
            p_1 = 1.0

        return p_1
