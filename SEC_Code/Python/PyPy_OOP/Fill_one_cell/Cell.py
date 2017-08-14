import math
import random

import Population
import Functions


def create_matrix(cols, rows):
    Matrix = [[0] * cols for _ in range(rows)]
    for i in range(rows):
        for j in range(cols):
            Matrix[i][j] = Cell(i, j, pop=False)
    return Matrix


class Cell(object):
    def __init__(self, i, j, pop):
        self.i = i
        self.j = j

        self.pop = False
        self.population = Population.Population(0, [], [])

    def prob_create(self, K, max_p_create):
        max_p_create = float(max_p_create)
        m = (max_p_create - 0) / (K - 0)
        Size = self.population.size
        p_create = Size * m
        return float(p_create)

    def create_population(self, rows, cols, Matrix, max_p_create, bot, K, r):
        population = '10'
        number_of_items_to_pick = 1
        p_create = self.prob_create(K, max_p_create)
        weights = [p_create, (1 - p_create)]
        create = [Functions.choice(population, weights) for i in range(number_of_items_to_pick)]
        if create[0] == '1':
            neighbors = [(x2, y2) for x2 in range(self.i - 1, self.i + 2)
                         for y2 in range(self.j - 1, self.j + 2)
                         if (-1 < self.i <= rows - 1
                             and -1 < self.j <= cols - 1
                             and (self.i != x2
                                  or self.j != y2)
                             and (0 <= x2 <= rows - 1)
                             and (0 <= y2 <= cols - 1))]

            Nlist_red = []
            for item in neighbors:
                i, j = item[0], item[1]
                if not Matrix[i][j].pop:
                    Nlist_red.append(item)

            if not Nlist_red:
                pass
            else:
                Nsam = random.randint(0, len(Nlist_red) - 1)
                i, j = Nlist_red[Nsam][0], Nlist_red[Nsam][1]
                Matrix[i][j].pop = True
                new_size = int(math.ceil(bot * self.population.size))
                new_locus_A = self.population.sample_population(self.population.locus_A, new_size)
                new_locus_B = self.population.sample_population(self.population.locus_B, new_size)
                Matrix[i][j].population = Population.Population(new_size, new_locus_A, new_locus_B)
        else:
            pass

    def Distance_calc(self, y, x, rows, cols, max_mig_rate, min_K, max_K):

        dist = (((y - self.i) ** 2) + ((x - self.j) ** 2)) ** (0.5)

        max_i = rows - 1
        max_j = cols - 1

        max_dis = (((max_j - x) ** 2) + ((max_i - y) ** 2)) ** (0.5)

        m_mig = (max_mig_rate - 0) / (max_dis - 0)

        mig_prop = max_mig_rate - m_mig * dist

        m_k = (min_K - max_K) / (max_dis - 0)

        K = max_K + m_k * dist

        return round(mig_prop, 4), int(math.ceil(K))

    def alleles_next_gen(self, rows, cols, pop_list, max_mig_rate, min_K, max_K):

        migration_weighted = []
        allele_weighted_A = []
        allele_weighted_B = []
        pop_size = []

        for pop in pop_list:
            y, x = pop[0], pop[1]
            if y == self.i and x == self.j:
                # print y, self.i, x, self.j, "SAME"
                pass
            else:
                # print y, self.i, x, self.j, "DIFFERENT"
                # print (self.i, self.j), (y, x), self.Distance_calc(y, x, rows, cols)[0]
                migration_weighted.append(self.Distance_calc(y, x, rows, cols, max_mig_rate, min_K, max_K)[0])
                allele_weighted_A.append(self.population.allele_freq(self.population.locus_A))
                allele_weighted_B.append(self.population.allele_freq(self.population.locus_B))
                pop_size.append(self.population.size)
        migration_weighted = sum(migration_weighted[g] * pop_size[g] / sum(pop_size) for g in range(len(migration_weighted)))
        allele_weighted_A = sum(allele_weighted_A[g] * pop_size[g] / sum(pop_size) for g in range(len(allele_weighted_A)))
        allele_weighted_B = sum(allele_weighted_B[g] * pop_size[g] / sum(pop_size) for g in range(len(allele_weighted_B)))

        pA1 = ((1 - migration_weighted) * self.population.allele_freq(self.population.locus_A)) + (migration_weighted * allele_weighted_A)
        pB1 = ((1 - migration_weighted) * self.population.allele_freq(self.population.locus_B)) + (migration_weighted * allele_weighted_B)

        # Ensure that negative probabilities do not occur.
        if pA1 < 0:
            pA1 = 0
        elif pA1 > 1.0:
            pA1 = 1.0
        if pB1 < 0:
            pB1 = 0
        elif pB1 > 1.0:
            pB1 = 1.0

        return pA1, pB1
