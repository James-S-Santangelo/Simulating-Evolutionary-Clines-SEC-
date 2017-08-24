import math
import random

import Functions


class Population(object):

    def __init__(self, size, locus_A, locus_B):
        self.size = size
        self.locus_A = locus_A
        self.locus_B = locus_B

    def pop_growth(self, size, K, r):
        K = float(K)
        new_size = self.size * K / (self.size + (K - self.size) * math.exp(-r))
        if new_size < K or new_size == K:
            return int(math.ceil(new_size))
        elif new_size > K:
            return int(math.floor(new_size))

    def phenotype(self, locus_A, locus_B):
        qA = 1 - self.allele_freq(self.locus_A)
        qB = 1 - self.allele_freq(self.locus_B)
        mut = qA ** 2 + qB ** 2 - (qA ** 2 * qB ** 2)
        WT = 1 - mut
        return WT

    @staticmethod
    def allele_freq(locus):
        sum = 0
        for i in locus:
            if i == 'A' or i == 'B':
                sum += 1
        p = sum / float(len(locus))
        return round(p, 4)

    @staticmethod
    def sample_population(locus, N):
        new_locus = [random.choice(locus) for _ in range(N)]
        return new_locus

    def sample_alleles(self, allele_freq, possibilities):
        num_of_items = self.size
        weights = [allele_freq, (1 - allele_freq)]
        new_locus = [Functions.choice(possibilities, weights) for i in range(num_of_items)]
        return new_locus
