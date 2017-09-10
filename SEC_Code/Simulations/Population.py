import math
import random

import Functions


class Population(object):
    """Controls populations assigned to cells in the landscape matrix"""

    def __init__(self, size, locus_A, locus_B):
        """Class constructor function

        Args:
            size (int): Size of the population. Equal to length of locus_A and locus_B
            locus_A (:obj:'list' of :obj:'str'): List storing dominant ('A') and recessive ('a') alleles at locus_A
            locus_B (:obj:'list' of :obj:'str'): List storing dominant ('B') and recessive ('b') alleles at locus_B
        """
        self.size = size
        self.locus_A = locus_A
        self.locus_B = locus_B

    def pop_growth(self, K, r):
        """Calculates new population size based on logistic growth model

        Args:
            K (int): Maximum size of population (i.e. carrying capacity)
            r (float): Instantaneous rate of population increase (i.e. growth rate)

        Returns:
            new_size (int): Rounded up new_size is less than or equal to K, rounded down otherwise.
        """
        K = float(K)
        new_size = self.size * K / (self.size + (K - self.size) * math.exp(-r))
        if new_size <= K:
            return int(math.ceil(new_size))
        elif new_size > K:
            return int(math.floor(new_size))

    def phenotype(self, locus_A, locus_B):
        """Calculates the frequency of the cyanogenic phenotype

        Args:
            locus_A (:obj:'list' of :obj:'str'): List storing dominant ('A') and recessive ('a') alleles at locus_A
            locus_B (:obj:'list' of :obj:'str'): List storing dominant ('B') and recessive ('b') alleles at locus_B

        Returns:
            WT (float): Frequency of wild-type (cyanogenic) phenotype.
        """
        qA = 1 - self.allele_freq(self.locus_A)
        qB = 1 - self.allele_freq(self.locus_B)
        mut = qA ** 2 + qB ** 2 - (qA ** 2 * qB ** 2)
        WT = 1 - mut
        return WT

    @staticmethod
    def allele_freq(locus):
        """Calculates frequency of allele

        Args:
            locus (:obj:'list' of :obj:'str'): List storing dominant and recessive alleles at locus. Can be either locus_A or locus_B

        Returns:
            p (float): Frequency of allele
        """
        sum = 0
        for i in locus:
            if i == 'A' or i == 'B':
                sum += 1
        p = sum / float(len(locus))
        return round(p, 4)

    @staticmethod
    def sample_population(locus, N):
        """Randomly samples alleles from a locus

        Args:
            locus (:obj:'list' of :obj:'str'): List storing dominant and recessive alleles at locus. Can be either locus_A or locus_B
            N (int): Number of alleles to sample

        Returns:
            new_locus (:obj:'list' of :obj:'str'): New list of length N containing dominant and recessive alleles sampled from locus.
        """
        new_locus = [random.choice(locus) for _ in range(N)]
        return new_locus

    def sample_alleles(self, allele_freq, possibilities):
        """Samples alleles from specified probability distribution

        Probability distribution comes from the predicted frequency of alleles
        in the next generation following migration and selection.

        Args:
            allele_freq (float): Frequency of dominant allele. Either 'A' or 'B'.
            possibilities (str): Possible alleles to sample. Either 'Aa' or 'Bb', where 'Aa' represents choices 'A' and 'a'.

        Returns:
            new_locus (:obj:'list' of :obj:'str'): New list of length N containing dominant and recessive alleles sampled from locus.
        """
        num_of_items = self.size
        weights = [allele_freq, (1 - allele_freq)]
        new_locus = [Functions.choice(possibilities, weights) for i in range(num_of_items)]
        return new_locus
