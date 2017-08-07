import math
import random

from Functions import cdf, choice

class Population( object ):

	def __init__( self, size, locus_A, locus_B ):
		self.size = size
		self.locus_A = locus_A
		self.locus_B = locus_B

	def pop_growth( self, size, K, r ):
		K = float( K )
		self.size = self.size * K / ( self.size + ( K - self.size ) * math.exp( -r ) )
		self.size = int( math.ceil( self.size ) )
		return self.size

	def phenotype( self, locus_A, locus_B ):
		qA = 1 - self.allele_freq( self.locus_A )
		qB = 1 - self.allele_freq( self.locus_B )
		mut= qA ** 2 + qB ** 2 - ( qA ** 2 * qB ** 2 )
		WT = 1 - mut
		return WT

	@staticmethod
	def allele_freq( locus ):
		sum = 0
		for i in locus:
			if i == 'A' or i == 'B':
				sum += 1
		p = sum / float(len(locus))
		return p

	@staticmethod
	def sample_population( locus, N ):
		new_locus = [ random.choice( locus ) for _ in range( N )]
		return new_locus

	def sample_alleles(self, allele, r, possibilities):
		number_of_items_to_pick = self.size # Number to sample. Corresponds to next generation's size.
		weights = [allele, (1 - allele)] # Sampling probabilities. Returned by 'alleles_next_gen'
		return [choice(possibilities, weights) for i in range(number_of_items_to_pick)]