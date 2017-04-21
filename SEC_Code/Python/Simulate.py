#!/usr/bin/env python

from Functions import *
from Parameters import *

def simulate():
	'''Generate 'sims' simulations of 'cline' function.

	Runs 'cline' function 'sims' times, re-initializing the simulator
	with each iteration. The results of each iteration are appended to
	the 'sim' dictionary, which is then exported as a csv to desired
	directory when all iterations have completed.

	Parameters:
	All parameters imported from Parameters.py

	Returns:
	csv with results to directory of choice.
	'''
	# Initilize 'sim' dictionary to store results and 'Distance_Dic' dictionary
	# to store distances and migration rates for simulations.
	sim = {}
	Distance_Dic = Distance_Mig(x_mat, y_mat, max_mig_rate)

	# Loop for 'sims' iterations
	for s in range(sims):

		# Initialize pop_counter, 'pops' and 'alleles' dictionaries and
		# Matrix, with population 1 initializes in Matrix
		pop_counter = [1]
		pops = OrderedDict({'1':[]}) # Re-initialize dictionary to store populations
		alleles = OrderedDict({'1':{'A':[],'B':[],'S':[N]}}) # Re-initialize dictionary to store allele lists
		Matrix = np.zeros((x_mat, y_mat), dtype = 'int')
		Matrix[0, 0] = 1

		# Initialize allele lists from which population 1 will sample its
		# alleles in generation 1.
		locus_A = (['A'] * int(N * pA) ) + (['a'] * int(round(N * qA)))  # Re-initialize initial allele lists.
		locus_B = (['B'] * int(N * pB) ) + (['b'] * int(round(N * qB)))

		# Run 'cline' function.
		cline(locus_A,locus_B, steps, N, max_p_create, pops, alleles, bot, Matrix, K, r, max_mig_rate, pop_counter, Distance_Dic) # Run cline function

		# Append results from 'cline' function 'sim' dictionary as new entry.
		# Keys correspond to interation.
		sim[s] = pops # Append results to 'sim' dictionary

	# Changed directory and write results to csv.
	os.chdir(sys.argv[3])
	write_to_csv(max_mig_rate, bot, sim)

if __name__ == '__main__':
	simulate() # Call function if script is run as main.



