#!/usr/local/bin/python2.7

from Functions_PyPy import *
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
	Distance_Dic = Distance_Mig(x_mat, y_mat, max_mig_rate)

	K_dict = varyK(Distance_Dic, max_K, min_K, x_mat, y_mat)

	print os.getpid(), (max_mig_rate, bot)

	# Open csv to which data will be written
	os.chdir(export_path)
	datestring = datetime.strftime(datetime.now(), '%Y%m%d')
	with open(datestring + "_SEC_Drift.Migration_(pA%.2f)(pB%.2f).csv" % (pA, pB), "wb") as f:
		writer = csv.writer(f)
		writer.writerow(["Sim","x","y","Population","Pop_size","Generation","pA","pB","Acyan", "Mig_rate", "K", "r", "max_p_create", "bot", "Mat.full"])

		# Loop for 'sims' iterations
		print_progress(0, sims, prefix='', suffix='', decimals=1, bar_length=100)
		for s in range(sims):

			sim = {}

			# Initialize pop_counter, 'pops' and 'alleles' dictionaries and
			# Matrix, with population 1 initializes in Matrix
			pop_counter = [1]
			pops = OrderedDict({'1':[]}) # Re-initialize dictionary to store populations
			alleles = OrderedDict({'1':{'A':[],'B':[],'S':[N]}}) # Re-initialize dictionary to store allele lists
			Matrix = create_matrix(x_mat, y_mat)
			Matrix[0][0] = 1

			# Initialize allele lists from which population 1 will sample its
			# alleles in generation 1.
			locus_A = (['A'] * int(N * pA) ) + (['a'] * int(round(N * qA)))  # Re-initialize initial allele lists.
			locus_B = (['B'] * int(N * pB) ) + (['b'] * int(round(N * qB)))

			# Run 'cline' function.
			cline(locus_A,locus_B, steps, N, max_p_create, pops, alleles, bot, Matrix, max_K, K_dict, r, max_mig_rate, pop_counter, Distance_Dic) # Run cline function

			# Append results from 'cline' function 'sim' dictionary as new entry.
			# Keys correspond to interation.
			sim[s] = pops # Append results to 'sim' dictionary

			write_to_csv(writer, sim)

			# Ensure sim is deleted to avoid overconsumption of RAM
			del sim

			print_progress(s + 1, sims, prefix='', suffix='', decimals=1, bar_length=100)

if __name__ == '__main__':
	simulate()