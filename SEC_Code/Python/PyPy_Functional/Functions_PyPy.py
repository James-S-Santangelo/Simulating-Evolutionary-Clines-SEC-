#!/usr/local/bin/python2.7

# Modules used throughout script
import random
from collections import OrderedDict
import csv
import time
from datetime import datetime
import os
import itertools
import math
import bisect
import sys as sys
from Parameters import x_mat, y_mat, max_mig_rate


def print_progress(iteration, total, prefix='', suffix='', decimals=1, bar_length=100):
	"""
	Call in a loop to create terminal progress bar
	@params:
		iteration   - Required  : current iteration (Int)
		total       - Required  : total iterations (Int)
		prefix      - Optional  : prefix string (Str)
		suffix      - Optional  : suffix string (Str)
		decimals    - Optional  : positive number of decimals in percent complete (Int)
		bar_length  - Optional  : character length of bar (Int)
	"""
	str_format = "{0:." + str(decimals) + "f}"
	percents = str_format.format(100 * (iteration / float(total)))
	filled_length = int(round(bar_length * iteration / float(total)))
	bar = '#' * filled_length + '-' * (bar_length - filled_length)

	sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),

	if iteration == total:
		sys.stdout.write('\n')
	sys.stdout.flush()

def sample_population_A(locus_A, N):
	'''Samples N alleles from locus_A allele list.

	Parameters:
	1. locus_A: List containing alleles 'a' and 'A'
	2. N: Number of alleles to sample from locus_A (i.e. population size)

	Returns:
	List containing N sampled alleles for locus A.
	'''
	new_locus_A = [random.choice(locus_A) for _ in range(N)]
	return new_locus_A

def sample_population_B(locus_B, N):
	'''Samples N alleles from locus_B allele list.

	Parameters:
	1. locus_B: List containing alleles 'b' and 'B'
	2. N: Number of alleles to sample from locus_B (i.e. population size)

	Returns:
	List containing N sampled alleles for locus B.
	'''
	new_locus_B = [random.choice(locus_B) for _ in range(N)]
	return new_locus_B

def create_matrix(x_mat, y_mat):
	Matrix = [[0] * y_mat for i in range(x_mat)]
	return Matrix

def migration_rate(Distance_Dic, max_mig_rate):
	'''Calculates migration rate between populations.

	Migration rate is calculated based on paiwise distances between
	populations. Migration rate decline linearly with increasing distance
	such that when distance = max(distance_in_matrix), migration rate = 0.
	The migration rate is appended as the second element in the 'Distance_Dic'
	dictionary returned by the 'Distance_Mig' function.

	Parameters:
	1. Distance_Dic: Dictionary retuned by the 'Distance_Mig' function where
	keys correspond to pairwise population location is the matrix and values
	are lists where the first element is the distance between populations and
	the second is the migration rate appended from this function.
	2. max_mig_rate: Desired migration rate when distance = 0

	Returns:
	Migration rate appended to 'Distance_Dic' dictionary.
	'''
	max_dis = max(Distance_Dic.values()) # Max distance between populations
	m = (max_mig_rate - 0)/(max_dis[0] - 0) # Slope. Assumes close to no migration at max distance. Relized migration at max distance may be slightly greater than 0 due to rounding.
	for Dkey, Dvalue in Distance_Dic.items():
		Mig_prop = max_mig_rate - m*Dvalue[0]
		Distance_Dic[Dkey].append(round(Mig_prop, 4)) # Append migration rate to 'Distance_Dic' dictionary

# Function to calculate distances between all pairwise combinations of cells in matrix. Also calls
# the migration rate between populations and appends this to the distance dictionary.
def Distance_Mig(x_mat, y_mat, max_mig_rate):
	'''Creates dictionary with pairwise distance and migratin rates between populations.

	Dictionary keys correspond to pairwise population location in the matrix and
	values are lists where the first element is the distance between populations and
	the second is the migration rate appended 'migration_rate' function.

	Parameters:
	1. x_mat: Number of columns in matrix
	2. y_mat: Number of rows in matrix
	max_mig_rate: Desired migration rate when distance = 0

	Returns:
	'Distance_Dic' dictionary.
	'''
	rows = [i for i in range(x_mat)]
	cols = [i for i in range(y_mat)]
	matelem = [(i,j) for i in rows for j in cols]
	dis = [[i,j] for i in matelem for j in matelem]
	for i in dis:
		if i[0][0] == i[1][0]: # Within rows
			dis1 = abs(i[1][1] - i[0][1])
			i.append(dis1)
		elif i[0][1] == i[1][1]: # Within columns
			dis2 = abs(i[1][0] - i[0][0])
			i.append(dis2)
		else: # Diagonals
			dis3 = (((i[1][1] - i[0][1])**2) + ((i[1][0] - i[0][0])**2))**(0.5)
			i.append(dis3)
	Distance_Dic = {'{0}.{1}'.format(key1, key2):[round(key3, 2)] for key1, key2, key3 in dis}
	migration_rate(Distance_Dic, max_mig_rate)
	return Distance_Dic

def cdf(weights):
	total = sum(weights)
	result = []
	cumsum = 0
	for w in weights:
		cumsum += w
		result.append(cumsum / total)
	return result

def choice(population, weights):
	assert len(population) == len(weights)
	cdf_vals = cdf(weights)
	x = random.random()
	idx = bisect.bisect(cdf_vals, x)
	return population[idx]

def varyK(Distance_Dic, max_K, min_K, x_mat, y_mat):

	max_dis = max(Distance_Dic.values())[0]
	m = (min_K - max_K) / (max_dis - 0)

	Matrix = create_matrix(x_mat, y_mat)
	Matrix_indices = [[(i, j) for j, v in enumerate(sublist)] for i, sublist in enumerate(Matrix)]
	Flat_indices = [item for sublist in Matrix_indices for item in sublist]
	start = (0, 0)

	groups = []
	for i in Flat_indices:
		con = '{0}.{1}'.format(start, i)

		if con in Distance_Dic.keys():
			groups.append([con, Distance_Dic[con][0]])

	groups = [[i[0].split(".")[1], i[1]] for i in groups]

	K_dict = {'{0}'.format(key1):[key2, max_K + m * key2] for key1, key2 in groups}

	return K_dict

def sample_alleles_A(pA1, Akey, Avalue, alleles, r, K_dict, Matrix):
	'''Sample alleles at locus A from infinite allele pool.

	Sample N alleles from probability distribution based on expected frequency
	of 'A' allele in next generation. N is the next generation population size based on
	logistic population growth. Expected frequency is a function of current allele
	frequency in population and mean frequency of migratory alleles, weighted by size
	of immigrant populations.

	Parameters:
	1. pA1: Probability of sampling allele 'A'. Returned by 'alleles_next_gen()'
	2. Akey: Index used to cycle through keys in alleles dictionary.
	Note keys correspond to populations. See 'cline' function.
	3. Avalue: Index used to cycle through values in alleles dictionary. See 'cline function.
	4. alleles: Dictionary used to stored lists of alleles population size for each population.
	Updated every generation.
	5. r: natural rate of increase
	6. K: Carrying capacity (i.e. maximum sustainable population size)

	Returns:
	list containing sampled 'A' alleles
	'''
	population = 'Aa' # Possible alleles to sample
	number_of_items_to_pick = pop_growth(r, Akey, Avalue, K_dict, Matrix) # Number to sample. Corresponds to next generation's size.
	weights = [pA1, (1 - pA1)] # Sampling probabilities. Returned by 'alleles_next_gen'
	return [choice(population, weights) for i in range(number_of_items_to_pick[0])] # Return list of newly sampled alleles. Becomes allele pool in the next generation

def sample_alleles_B(pB1, Akey, Avalue, alleles, r, K_dict, Matrix):
	'''Sample alleles at locus B from infinite allele pool.

	Sample N alleles from probability distribution based on expected frequency
	of 'B' allele in next generation. N is the next generation population size based on
	logistic population growth. Expected frequency is a function of current allele
	frequency in population and mean frequency of migratory alleles, weighted by size
	of immigrant populations.

	Parameters:
	1. pB1: Probability of sampling allele 'B'. Returned by 'alleles_next_gen()'
	2. Akey: Index used to cycle through keys in alleles dictionary.
	Note keys correspond to populations. See 'cline' function.
	3. Avalue: Index used to cycle through values in alleles dictionary. See 'cline function.
	4. alleles: Dictionary used to stored lists of alleles population size for each population.
	Updated every generation.
	5. r: natural rate of increase
	6. K: Carrying capacity (i.e. maximum sustainable population size)

	Returns:
	list containing sampled 'B' alleles
	'''
	population = 'Bb' # Possible alleles to sample
	number_of_items_to_pick = pop_growth(r, Akey, Avalue, K_dict, Matrix) # Number to sample. Corresponds to next generation's size.
	weights = [pB1, (1 - pB1)] # Sampling probabilities. Returned by 'alleles_next_gen'
	return [choice(population, weights) for i in range(number_of_items_to_pick[0])] # Return list of newly sampled alleles. Becomes allele pool in the next generation

def matrix_full(Matrix):
	'''Check if Matrix has been filled with populations

	Parameters:
	1. Matrix: Matrix: m x n dimensional matrix, initialized at the outset of each simulation.

	Returns:
	1. 0 if Matrix still has empty cells
	2. 1 if Matrix is full
	'''
	if any(0 in j for j in Matrix):
		return 0
	else:
		return 1

def alleles_next_gen(Akey, pop_list, alleles, Matrix, Distance_Dic):
	'''Determines expected frequency of 'A' and 'B' alleles in next generation.

	Computes the probability of sampling an 'A' or 'B' allele in the
	next generation. Based on current allele frequency and mean frequency of
	allele across all populations from which alleles can immigrate. Frequency
	of alleles from immigrant populations are weighted by their size so that larger
	populations contribute proportionally more alleles. See README.

	Parameters:
	1. Akey: Index used to cycle through keys in alleles dictionary.
	Note keys correspond to populations. See 'cline' function.
	2. pop_list: list containing all current populations in existence
	3. alleles: Dictionary used to stored lists of alleles population size for each population.
	Updated every generation.
	4. Matrix: m x n dimensional matrix, initialized at the outset of each simulation.
	Contains empty cells that may become filled with populations.
	5. Distance_Dic: Dictionary containing distance and migration rates between populations.

	Returns:
	1. pA1: Probability of sampling 'A' allele in next generation.
	2. pA2: Probability of sampling 'B' allele in next generation.
	'''
	to_pop = [(i, sublist.index(int(Akey))) for i, sublist in enumerate(Matrix) if int(Akey) in sublist][0] # Location of focal population in matrix.
	migration_weighted = [] # List holding migration rates
	allele_weighted_A = [] # List holding frequency of 'A' alleles
	allele_weighted_B = [] # List holding frequency of 'B' alleles
	Size = [] # List holding population sizes
	for i in pop_list:
		if Akey == i:
			pass
		else:
			from_pop = [(j, sublist.index(int(i))) for j, sublist in enumerate(Matrix) if int(i) in sublist][0] # Location of source population in matrix
			con = str(to_pop) + '.' + str(from_pop) # Create concatenated string from current focal population and source population
			migration_weighted.append(Distance_Dic[con][1]) # Append migration rate to list
			allele_weighted_A.append(allele_freq(alleles[i]['A'])) # Append alleles 'A' frequency to list
			allele_weighted_B.append(allele_freq(alleles[i]['B'])) # Append allele 'B' frequency to list
			Size.append(alleles[i]['S'][0]) # Append population size to list
	migration_weighted = sum(migration_weighted[g] * Size[g] / sum(Size) for g in range(len(migration_weighted))) # Weighted migration rate
	allele_weighted_A = sum(allele_weighted_A[g] * Size[g] / sum(Size) for g in range(len(allele_weighted_A))) # Weighted allele 'A'
	allele_weighted_B = sum(allele_weighted_B[g] * Size[g] / sum(Size) for g in range(len(allele_weighted_B))) # Weighted allele 'B'
	pA1 = ((1 - migration_weighted) * allele_freq(alleles[Akey]['A'])) + (migration_weighted * allele_weighted_A) # Probability of sampling 'A' in next generation
	pB1 = ((1 - migration_weighted) * allele_freq(alleles[Akey]['B'])) + (migration_weighted * allele_weighted_B) # Probability of sampling 'B' in next generation

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

def allele_freq(locus):
	'''Calculate frequency of dominant (i.e. 'A' or 'B') allele.

	Parameters:
	1. locus: List containing alleles

	Returns:
	p: allele frequency as float.
	'''
	sum = 0
	for i in locus:
		if i == 'A' or i == 'B':
			sum += 1
	p = sum / float(len(locus))
	return p

def pop_growth(r, Akey, Avalue, K_dict, Matrix):
	'''Calculates population size in next generation based on logistic model.

	Parameters:
	1. r: natural rate of increase
	2. Akey: Index used to cycle through keys in alleles dictionary.
	Note keys correspond to populations. See 'cline' function.
	3. Avalue: Index used to cycle through values in alleles dictionary. See 'cline function.
	4. K: Carrying capacity (i.e. maximum sustainable population size)

	Returns:
	Next generation's population size as interger (rounded up).
	'''
	pop = str([(i, sublist.index(int(Akey))) for i, sublist in enumerate(Matrix) if int(Akey) in sublist][0])

	size = Avalue['S'][0] # Retrieve size of population. 'Akey' allows indexing of alleles dictionary in cline function
	K = K_dict[pop][1]
	new_size = size * K/(size + (K - size) * math.exp(-r)) # Calculates the proportional reduction of population growth rate based on desired carrying capacity ('K'). At 'K', growth rate = 1 = no change
	return [int(math.ceil(new_size))]

def bottle(bot, Akey, Avalue):
	'''Calcualates number of alleles to sample based on bottleneck proportion

	Parameters:
	1. bot: Desired bottleneck proportion (i.e. size of bottleneck upon creation of new
	population where 1.0 = no bottleneck)
	2. key: Index used to cycle through keys in alleles dictionary.
	Note keys correspond to populations. See 'cline' function.
	3. Avalue: Index used to cycle through values in alleles dictionary. See 'cline function.

	Returns:
	Number of alleles to sample as integer (rounded up).
	'''
	return int(math.ceil(bot * Avalue['S'][0]))

def prob_create(max_K, max_p_create, Akey, Avalue):
	'''Calculates probability of creating new population based on current size.

	Based on desired maximum probability of population creation, determines
	the probability that a population will create a new one based on its current
	size. The probability declines linearly with decreasing size such that the
	probability of creation is greatest when a population is at carrying capacity (K).

	Parameters:
	1. K: Carrying capacity
	2. max_p_create: maximum probability of creating a new population
	3. Akey: Index used to cycle through keys in alleles dictionary.
	Note keys correspond to populations. See 'cline' function.
	4. Avalue: Index used to cycle through values in alleles dictionary. See 'cline function.

	Returns:
	probability of creating a new population as float.
	'''
	m = (max_p_create - 0)/(max_K - 0) # Slope. Assumes close to no migration at max distance. Relized migration at max distance may be slightly greater than 0 due to rounding.
	Size = Avalue['S'][0]
	p_create = Size * m
	return float(p_create)

# Create new population as empty list and add to 'pops' dictionary. Also create four lists of alleles
# sampled from pool of alleles from population that generated the new one. Alleles added to 'alleles'
# dictionary. Also adds population to matrix. This function first evaluates whether a population will
# be created then randomly selects a vacant neighboring cell where population will go. If no cells are
# vacant, the function passes.
def create_population(max_p_create, max_K, Akey, Avalue, pops, alleles, bot, Matrix, pop_counter):
	'''Creates new populations and initializes them in existing dictionaries

	Create new population as empty list and add to 'pops' dictionary. Also create
	two lists of alleles sampled from pool of alleles from population that generated the
	new one. Allele lists added to 'alleles' dictionary. Also adds population to matrix.
	This function first evaluates whether a population will be created then randomly
	selects a vacant neighboring cell where population will go. If no cells are
	vacant, the function passes.

	Parameters:
	1. max_p_create: maximum probability of creating a new population
	2. Akey: Index used to cycle through keys in alleles dictionary.
	Note keys correspond to populations. See 'cline' function.
	3. Avalue: Index used to cycle through values in alleles dictionary. See 'cline function.
	4. pops: Dictionary containing information (e.g. allele frequencies) of each population
	Updated at the end of every generation.
	5. alleles: Dictionary used to stored lists of alleles population size for each population.
	Updated every generation.
	6. bot: Desired bottleneck proportion.
	7. Matrix: m x n dimensional matrix, initialized at the outset of each simulation.
	Contains empty cells that may become filled with populations.
	8. pop_counter: List with single element corresponding to running total of number
	of populations currently in existence.

	Returns:
	1. New population added to 'pops' dictionary and to Matrix
	2. Two allele lists ('A' and 'B') and population size ('S') added to 'alleles' dictionary
	'''
	if not alleles['1']['A']: #If there are no alleles for first population, pass. Only valid for first iteration when the populations have yet to be initialized
		#print 'There are no populations from which to sample!!'
		pass
	else:
		population = '10'
		number_of_items_to_pick = 1
		p_create = prob_create(max_K, max_p_create, Akey, Avalue)
		weights = [p_create, (1 - p_create)]
		create = [choice(population, weights) for i in range(number_of_items_to_pick)]
		if create[0] == '1': #If a '1' is sampled, create population
			pop_index = [(i, sublist.index(int(Akey))) for i, sublist in enumerate(Matrix) if int(Akey) in sublist][0]
			x, y = pop_index[0], pop_index[1]
			X, Y = (x_mat - 1), (y_mat - 1)
			# Create list containijng all neighboring cells
			Nlist = [(x2, y2) for x2 in range(x-1, x+2)
				for y2 in range(y-1, y+2)
					if (-1 < x <= X and
					-1 < y <= Y and
					(x != x2 or y != y2) and
					(0 <= x2 <= X) and
					(0 <= y2 <= Y))]
			# Reduced list containing only neighboring cells that lack a population
			Nlist_red = []
			for item in Nlist:
				i, j = item[0], item[1]
				if Matrix[i][j] == 0:
					Nlist_red.append(item)
			# If all neighboring cells are occupied, pass
			if not Nlist_red:
				pass
			else:
			# Otherwise, select and empty neighboring cell at random and place new population
				Nsam = random.randint(0, len(Nlist_red) - 1)
				i, j = Nlist_red[Nsam][0], Nlist_red[Nsam][1]
				pop_counter[0] += 1 #Increment 'pop_counter' by 1 if population is being created.
				#Add two lists to 'alleles' dictionary ('A' and 'B'). Naming: 'Pop.number'
				alleles['{0}'.format(pop_counter[0])] = {'A':sample_population_A(Avalue['A'], bottle(bot, Akey, Avalue)),'B':sample_population_B(Avalue['B'], bottle(bot, Akey, Avalue)), 'S':[bottle(bot, Akey, Avalue)]}
				pops['{0}'.format(pop_counter[0])] = [] #Empty list for new population. Naming same as alleles.
				Matrix[i][j] = pop_counter[0]
		else:
			pass

def phenotype(pA, pB):
	'''Calculates the frequency of recessive phenotype

	Modeled based on cyanogenesis system in white clover (T. repens).
	Recessive phenotype regers to phenotype lacking the production of
	hydrogen cyanide. With respect to genotypes, this corresponds to any
	genotype lacking a functional (i.e. dominant, 'A' or 'B') allele at
	either locus. Thus, only A- B- plants produce cyanide.

	Parameters:
	1. pA: Frequency of 'A' allele
	2. pB: Frequency of 'B' allele

	Returns:
	Frequency of recessive (i.e. acyanogenic) phenotype as float.
	'''
	qA = 1 - pA
	qB = 1 - pB
	mut= qA ** 2 + qB ** 2 - (qA ** 2 * qB ** 2)
	WT = 1 - mut
	return mut # Frequency of acyanogenic phenotype

# Cline function. Every generation, alleles are exchanged among populations. Populations follow
# logistic population growth. Ever generation, every population has some probability (p) of generating a new population, with alleles
# sampled from the population that created it.
def cline(locus_A, locus_B, steps, N, max_p_create, pops, alleles, bot, Matrix, max_K, K_dict, r, max_mig_rate, pop_counter, Distance_Dic, min_K):
	'''Generates cline based on number of generations.

	Loops over all populations for desired number of generations (i.e. steps).
	Every generation, calls necessary functions to determine new population size,
	allele frequency and created new populations. Does this for every population
	in existence. Finally, summarizes allele and phenotype frequencies of every
	populations in 'pops' dictionary.

	Parameters:
	1. locus_A: List containing alleles 'a' and 'A'
	2. locus_B: List containing alleles 'b' and 'B'
	3. steps: Number of generations
	4. N: Number of alleles to sample (i.e. population size). In this case, starting population size.
	5. max_p_create: Maximum desired probability of creating a new population
	6. pops: Dictionary containing information (e.g. allele frequencies) of each population
	Updated at the end of every generation.
	7. alleles: Dictionary used to stored lists of alleles population size for each population.
	Updated every generation.
	8. bot: Desired bottleneck proportion.
	9. Matrix: m x n dimensional matrix, initialized at the outset of each simulation.
	Contains empty cells that may become filled with populations.
	10. K: Desired carrying capacity.
	11. r: natural rate of increase
	12. max_mig_rate: Desired migration rate when distance = 0
	13. pop_counter: List with single element corresponding to running total of number
	of populations currently in existence.
	14. Distance_Dic: Dictionary containing distance and migration rates between populations.

	Returns:
	Summary stats (e.g. allele and phenotype frequency) appended to 'pops' dictionary.
	'''
	for i in range(steps):
		pop_list = pops.keys()
		for Akey, Avalue in alleles.items():
			if 'A' and 'B' in Avalue.keys():
				if not Avalue['A'] and not Avalue['B']:
					#If allele lists are empty, sample from list of initial allele frequencies. Only used for first generation
					Avalue['S'] = [N]
					Avalue['A'] = (sample_population_A(locus_A, N))
					Avalue['B'] = (sample_population_B(locus_B, N))
				else:
					#If allele lists are not empty, sample from previously sampled set of alleles.
						Avalue['S'] = pop_growth(r, Akey, Avalue, K_dict, Matrix)
						Avalue['A'] = sample_alleles_A(alleles_next_gen(Akey, pop_list, alleles, Matrix, Distance_Dic)[0], Akey, Avalue, alleles, r, K_dict, Matrix)
						Avalue['B'] = sample_alleles_B(alleles_next_gen(Akey, pop_list, alleles, Matrix, Distance_Dic)[1], Akey, Avalue, alleles, r, K_dict, Matrix)
			create_population(max_p_create, max_K, Akey, Avalue, pops, alleles, bot, Matrix, pop_counter) #Create population. Alleles will be sampled (see above). Population is currently empty list
		for Akey, Avalue in alleles.items():
			#Calculate allele and phenotype frequencies for every population, including newly created ones.
			pA = allele_freq(Avalue['A'])
			pB = allele_freq(Avalue['B'])
			pop_index = [(j, sublist.index(int(Akey))) for j, sublist in enumerate(Matrix) if int(Akey) in sublist][0]
			pops[Akey].append([pop_index[0], pop_index[1], Avalue['S'][0], i, round(pA, 3), round(pB, 3), round(phenotype(pA, pB), 3), max_mig_rate, K_dict[str(pop_index)][1], round(r, 3), max_p_create, bot, matrix_full(Matrix), max_K, min_K])
	return pops

def write_to_csv(writer, sim):
	'''Writes 'sim' dictionary to csv.

	Parameters:
	1. max_mig_rate: Desired migration rate when distance = 0
	2. bot: Desired bottleneck proportion.
	3. Dictionary containing key results (i.e. 'pops' dictionary) for every simulation.

	Returns:
	csv with results exported to directory of choice.
	'''
	for i in sim.keys():
		for j, x in sim[i].items():
			for z in x:
				writer.writerow([i, z[0], z[1], j, z[2], z[3], z[4], z[5], z[6], z[7], z[8], z[9], z[10], z[11], z[12], z[13], z[14]])

