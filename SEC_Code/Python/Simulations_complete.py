#!/usr/bin/env python

# Modules used throughout script
import random
from collections import OrderedDict
import csv
import time
from datetime import datetime
import os
import itertools
import math
import numpy as np
import pandas as pd
from numpy.random import choice
import glob
import sys as sys

# Randomly sample 'N' alleles from lists containing alleles for locus A
def sample_population_A(locus_A, N):
	'''
	## PARAMETERS ##

	locus_A: List containing alleles 'a' and 'A'

	N: Number of alleles to sample from locus_A (i.e. population size)

	## USED IN FUNCTIONS ##

	1. migrate_A
	2. migrate_B
	3. create_population_2
	4. cline
	'''
	new_locus_A = [random.choice(locus_A) for _ in range(N)]
	return new_locus_A

# Randomly sample 'N' alleles from lists containing alleles for locus B
def sample_population_B(locus_B, N):
	'''
	## PARAMETERS ##

	locus_B: List containing alleles 'b' and 'B'

	N: Number of alleles to sample from locus_B (i.e. population size)

	## USED IN FUNCTIONS ##

	1. migrate_A
	2. migrate_B
	3. create_population_2
	4. cline
	'''
	new_locus_B = [random.choice(locus_B) for _ in range(N)]
	return new_locus_B

# Function used to determine migration rates based on distances between all populations.
# Linear decrease in migration rate with increasing distance. Rate of decreased based on slope ('m')
# which depends on the maximum distance between populations and the desired 'max_mig_rate' (i.e. migration rate when distace = 0)
def migration_rate(Distance_Dic, max_mig_rate):
	'''
	## PARAMETERS ##

	Distance_Dic: Dictionary containing distances between all pairwise combinations of cells in 'Matrix'.
	Will also contain migration rates after use of this function

	max_mig_rate: Desired migration rate when distance = 0

	## USED IN FUNCTIONS ##

	1. Distance_Mig
	'''
	max_dis = max(Distance_Dic.values()) # Max distance between populations
	m = (max_mig_rate - 0)/(max_dis[0] - 0) # Slope. Assumes close to no migration at max distance. Relized migration at max distance may be slightly greater than 0 due to rounding.
	for Dkey, Dvalue in Distance_Dic.items():
		Mig_prop = max_mig_rate - m*Dvalue[0]
		Distance_Dic[Dkey].append(round(Mig_prop, 3)) # Append migration rate to 'Distance_Dic' dictionary

# Function to calculate distances between all pairwise combinations of cells in matrix. Also calls
# the migration rate between populations and appends this to the distance dictionary.
def Distance_Mig(max_mig_rate):
	'''
	## PARAMETERS ##

	Matrix: m x n dimensional matrix, initialized at the outset of each simulation.
	Contains empty cells that may become filled with populations.

	max_mig_rate: Desired migration rate when distance = 0

	## USED IN FUNCTIONS ##

	1. Simulate
	'''
	Matrix = np.zeros((x_mat, y_mat), dtype = 'int')
	rows = [i for i in range(Matrix.shape[0])]
	cols = [i for i in range(Matrix.shape[1])]
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
	global Distance_Dic
	Distance_Dic = {'{0}.{1}'.format(key1, key2):[round(key3, 2)] for key1, key2, key3 in dis}
	migration_rate(Distance_Dic, max_mig_rate)

# Function used to sample alleles based on infinite allele pool model. Allele 'A'.
def sample_alleles_A(pA1, Akey, Avalue, alleles, r, K):
	'''
	## PARAMETERS ##

	pA1: Probability of sampling allele 'A'. Returned by 'alleles_next_gen()'

	Akey: Index used to cycle through keys in alleles dictionary.
	Note keys correspond to populations. See 'cline' function.

	Avalue: Index used to cycle through values in alleles dictionary. See 'cline function.

	alleles: Dictionary used to stored lists of alleles population size for each population.
	Updated every generation.

	r: natural rate of increase

	K: Carrying capacity (i.e. maximum sustainable population size)

	## USED IN FUNCTIONS ##

	1. alleles_next_gen
	'''
	list_of_candidates = ['A','a'] # Possible alleles to sample
	number_of_items_to_pick = pop_growth(r, Akey, Avalue, K) # Number to sample. Corresponds to next generation's size.
	probability_distribution= [pA1, (1 - pA1)] # Sampling probabilities. Returned by 'alleles_next_gen'
	draw = choice(list_of_candidates, number_of_items_to_pick, p = probability_distribution) # Sample alleles
	return list(draw) # Return list of newly sampled alleles. Becomes allele pool in the next generation

# Same as above but for allele 'B'
def sample_alleles_B(pB1, Akey, Avalue, alleles, r, K):
	list_of_candidates = ['B','b']
	number_of_items_to_pick = pop_growth(r, Akey, Avalue, K)
	probability_distribution = [pB1, (1 - pB1)]
	draw = choice(list_of_candidates, number_of_items_to_pick, p = probability_distribution)
	return list(draw)

def alleles_next_gen(Akey, pop_list, alleles, Matrix):
	'''
	## PARAMETERS ##

	Akey: Index used to cycle through keys in alleles dictionary.
	Note keys correspond to populations. See 'cline' function.

	pop_list: list containing all current populations in existence

	alleles: Dictionary used to stored lists of alleles population size for each population.
	Updated every generation.

	Matrix: m x n dimensional matrix, initialized at the outset of each simulation.
	Contains empty cells that may become filled with populations.

	## USED IN FUNCTIONS ##

	1. cline
	'''
	to_pop = (np.where(Matrix == int(Akey))[0][0], np.where(Matrix == int(Akey))[1][0]) # Location of focal population in matrix.
	migration_weighted = [] # List holding migration rates
	allele_weighted_A = [] # List holding frequency of 'A' alleles
	allele_weighted_B = [] # List holding frequency of 'B' alleles
	Size = [] # List holding population sizes
	for i in pop_list:
		if Akey == i:
			pass
		else:
			from_pop = (np.where(Matrix == int(i))[0][0], np.where(Matrix == int(i))[1][0]) # Location of source population in matrix
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
	return pA1, pB1

# From list containing alleles, calculate the frequency of 'A' or 'B' allele.
def allele_freq(locus):
	'''
	## PARAMETERS ##

	locus: List containing alleles

	## USED IN FUNCTIONS ##

	1. cline
	'''
	p = sum(1 * i.isupper() for i in locus)/float(len(locus))
	return p

# Return maximum rate of increase from final desired per capita growth rate. Assumes initial population size is 1
# and mating occurs over 1 generation
def rate_of_increase(Pf):
	return math.log(Pf)

# Function for logistic population growth. Takes current population size (from alleles dictionary)
# as input and return new population size.
def pop_growth(r, Akey, Avalue, K):
	'''
	## PARAMETERS ##

	r: natural rate of increase

	Akey: Index used to cycle through keys in alleles dictionary.
	Note keys correspond to populations. See 'cline' function.

	Avalue: Index used to cycle through values in alleles dictionary. See 'cline function.

	K: Carrying capacity (i.e. maximum sustainable population size)

	## USED IN FUNCTIONS ##

	1. cline
	'''
	size = Avalue['S'][0] # Retrieve size of population. 'Akey' allows indexing of alleles dictionary in cline function
	K = float(K)
	new_size = size * K/(size + (K - size) * math.exp(-r)) # Calculates the proportional reduction of population growth rate based on desired carrying capacity ('K'). At 'K', growth rate = 1 = no change
	return [int(round(new_size))]

# Simple bottleneck function.
def bottle(bot, Akey, Avalue):
	'''
	## PARAMETERS ##

	bot: Desired bottleneck proportion

	Akey: Index used to cycle through keys in alleles dictionary.
	Note keys correspond to populations. See 'cline' function.

	Avalue: Index used to cycle through values in alleles dictionary. See 'cline function.

	## USED IN FUNCTIONS ##

	1. create_population_2
	'''
	return int(math.ceil(bot * Avalue['S'][0]))

# Linear function that takes size of population and returns the probability of creating a new
# population.
def prob_create(K, max_p_create, Akey, Avalue):
	'''
	## PARAMETERS ##

	K: Carrying capacity

	max_p_create: maximum probability of creating a new population

	Akey: Index used to cycle through keys in alleles dictionary.
	Note keys correspond to populations. See 'cline' function.

	Avalue: Index used to cycle through values in alleles dictionary. See 'cline function.

	## USED IN FUNCTIONS ##

	1. create_population
	'''
	m = (max_p_create - 0)/(K - 0) # Slope. Assumes close to no migration at max distance. Relized migration at max distance may be slightly greater than 0 due to rounding.
	Size = Avalue['S'][0]
	p_create = Size * m
	return float(p_create)

# Create new population as empty list and add to 'pops' dictionary. Also create four lists of alleles
# sampled from pool of alleles from population that generated the new one. Alleles added to 'alleles'
# dictionary. Also adds population to matrix. This function first evaluates whether a population will
# be created then randomly selects a vacant neighboring cell where population will go. If no cells are
# vacant, the function passes.
def create_population(max_p_create, K, Akey, Avalue, pops, alleles, bot, Matrix):
	'''
	## PARAMETERS ##

	max_p_create: maximum probability of creating a new population

	Akey: Index used to cycle through keys in alleles dictionary.
	Note keys correspond to populations. See 'cline' function.

	Avalue: Index used to cycle through values in alleles dictionary. See 'cline function.

	pops: Dictionary containing information (e.g. allele frequencies) of each population
	Updated at the end of every generation.

	alleles: Dictionary used to stored lists of alleles population size for each population.
	Updated every generation.

	bot: Desired bottleneck proportion.

	Matrix: m x n dimensional matrix, initialized at the outset of each simulation.
	Contains empty cells that may become filled with populations.

	## USED IN FUNCTIONS ##

	1. cline
	'''
	if not alleles['1']['A']: #If there are no alleles for first population, pass. Only valid for first iteration when the populations have yet to be initialized
		#print 'There are no populations from which to sample!!'
		pass
	else:
		list_of_candidates = ['1','0']
		number_of_items_to_pick = 1
		p_create = prob_create(K, max_p_create, Akey, Avalue)
		probability_distribution = [p_create, (1 - p_create)]
		create = list(choice(list_of_candidates, number_of_items_to_pick, p = probability_distribution))
		if create[0] == '1': #If a '1' is sampled, create population
			x, y = np.where(Matrix == int(Akey))[0][0], np.where(Matrix == int(Akey))[1][0]
			X, Y = (Matrix.shape[0] - 1), (Matrix.shape[1] - 1)
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
				if Matrix[i, j] == 0:
					Nlist_red.append(item)
			# If all neighboring cells are occupied, pass
			if not Nlist_red:
				pass
			else:
			# Otherwise, select and empty neighboring cell at random and place new population
				Nsam = random.randint(0, len(Nlist_red) - 1)
				i, j = Nlist_red[Nsam][0], Nlist_red[Nsam][1]
				global pop_counter #global variable that tracks the number of populations created
				pop_counter += 1 #Increment 'pop_counter' by 1 if population is being created.
				#Add two lists to 'alleles' dictionary ('A' and 'B'). Naming: 'Pop.number'
				alleles['{0}'.format(pop_counter)] = {'A':sample_population_A(Avalue['A'], bottle(bot, Akey, Avalue)),'B':sample_population_B(Avalue['B'], bottle(bot, Akey, Avalue)), 'S':[bottle(bot, Akey, Avalue)]}
				pops['{0}'.format(pop_counter)] = [] #Empty list for new population. Naming same as alleles.
				Matrix[i, j] = pop_counter
		else:
			pass

# Given the frequencies of 'A' and 'B' alleles, return the frequency of the 'acyanogenic' phenotype (i.e. recessive
# at either the A locus, B locus, or both)
def phenotype(pA, pB):
	'''
	## PARAMETERS ##

	pA: Frequency of 'A' allele

	pB: Frequency of 'B' allele

	## USED IN FUNCTIONS ##

	1. cline
	'''
	qA = 1 - pA
	qB = 1 - pB
	mut= qA ** 2 + qB ** 2 - (qA ** 2 * qB ** 2)
	WT = 1 - mut
	return mut # Frequency of acyanogenic phenotype

# Cline function. Every generation, alleles are exchanged among populations. Populations follow
# logistic population growth. Ever generation, every population has some probability (p) of generating a new population, with alleles
# sampled from the population that created it.
def cline(locus_A, locus_B, steps, N, max_p_create, pops, alleles, bot, Matrix, K, r, max_mig_rate):
	'''
	## PARAMETERS ##

	locus_A: List containing alleles 'a' and 'A'

	locus_B: List containing alleles 'b' and 'B'

	steps: Number of generations

	N: Number of alleles to sample (i.e. population size). In this case, starting population size.

	p: Desired probability of creating a new population

	pops: Dictionary containing information (e.g. allele frequencies) of each population
	Updated at the end of every generation.

	alleles: Dictionary used to stored lists of alleles population size for each population.
	Updated every generation.

	bot: Desired bottleneck proportion.

	Matrix: m x n dimensional matrix, initialized at the outset of each simulation.
	Contains empty cells that may become filled with populations.

	K: Desired carrying capacity.

	r: natural rate of increase

	max_mig_rate: Desired migration rate when distance = 0

	## USED IN FUNCTIONS ##

	1. simulate
	'''
	for i in range(steps):
		pop_list = pops.keys()
		for Akey, Avalue in alleles.items():
			if Akey in pops.keys():
				if 'A' and 'B' in Avalue.keys():
					if not Avalue['A'] and not Avalue['B']:
						#If allele lists are empty, sample from list of initial allele frequencies. Only used for first generation
						Avalue['S'] = [N]
						Avalue['A'] = (sample_population_A(locus_A, N))
						Avalue['B'] = (sample_population_B(locus_B, N))
					else:
						#If allele lists are not empty, sample from previously sampled set of alleles.
						Avalue['S'] = pop_growth(r, Akey, Avalue, K)
						Avalue['A'] = sample_alleles_A(alleles_next_gen(Akey, pop_list, alleles, Matrix)[0], Akey, Avalue, alleles, r, K)
						Avalue['B'] = sample_alleles_B(alleles_next_gen(Akey, pop_list, alleles, Matrix)[1], Akey, Avalue, alleles, r, K)
				create_population(max_p_create, K, Akey, Avalue, pops, alleles, bot, Matrix) #Create population. Alleles will be sampled (see above). Population is currently empty list
		for Akey, Avalue in alleles.items():
			#Calculate allele and phenotype frequencies for every population, including newly created ones.
			pA = allele_freq(Avalue['A'])
			pB = allele_freq(Avalue['B'])
			pops[Akey].append([np.where(Matrix == int(Akey))[0][0], np.where(Matrix == int(Akey))[1][0], Avalue['S'][0], i, pA, pB, phenotype(pA, pB), max_mig_rate, K, r, max_p_create, bot])
	return pops


# Using the functions defined above, 'simulate' performs 'sims' iterations of the cline function -- simulating
# the combined effects of drift and migration in a spatially explicit framework  -- each time storing the results.
def simulate(pA, pB, steps, N, sims, max_mig_rate, K, bot, max_p_create, r, x_mat, y_mat):
	'''
	## PARAMETERS ##

	pA: Initial frequency of 'A' alleles.

	pB: Initial frequency of 'B' alleles.

	steps: Number of generations

	N: Number of alleles to sample (i.e. population size). In this case, starting population size.

	sims: Number of iterations.
	'''
	qA = 1 - pA # Frequency of 'a' allele
	qB = 1 - pB
	r = float(r)
	max_p_create = float(max_p_create)
	# Make the two lists based on the allele frequency to represent the initial population
	locus_A = (['A'] * int(N * pA) ) + (['a'] * int(round(N * qA)) ) # [A,A,A,A,a,a,a,a,....]
	locus_B = (['B'] * int(N * pB) ) + (['b'] * int(round(N * qB)) )
	####### sims simulations #####################
	# We will simulate 'steps' iterations of resampling this population to simulate drift
	# We will then repeat that simulation of 'steps' iterations 1000 times to get a mean
	##############################################
	for s in range(sims):
		pops = OrderedDict({'1':[]}) # Re-initialize dictionary to store populations
		alleles = OrderedDict({'1':{'A':[],'B':[],'S':[N]}}) # Re-initialize dictionary to store allele lists
		Matrix = np.zeros((x_mat, y_mat), dtype = 'int')
		Matrix[0, 0] = 1
		global pop_counter # Reset population counter
		pop_counter = 1
		# reset the population for each iteration. I don't actually think this is necessary
		locus_A = (['A'] * int(N * pA) ) + (['a'] * int(round(N * qA)))  # Re-initialize initial allele lists.
		locus_B = (['B'] * int(N * pB) ) + (['b'] * int(round(N * qB)))
		cline(locus_A,locus_B, steps, N, max_p_create, pops, alleles, bot, Matrix, K, r, max_mig_rate) # Run cline function
		global sim
		sim[s] = pops # Append results to global 'sim' dictionary

def write_to_csv(sim):
	datestring = datetime.strftime(datetime.now(), '%Y%m%d')
	DataFrame = []
	Colnames = ["Sim","x","y","Population","Pop_size","Generation","pA","pB","Acyan", "Mig_rate", "K", "r", "max_p_create", "bot"]
	for i in sim.keys():
		for j, x in sim[i].items():
			for z in x:
				DataFrame.append([i, z[0], z[1], j, z[2], z[3], z[4], z[5], z[6], z[7], z[8], z[9], z[10], z[11]])
	Test = pd.DataFrame(DataFrame, columns = Colnames)
	Test.to_csv(datestring + "_SEC_Drift.Migration.1D(m%.2f)(bot%.2f)_Complete.csv" % (max_mig_rate, bot))


# DEFINE LOCAL PARAMETERS

pA = 0.5 # Initial frequency of 'A' allele
pB = 0.5 # Initial frequency of 'B' allele
steps = 20 # Number of generations
N =  100 # Initial populations size (i.e. starting size of first initialized population)
sims = 2 # Number of simulations
max_mig_rate = float(sys.argv[1]) # Maximum migration rate. Declines linearly with distance. See 'migration_rate' function in Functions.py
K = 100 # Carrying capacity
x_mat = 1 # Number of columns in matrix (i.e. landscape)
y_mat = 20 # Number of rows in matrix (i.e. landscape)
bot = float(sys.argv[2]) # Proportion of alleles sampled upon creation of new populations (i.e. bottleneck proportion)
max_p_create = float(1) # Maximum probability of creating a new population. Decreases with population size. See 'prob_create' function in Functions.py
r = rate_of_increase(Pf = 2) # Natural rate of increase. 'Pf' is desired instantaneous population growth rate (e.g. Pf = 2 is a doubling every generation)
qA = 1 - pA # Frequency of 'a' allele
qB = 1 - pB # Frequency of 'b' allele

# DEFINE GLOBAL PARAMETER

pop_counter = 1 # Running count of number of populations created. Used for assingning populations to matrix

sim = {} # Disctionary to store results of simulations

Distance_Mig(max_mig_rate) # Returns a dictionary with all pairwise combinations of populations in matrix, their distance from one another, and the associated migration rate.
os.chdir(sys.argv[3])
simulate(pA, pB, steps, N, sims, max_mig_rate, K, bot, max_p_create, r, x_mat, y_mat)
write_to_csv(sim)

