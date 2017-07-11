#!/usr/local/bin/python2.7

import math
import sys as sys

pA = float(0.5) # Initial frequency of 'A' allele

pB = float(0.5) # Initial frequency of 'B' allele

steps = 40 # Number of generations

N =  10 # Initial population size (i.e. starting size of first initialized population)

sims = 15 # Number of simulations

max_mig_rate = float(0.5) # Maximum migration rate. Declines linearly with distance. See 'migration_rate' function in Functions.py

# Carrying capacity
max_K = 1000
min_K = 1000 

x_mat = 1 # Number of columns in matrix (i.e. landscape)

y_mat = 10 # Number of rows in matrix (i.e. landscape)

bot = float(0.5) # Proportion of alleles sampled upon creation of new populations (i.e. bottleneck proportion)

max_p_create = float(1) # Maximum probability of creating a new population. Decreases with population size. See 'prob_create' function in Functions.py

r = float(math.log(2)) # Natural rate of increase. Enter desired instantaneous population growth rate (e.g. 2 is a doubling every generation).

qA = 1 - pA # Frequency of 'a' allele

qB = 1 - pB # Frequency of 'b' allele

export_path = "/Users/jamessantangelo/Desktop/CSV/" # Path where final dataset will be exported

