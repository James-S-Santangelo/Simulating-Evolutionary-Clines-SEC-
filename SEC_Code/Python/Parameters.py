#!/usr/local/bin/python2.7

import math
import sys as sys

pA = float(sys.argv[1]) # Initial frequency of 'A' allele

pB = float(sys.argv[2]) # Initial frequency of 'B' allele

steps = 30 # Number of generations

N =  10 # Initial population size (i.e. starting size of first initialized population)

sims = 5 # Number of simulations

<<<<<<< HEAD
max_mig_rate = float(0.0) # Maximum migration rate. Declines linearly with distance. See 'migration_rate' function in Functions.py
=======
max_mig_rate = float(0.5) # Maximum migration rate. Declines linearly with distance. See 'migration_rate' function in Functions.py
>>>>>>> Removed all numpy functionality. Added custom functions for weighted probability sampling

K = 1000 # Carrying capacity

x_mat = 1 # Number of columns in matrix (i.e. landscape)

y_mat = 10 # Number of rows in matrix (i.e. landscape)

<<<<<<< HEAD
bot = float(1.0) # Proportion of alleles sampled upon creation of new populations (i.e. bottleneck proportion)
=======
bot = float(0.5) # Proportion of alleles sampled upon creation of new populations (i.e. bottleneck proportion)
>>>>>>> Removed all numpy functionality. Added custom functions for weighted probability sampling

max_p_create = float(1) # Maximum probability of creating a new population. Decreases with population size. See 'prob_create' function in Functions.py

r = float(math.log(2)) # Natural rate of increase. Enter desired instantaneous population growth rate (e.g. 2 is a doubling every generation).

qA = 1 - pA # Frequency of 'a' allele

qB = 1 - pB # Frequency of 'b' allele

export_path = "/Users/jamessantangelo/Documents/Academia/Doctorate_PhD/Projects/SEC_Simulating.evolutionary.clines/SEC_Git/SEC_Code/Python" # Path where final dataset will be exported

