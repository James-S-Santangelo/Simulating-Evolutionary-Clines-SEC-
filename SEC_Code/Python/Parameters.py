#!/usr/local/bin/python2.7

import math
import sys as sys

pA = 0.5 # Initial frequency of 'A' allele

pB = 0.5 # Initial frequency of 'B' allele

steps = 20 # Number of generations

N =  100 # Initial population size (i.e. starting size of first initialized population)

sims = 2 # Number of simulations

max_mig_rate = float(sys.argv[1]) # Maximum migration rate. Declines linearly with distance. See 'migration_rate' function in Functions.py

K = 100 # Carrying capacity

x_mat = 1 # Number of columns in matrix (i.e. landscape)

y_mat = 10 # Number of rows in matrix (i.e. landscape)

bot = float(sys.argv[2]) # Proportion of alleles sampled upon creation of new populations (i.e. bottleneck proportion)

max_p_create = float(1) # Maximum probability of creating a new population. Decreases with population size. See 'prob_create' function in Functions.py

r = float(math.log(2)) # Natural rate of increase. Enter desired instantaneous population growth rate (e.g. 2 is a doubling every generation).

qA = 1 - pA # Frequency of 'a' allele

qB = 1 - pB # Frequency of 'b' allele

export_path = sys.argv[3] # Path where final dataset will be exported
