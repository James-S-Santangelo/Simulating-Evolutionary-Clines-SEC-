import os
from datetime import datetime
import csv
import math

import Cell
import Functions
import Population

# Initial allele frequncies
pA = float(0.5)
pB = float(0.5)
qA = 1 - pA
qB = 1 - pB
# Number of generations
steps = 50
# Size of starting population
N = 100
# Number of simulations
sims = 5
# Maximum migration rate.
max_mig_rate = float(0.5)
# Carrying capacity
max_K = 1000
min_K = 10
# Number of rows and columns to be used in Matrix
rows = 1
cols = 50
# Proportion of alleles sampled upon creation of new populations
bot = float(0.5)
# Maximum probability of creating a new population.
max_p_create = float(0.1)
# Natural rate of increase.
r = float(math.log(1.5))
# Path where final dataset will be exported
export_path = "/Users/jamessantangelo/Desktop/CSV/"


def simulate():

    print os.getpid()

    os.chdir(export_path)
    datestring = datetime.strftime(datetime.now(), '%Y%m%d')

    with open(datestring + "_SEC_Drift.Migration_(pA%.2f)(pB%.2f).csv" % (pA, pB), "wb") as f:
        writer = csv.writer(f)
        writer.writerow(["Sim", "x", "y", "Generation", "pA", "pB", "Cyan",
                         "p_create", "K", "r", "bot", "Mig_rate", "Mat_full",
                         "Pop_size",
                         "min_K", "max_K"])

        Functions.print_progress(0, sims, prefix='', suffix='', decimals=1, bar_length=100)
        for s in range(sims):

            results = []

            locus_A = (['A'] * int(N * pA)) + (['a'] * int(round(N * qA)))
            locus_B = (['B'] * int(N * pB)) + (['b'] * int(round(N * qB)))

            Matrix = Cell.create_matrix(cols, rows)
            Matrix[0][0].pop = True
            Matrix[0][0].population = Population.Population(N, locus_A, locus_B)

            Functions.cline(s, results, rows, cols, steps, pA, pB, Matrix, max_p_create, bot, min_K, max_K, r, max_mig_rate)

            Functions.write_to_csv(writer, results)

            del results

            Functions.print_progress(s + 1, sims, prefix='', suffix='', decimals=1, bar_length=100)


if __name__ == '__main__':
    simulate()
