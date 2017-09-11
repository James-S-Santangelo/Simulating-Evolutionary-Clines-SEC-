import os
from datetime import datetime
import csv
import math

from Cell import Cell
from Population import Population
import Functions

pA = 0.5
pB = 0.5
qA = 1 - pA
qB = 1 - pB
"""float: Frequency of 'A', 'B', 'a', and 'b' alleles"""

steps = 50
"""int: Number of generations"""

N = 1000
"""int: Size of starting population."""

sims = 1
"""int: Number of simulations"""

max_mig_rate = 0.0
"""float: Maximum migration rate between any two populations"""

max_K = 1000
min_K = 10
"""int: Maximum and minimum carying capacity of cells across the matrix"""

num_rows = 1
num_cols = 5
"""int: Number of rows and number of columns in landscape matrix"""

bot_prop = 1.0
"""float: Bottleneck proportion

Proportion of alleles sampled when new population is created
"""

max_p_create = 0
"""float: Maximum probability of creating a new population"""

r = math.log(1.5)
"""float: Instantaneous rate of population increase"""

export_path = "/Users/jamessantangelo/Desktop/CSV/"
"""str: Path to where dataset should be exported"""


def simulate():
    """Main function. Runs 'sims' iterations of 'cline' across landscape

    Args:
        None

    Returns:
        None
    """
    print os.getpid()

    os.chdir(export_path)
    datestring = datetime.strftime(datetime.now(), '%Y%m%d')

    with open(datestring + "Drift.Migration_AllFill(pA%.2f)(pB%.2f).csv" % (pA, pB), "wb") as f:
        writer = csv.writer(f)
        writer.writerow(["Sim", "x", "y", "Generation", "pA", "pB", "Cyan",
                         "p_create", "K", "r", "bot", "Mig_rate", "Mat_full",
                         "Pop_size",
                         "min_K", "max_K"])

        Functions.print_progress(0, sims, prefix='', suffix='', decimals=1, bar_length=100)
        for s in range(sims):

            results = []

            Matrix = Cell.initialize_matrix(num_cols, num_rows)
            for i in range(num_rows):
                for j in range(num_cols):
                    Matrix[i][j].pop = True
                    K = Matrix[i][j].real_K(num_rows, num_cols, min_K, max_K)
                    locus_A = (['A'] * int(K * pA)) + (['a'] * int(round(K * qA)))
                    locus_B = (['B'] * int(K * pB)) + (['b'] * int(round(K * qB)))
                    Matrix[i][j].population = Population(K, locus_A, locus_B)

            Functions.cline(s, results, num_rows, num_cols, steps, pA, pB, Matrix, max_p_create, bot_prop, min_K, max_K, r, max_mig_rate)

            Functions.write_to_csv(writer, results)

            del results

            Functions.print_progress(s + 1, sims, prefix='', suffix='', decimals=1, bar_length=100)


if __name__ == '__main__':
    simulate()