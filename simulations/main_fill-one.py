import os
from datetime import datetime
import csv
import math

from simulations.cell import Cell
from simulations.population import Population
from simulations import functions

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

    with open(datestring + "Drift.Migration_OneFill(pA%.2f)(pB%.2f).csv" % (pA, pB), "wb") as f:
        writer = csv.writer(f)
        writer.writerow(["Sim", "x", "y", "Generation", "pA", "pB", "Cyan",
                         "p_create", "K", "r", "bot", "Mig_rate", "Mat_full",
                         "Pop_size",
                         "min_K", "max_K"])

        functions.print_progress(0, sims, prefix='', suffix='', decimals=1, bar_length=100)
        for s in range(sims):

            results = []

            locus_A = (['A'] * int(N * pA)) + (['a'] * int(round(N * qA)))
            locus_B = (['B'] * int(N * pB)) + (['b'] * int(round(N * qB)))

            Matrix = Cell.initialize_matrix(num_cols, num_rows)
            Matrix[0][0].pop = True
            Matrix[0][0].population = Population(N, locus_A, locus_B)

            functions.cline(s, results, num_rows, num_cols, steps, pA, pB, Matrix, max_p_create, bot_prop, min_K, max_K, max_mig_rate)

            functions.write_to_csv(writer, results)

            del results

            functions.print_progress(s + 1, sims, prefix='', suffix='', decimals=1, bar_length=100)


if __name__ == '__main__':
    simulate()
