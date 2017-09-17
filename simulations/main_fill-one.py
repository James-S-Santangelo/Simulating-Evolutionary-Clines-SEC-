import os
from datetime import datetime
import csv
import sys

from simulations.cell import Cell
from simulations.population import Population
from simulations import functions


def simulate():
    """Main function. Runs 'sims' iterations of 'cline' across landscape

    Args:
        None

    Returns:
        None
    """

    sims = 1
    """int: Number of simulations"""

    steps = 50
    """int: Number of generations"""

    export_path = "/Users/jamessantangelo/Desktop/CSV/"
    """str: Path to where dataset should be exported"""

    print os.getpid()

    os.chdir(export_path)
    datestring = datetime.strftime(datetime.now(), '%Y%m%d')

    with open(datestring + "Drift.Migration_OneFill(pA%.2f)(pB%.2f).csv" % (Population.pA, Population.pB), "wb") as f:
        writer = csv.writer(f)
        writer.writerow(["Sim", "x", "y", "Generation", "pA", "pB", "Cyan",
                         "p_create", "K", "r", "bot", "Mig_rate", "Mat_full",
                         "Pop_size",
                         "min_K", "max_K"])

        functions.print_progress(0, sims, prefix='', suffix='', decimals=1, bar_length=100)
        for s in range(sims):

            results = []

            N = 1000
            pA = Population.pA
            pB = Population.pB
            qA = 1 - pA
            qB = 1 - pB
            locus_A = (['A'] * int(N * pA)) + (['a'] * int(round(N * qA)))
            locus_B = (['B'] * int(N * pB)) + (['b'] * int(round(N * qB)))

            Matrix = Cell.initialize_matrix()
            Matrix[0][0].pop = True
            Matrix[0][0].population = Population(N, locus_A, locus_B)

            functions.cline(s, results, steps, Matrix)

            functions.write_to_csv(writer, results)

            del results

            functions.print_progress(s + 1, sims, prefix='', suffix='', decimals=1, bar_length=100)


if __name__ == '__main__':
    simulate()
