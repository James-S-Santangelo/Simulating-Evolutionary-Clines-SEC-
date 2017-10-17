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

    sims = 1000
    """int: Number of simulations"""

    steps = 250
    """int: Number of generations"""

    export_path = "/scratch/research/projects/trifolium/SEC_Simulation.Evolutionary.Clines/SEC_Data/Drift.Migration/1D/AlleleFreq_Vary/all-fill/"
    """str: Path to where dataset should be exported"""

    print os.getpid()

    os.chdir(export_path)
    datestring = datetime.strftime(datetime.now(), '%Y%m%d')

    with open(datestring + "_" +"AlleleFreqVary_AllFill(pA%.2f)(pB%.2f).csv" % (Population.pA, Population.pB), "wb") as f:
        writer = csv.writer(f)
        writer.writerow(["Sim", "x", "y", "Generation", "pA", "pB", "Cyan",
                         "p_create", "K", "r", "bot", "Mig_rate", "Mat_full",
                         "Pop_size",
                         "min_K", "max_K"])

        functions.print_progress(0, sims, prefix='', suffix='', decimals=1, bar_length=100)
        for s in range(sims):

            results = []

            Matrix = Cell.initialize_matrix()
            pA = float(Population.pA)
            pB = float(Population.pB)
            qA = 1 - pA
            qB = 1 - pB
            for i in range(Cell.num_rows):
                for j in range(Cell.num_cols):
                    Matrix[i][j].pop = True
                    K = Matrix[i][j].real_K()
                    locus_A = (['A'] * int(K * pA)) + (['a'] * int(round(K * qA)))
                    locus_B = (['B'] * int(K * pB)) + (['b'] * int(round(K * qB)))
                    Matrix[i][j].population = Population(K, locus_A, locus_B)

            functions.cline(s, results, steps, Matrix)

            functions.write_to_csv(writer, results)

            del results

            functions.print_progress(s + 1, sims, prefix='', suffix='', decimals=1, bar_length=100)


if __name__ == '__main__':
    simulate()
