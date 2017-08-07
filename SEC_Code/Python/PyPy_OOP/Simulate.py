import os
from datetime import datetime
import csv

from Cell import Cell
from Population import Population
from Parameters import *
from Functions import *

def simulate():

    print os.getpid()

    os.chdir(export_path)
    datestring = datetime.strftime(datetime.now(), '%Y%m%d')

    with open(datestring + "_SEC_Drift.Migration_(pA%.2f)(pB%.2f).csv" % (pA, pB), "wb") as f:
        writer = csv.writer(f)
        writer.writerow(["Sim","x","y","Generation","pA","pB","Cyan","p_create","K", "r", "bot", "Mig_rate", "Mat_full", "Pop_size", "min_K", "max_K"])

        print_progress(0, sims, prefix='', suffix='', decimals=1, bar_length=100)
        for s in range(sims):

            results = []

            locus_A = ( [ 'A' ] * int( N * pA ) ) + ( [ 'a' ] * int( round( N * qA ) ) )
            locus_B = ( [ 'B' ] * int( N * pB ) ) + ( [ 'b' ] * int( round( N * qB ) ) )

            Matrix = create_matrix(cols, rows)
            Matrix[0][0].pop = True
            Matrix[0][0].population = Population(N, locus_A, locus_B)

            cline(s, results, rows, cols, steps, pA, pB, Matrix, max_p_create, bot, min_K, max_K, r)

            write_to_csv(writer, results)

            del results

            print_progress(s + 1, sims, prefix='', suffix='', decimals=1, bar_length=100)


if __name__ == '__main__':
    simulate()