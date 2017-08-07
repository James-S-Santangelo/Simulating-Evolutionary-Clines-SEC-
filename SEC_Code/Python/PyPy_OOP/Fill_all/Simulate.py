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

	with open(datestring + "_SEC_Kvary_Migration_(m%.4f).csv" % (max_mig_rate), "wb") as f:
		writer = csv.writer(f)
		writer.writerow(["Sim","x","y","Generation","pA","pB","Cyan","K", "Mig_rate", "Pop_size", "min_K", "max_K"])

		print_progress(0, sims, prefix='', suffix='', decimals=1, bar_length=100)
		for s in range(sims):

			results = []

			Matrix = create_matrix(cols, rows)
			for i in range(rows):
				for j in range(cols):
					Matrix[i][j].pop = True
					K = Matrix[i][j].Distance_calc(0, 0, rows, cols, max_mig_rate, min_K, max_K)[1]
					locus_A = ( [ 'A' ] * int( K * pA ) ) + ( [ 'a' ] * int( round( K * qA ) ) )  # Re-initialize initial allele lists.
					locus_B = ( [ 'B' ] * int( K * pB ) ) + ( [ 'b' ] * int( round( K * qB ) ) )
					Matrix[i][j].population = Population(K, locus_A, locus_B)

			cline(s, results, rows, cols, steps, pA, pB, Matrix, max_p_create, bot, min_K, max_K, r)

			write_to_csv(writer, results)

			del results

			print_progress(s + 1, sims, prefix='', suffix='', decimals=1, bar_length=100)


if __name__ == '__main__':
	simulate()