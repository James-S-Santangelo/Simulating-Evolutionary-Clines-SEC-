import sys
import random
import bisect

import Cell
from Parameters import rows, cols, max_mig_rate

def matrix_full(Matrix):

	pop = [Matrix[i][j].pop for i in range(rows) for j in range(cols)]

	if all(pop):
		return 1
	else:
		return 0

def print_progress(iteration, total, prefix='', suffix='', decimals=1, bar_length=100):
	"""
	Call in a loop to create terminal progress bar
	@params:
		iteration   - Required  : current iteration (Int)
		total       - Required  : total iterations (Int)
		prefix      - Optional  : prefix string (Str)
		suffix      - Optional  : suffix string (Str)
		decimals    - Optional  : positive number of decimals in percent complete (Int)
		bar_length  - Optional  : character length of bar (Int)
	"""
	str_format = "{0:." + str(decimals) + "f}"
	percents = str_format.format(100 * (iteration / float(total)))
	filled_length = int(round(bar_length * iteration / float(total)))
	bar = '#' * filled_length + '-' * (bar_length - filled_length)

	sys.stdout.write('\r%s |%s| %s%s %s' % (prefix, bar, percents, '%', suffix)),

	if iteration == total:
		sys.stdout.write('\n')
	sys.stdout.flush()

def create_matrix(cols, rows):
	Matrix = [[0] * cols for i in range(rows)]
	for i in range(rows):
		for j in range(cols):
			Matrix[i][j] = Cell.Cell(i, j, pop = False)
	return Matrix

def cdf(weights):
	total = sum(weights)
	result = []
	cumsum = 0
	for w in weights:
		cumsum += w
		result.append(cumsum / total)
	return result

def choice(population, weights):
	assert len(population) == len(weights)
	cdf_vals = cdf(weights)
	x = random.random()
	idx = bisect.bisect(cdf_vals, x)
	return population[idx]

def cline(s, results, rows, cols, steps, pA, pB, Matrix, max_p_create, bot, min_K, max_K, r):

	for step in range(steps):

		pop_list = [(i, j) for j in range(cols) for i in range(rows) if Matrix[i][j].pop]
		mat_full = matrix_full(Matrix)

		for pop in pop_list:


			i, j = pop[0], pop[1]

			size = Matrix[i][j].population.size
			pA = Matrix[i][j].population.allele_freq(Matrix[i][j].population.locus_A)
			pB = Matrix[i][j].population.allele_freq(Matrix[i][j].population.locus_B)
			phen = Matrix[i][j].population.phenotype(pA, pB)
			K = Matrix[i][j].Distance_calc(0, 0, rows, cols, max_mig_rate, min_K, max_K)[1]

			results.append([s, i, j, step, round(pA, 3), round(pB, 3),
						round(phen, 3), max_p_create, K, round(r, 3),
						bot, max_mig_rate, mat_full, size, min_K, max_K])

			if step == 0:

				K = max_K
				pass

			else:

				K = Matrix[i][j].Distance_calc(0, 0, rows, cols, max_mig_rate, min_K, max_K)[1]
				Matrix[i][j].population.size = Matrix[i][j].population.pop_growth(Matrix[i][j].population.size, K, r)

				pA1 = Matrix[i][j].alleles_next_gen(rows, cols, pop_list, max_mig_rate, min_K, max_K)[0]
				pB1 = Matrix[i][j].alleles_next_gen(rows, cols, pop_list, max_mig_rate, min_K, max_K)[1]

				Matrix[i][j].population.locus_A = Matrix[i][j].population.sample_alleles(pA1, r, 'Aa')
				Matrix[i][j].population.locus_B = Matrix[i][j].population.sample_alleles(pB1, r, 'Bb')

			Matrix[i][j].create_population(rows, cols, Matrix, max_p_create, bot, K, r)


def write_to_csv(writer, results):

	for z in results:

		writer.writerow([z[0], z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8], z[9], z[10], z[11], z[12], z[13], z[14], z[15]])