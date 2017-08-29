import sys
import random
import bisect


def matrix_full(Matrix):

    num_rows = len(Matrix)
    num_cols = len(Matrix[0])

    pop = [Matrix[i][j].pop for i in range(num_rows) for j in range(num_cols)]

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


def cdf(weights):
    total = sum(weights)
    result = []
    cumsum = 0
    for w in weights:
        cumsum += w
        result.append(cumsum / total)
    return result


def choice(possibilities, weights):
    assert len(possibilities) == len(weights)
    cdf_vals = cdf(weights)
    x = random.random()
    idx = bisect.bisect(cdf_vals, x)
    return possibilities[idx]


def cline(s, results, num_rows, num_cols, steps, pA, pB, Matrix, max_create_prob, bot_prop, min_K, max_K, r, max_mig_rate):

    for step in range(steps):

        pop_list = [(i, j) for j in range(num_cols) for i in range(num_rows) if Matrix[i][j].pop]

        mat_full = matrix_full(Matrix)

        for pop in pop_list:

            i, j = pop[0], pop[1]

            create_results(results, Matrix, i, j, num_rows, num_cols, max_mig_rate, min_K, max_K, bot_prop, r, mat_full, s, step, max_create_prob)

            K = Matrix[i][j].real_K(num_rows, num_cols, min_K, max_K)

            Matrix[i][j].population.size = Matrix[i][j].population.pop_growth(K, r)

            pA1 = Matrix[i][j].alleles_next_gen(pop_list, Matrix, num_rows, num_cols, max_mig_rate)[0]
            pB1 = Matrix[i][j].alleles_next_gen(pop_list, Matrix, num_rows, num_cols, max_mig_rate)[1]

            Matrix[i][j].population.locus_A = Matrix[i][j].population.sample_alleles(pA1, 'Aa')
            Matrix[i][j].population.locus_B = Matrix[i][j].population.sample_alleles(pB1, 'Bb')

            Matrix[i][j].create_population(num_rows, num_cols, Matrix, max_create_prob, bot_prop, K)


def create_results(results, Matrix, i, j, num_rows, num_cols, max_mig_rate, min_K, max_K, bot_prop, r, mat_full, s, step, max_create_prob):

    size = Matrix[i][j].population.size
    pA = Matrix[i][j].population.allele_freq(Matrix[i][j].population.locus_A)
    pB = Matrix[i][j].population.allele_freq(Matrix[i][j].population.locus_B)
    phen = Matrix[i][j].population.phenotype(pA, pB)
    K = Matrix[i][j].real_K(num_rows, num_cols, min_K, max_K)

    results.append([s, i, j, step, round(pA, 3), round(pB, 3),
                    round(phen, 3), max_create_prob, K, round(r, 3),
                    bot_prop, max_mig_rate, mat_full, size, min_K, max_K])


def write_to_csv(writer, results):

    for z in results:

        writer.writerow([z[0], z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8], z[9], z[10], z[11], z[12], z[13], z[14], z[15]])
