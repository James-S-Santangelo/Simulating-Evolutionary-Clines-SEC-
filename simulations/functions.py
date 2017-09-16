import sys
import random
import bisect

from simulations.cell import Cell


def matrix_full(Matrix):
    """Checks whether every every cell in landscape matrix contains a population

    Args:
        Matrix (:obj:'list' of :obj:'int'): 2D array storing instance of Cell at every position.

    Returns:
        (int): 1 if Matrix is full, 0 if it still contains empty cells.
    """
    num_rows = len(Matrix)
    num_cols = len(Matrix[0])
    pop = [Matrix[i][j].pop for i in range(num_rows) for j in range(num_cols)]
    if all(pop):
        return 1
    else:
        return 0


def print_progress(iteration, total, prefix='', suffix='', decimals=1, bar_length=100):
    """Call in a loop to create terminal progress bar

    Args:
        iteration (int): current iteration
        total: (int): total iterations
        prefix (str): prefix string
        suffix (str): suffix string
        decimals (int): positive number of decimals in percent complete
        bar_length (int): character length of bar

    Returns:
        None: Prints progress bar to terminal
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
    """Calculates cumulative probability distribution from list of values

    Args:
        weights (:obj:'list' of :obj:'float'): list containing weighting factors. Must sum to 1.

    Returns:
        Result (:obj:'list' of :obj:'float'): list where each element is the probability of sampling the ith element.
    """
    total = sum(weights)
    result = []
    cumsum = 0
    for w in weights:
        cumsum += w
        result.append(cumsum / total)
    return result


def choice(possibilities, weights):
    """Weighted random sampling

    Args:
        possibilities (str): Possible choices to sample where each position of string is a different choice
        weights (:obj:'list' of :obj:'float'): list containing weighting factors. Must sum to 1.

    Returns:
        possibilities[idx] (str); Random choice from possibilities based on weights provided.
    """
    assert len(possibilities) == len(weights)
    cdf_vals = cdf(weights)
    x = random.random()
    idx = bisect.bisect(cdf_vals, x)
    return possibilities[idx]


def cline(s, results, steps, pA, pB, Matrix):
    """Does all of the heavy lifting in simulations

    Uses most core functions and methods from other modules to excercise both within- and between population dynamics. This includes population growth, migration, selection, and population creation.

    Args:
        s (int): current iteration of simulations
        results (:obj:'list' of :obj:'int' or 'float'): list storing all information of populations each generation and parameters used in simulations
        steps (int): number of generations
        pA (float): frequency of 'A' allele
        pB (float): frequency of 'B' allele
        Matrix (:obj:'list' of :obj:'int'): 2D array storing instance of Cell at every position.

    Returns:
        None: Simply generates results and appends them to results list.
    """
    for step in range(steps):

        pop_list = [(i, j) for j in range(Cell.num_cols) for i in range(Cell.num_rows) if Matrix[i][j].pop]

        mat_full = matrix_full(Matrix)

        for pop in pop_list:

            i, j = pop[0], pop[1]

            create_results(results, Matrix, i, j, mat_full, s, step)

            K = Matrix[i][j].real_K()

            Matrix[i][j].population.size = Matrix[i][j].population.pop_growth(K)

            pA1 = Matrix[i][j].alleles_next_gen(pop_list, Matrix)[0]
            pB1 = Matrix[i][j].alleles_next_gen(pop_list, Matrix)[1]

            Matrix[i][j].population.locus_A = Matrix[i][j].population.sample_alleles(pA1, 'Aa')
            Matrix[i][j].population.locus_B = Matrix[i][j].population.sample_alleles(pB1, 'Bb')

            Matrix[i][j].create_population(Matrix, K)


def create_results(results, Matrix, i, j, mat_full, s, step):
    """Appends all population statistics and global parameters to results list

        Args:
        results (:obj:'list' of :obj:'int' or 'float'): list storing all information of populations each generation and parameters used in simulations
        Matrix (:obj:'list' of :obj:'int'): 2D array storing instance of Cell at every position.
        i (int): Row number of cell in landscape matrix
        j (int): Column nujmber of cell in landscape matrix
        mat_full (int): 1 if matrix is full, 0 if there are sill empty cells.
        (float): frequency of 'A' allele
        s (int): current iteration of simulations
        step (int): current generation

    Returns:
        None: Appends them to results list.
    """
    r = Matrix[0][0].population.r
    max_create_prob = Cell.max_create_prob
    bot_prop = Cell.bot_prop
    max_mig_rate = Cell.max_mig_rate
    min_K = Cell.min_K
    max_K = Cell.max_K

    size = Matrix[i][j].population.size
    pA = Matrix[i][j].population.allele_freq(Matrix[i][j].population.locus_A)
    pB = Matrix[i][j].population.allele_freq(Matrix[i][j].population.locus_B)
    phen = Matrix[i][j].population.phenotype(pA, pB)
    K = Matrix[i][j].real_K()

    results.append([s, i, j, step, round(pA, 3), round(pB, 3),
                    round(phen, 3), max_create_prob, K, round(r, 3),
                    bot_prop, max_mig_rate, mat_full, size, min_K, max_K])


def write_to_csv(writer, results):
    """Writes results list to csv

    Aregs:
        writer (:obj:'class method'): writer method from csv module
        results (:obj:'list' of :obj:'int' or 'float'): list storing all information of populations each generation and parameters used in simulations

    Returns:
        None: Simply writes results to csv.
    """
    for z in results:

        writer.writerow([z[0], z[1], z[2], z[3], z[4], z[5], z[6], z[7], z[8], z[9], z[10], z[11], z[12], z[13], z[14], z[15]])
