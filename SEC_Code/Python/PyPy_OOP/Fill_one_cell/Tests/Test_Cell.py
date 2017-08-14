from Fill_one_cell import Cell
from Fill_one_cell import Population

rows = 1
cols = 10
N = 10
pA = 0.5
pB = 0.5
qA = 1 - pA
qB = 1 - pB
locus_A = (['A'] * int(N * pA)) + (['a'] * int(round(N * qA)))
locus_B = (['B'] * int(N * pB)) + (['b'] * int(round(N * qB)))

Matrix = Cell.create_matrix(cols, rows)

Matrix[0][0].pop = True
Matrix[0][0].population = Population.Population(N, locus_A, locus_B)


def test_init():
    assert Matrix[0][0].i == 0


def test_prob_create():
    assert Matrix[0][0].prob_create(1000, 0) == 0

    assert Matrix[0][0].prob_create(1000, 0.5) == 0.005
    assert Matrix[0][0].prob_create(10, 0.5) == 0.5

    assert Matrix[0][0].prob_create(1000, 1) == 0.01


def test_create_population():
    Matrix[0][0].create_population(rows, cols, Matrix, 1.0, 1.0, 10, 1.5)
    assert Matrix[0][1].pop
    assert Matrix[0][1].population.size == 10
