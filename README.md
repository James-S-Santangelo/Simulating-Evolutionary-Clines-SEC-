[![Build Status](https://travis-ci.org/James-S-Santangelo/Simulating-Evolutionary-Clines-SEC-.svg?branch=master)](https://travis-ci.org/James-S-Santangelo/Simulating-Evolutionary-Clines-SEC-)

**Example raw datasets and summary datasets can be found [HERE](https://www.dropbox.com/sh/nwqfc1pkfqg2hsn/AAA2vTDkLhjsSNoBVRfpZr8Ka?dl=0). Unzip the folders and read the README.md files for details on the datasets. All datasets will be publically available on DataDryad following publication.**

**A preprint of the manuscript, which has been accepted at _Proceeding of the Royal Society B_, is available on [bioRxiv](https://www.biorxiv.org/content/early/2018/03/27/289777).**

# Simulating Evolutionary Clines (SEC)
### Author: James Santangelo

## Background

Urban environments offer the opportunity to study the role of adaptive and non-adaptive evolutionary processes on an unprecedented scale. While the presence of parallel clines in heritable phenotypic traits is often considered strong evidence for the role of natural selection, non-adaptive evolutionary processes can also generate clines, and this may be more likely when traits have a non-additive genetic basis due to epistasis. In this paper, we use spatially-explicit simulations modelled according to the cyanogenesis (HCN) polymorphism in white clover (Trifolium repens) to examine the formation of phenotypic clines along urbanization gradients under varying levels of drift, gene flow and selection. HCN results from an epistatic interaction between two Mendelian-inherited loci. Our results demonstrate that the genetic architecture of this trait makes natural populations susceptible to decreases in HCN frequencies via drift. Gradients in the strength of drift across a landscape resulted in phenotypic clines with lower frequencies of HCN in strongly drifting populations, giving the misleading appearance of deterministic adaptive changes in the phenotype. Studies of heritable phenotypic change in urban populations should generate null models of phenotypic evolution based on the genetic architecture underlying focal traits prior to invoking selection’s role in generating adaptive differentiation.


## How to use the code

simulations/ contains all the code necessary to run the simulations. Simulation are run using [PyPy](https://pypy.org/) but can alternatively be run using python 2.7, although the runtime will be substantially longer. The simulations can be run as follows:

1. Clone the Github repository
2. Run `pip install -e .` This will install an editable version of the simulations that will auto-update when changes are made to the Github repository. Note this will install the simulations as a python package system-wide. To avoid this, install the package within a virtual environment.
3. The only scripts necessary for running the simulations are the _main_fill-\*.py_ scripts. The choice of which to use to run the simulations depends on the desired colonization scenario: _main_fill-all.py_ initializes all populations at carrying capacity whereas _main_fill-one.py_ initializes a single population at carrying capacity and includes colonization through serial founder effects.
    * Parameters are set in cell.py, population.py and _main_fill-\*.py_. Explanations of parameters are included in the relevant scripts.
    * Datasets will be exported in the directory specified in _main_fill-\*.py_ and named according to the filename specified in the same script.

While the above approach is useful for generating small amounts of data and testing code functionality, it does not provide sufficient data to address the above questions. For this we need to run multiple iterations (e.g. 1000 simulations) of varying parameter combinations. This is a computationally intensive process that requires the use of high performance computing cluster. To do this, I recommend the following approach be taken. Note this assumes you have access to a cluster.

1. Ensure you have GNU Parallel installed. Installation instructions can be found [here](https://www.gnu.org/software/parallel/).
2. Ensure the parameters you would like to vary accept arguments from the command line (e.g. using `sys.argv[1]`, `sys.argv[2]`, etc.). You can also modify the name of the exported datasets in _main_fill-\*.py_ to automatically name datasets according to the parameters being varied.
3. Run the simulations in parallel using the following syntax
	`parallel -u -v pypy Simulate_*.py ::: x1 x2 x3 ::: y1 y2 y3`
	where x1, x2, x3 refer to different values of the first parameter of interest and y1, y2, y3 refer to the sencond parameter of interest. More than 2 parameters can be varied. However beware, the above syntax will perform all pairiwise combinations of specified parameter values.
4. See all exported datasets in specified directory (one for each parameter combination).
