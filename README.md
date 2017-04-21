# Simulating Evolutionary Clines (SEC)
### Author: James Santangelo

## Background

Phenotypic clines are the gradual change in the phenotype (or genotype) of a species over a geographical area and are often the result of species "tracking" some environmental factor. Given this association bewteen environmental variables and phenotype (or genotype) frequencies in populations, clines have often been interpreted as strong evidence of adaptive evolution. However, demographic (e.g. gene flow) and neutral (e.g. genetic drift) processes may also generate covariance between traits and environmental gradients (e.g. isolation by distance, spatially restricted gene flow). Ideally. neutral and demographic processes should be ruled out before invoking the role of natural selection in generating clines.

The formation of clines via neutral and demographic processes may be more likely in cases where the phenotype is the result of epistatic interactionsS among multiple genes. A useful model in this sense is the cyanogenesis polymorphism in white clover (*Trifolium repens*) whereby hydrogen cyanide (HCN) — a potent anti-herbivore defence — is controlled by two independently segregating Mendelian loci; plants require a functional (i.e. dominant) allele at each locus to produce HCN. Because of this, only 25% of genotypes (i.e. A- B-) produce HCN whereas all other genotypes (i.e. A- b-; a- B-; a- b-) do not, and both cyanogenic (HCN+) and ayanogenic (HCN-) plants co-occur in natural populations.

Over the past two years, collegues and I have collected approximately 12,000 *T. repens* plants from 723 populations spanning urban-rural gradients in the 16 largest cities in eastern North America. In each population we determined the proportion of plants producing HCN and examined how this varied across the urban-rural transect. In 8 of the 16 cities, we find significantly less HCN inside of cities than in surrounding rural regions. In other words there is a phenotypic cline whereby the frequency of HCN increases linearly with increasing distance from the urban core. While we have reson to believe selection plays some role in generating these cline, my goal here is to determine whether non-adaptive processes may play a part in their formation as well.

I am conducting a series of spatially explicit simulations in Python to examine the importance of both adaptive (i.e. selection) and non-adaptive (e.g. drift, gene-flow, colonization history) processes in generating phenotypic clines across space with an emphasis on multi-locus phenotypes.Simple simulations have shown that under a stepping-stone model of evolution — whereby a new population is created every generation by sampling from the previously generated population — the frequency of acyanogenesis can increase across space via genetic drift alone. This would be analogous to clover population starting in rural regions and gradually moving into the city, resulting in the observed patterns of less HCN in urban cores. While this model is useful as a starting point, it lacks biological realism as it does not allow for variations in population size, migration, colonization history, selection, or any other factors that may be important in structuring clines.

Ultimately, these simulations will provide me with a theoretical foundation upon which I can more clearly interpret the clinal patterns we are seeing in HCN across cities. I hope to address the following specific questions: (1) Can urban-rural clines in cyanogenesis be explained by neutral and demographic processes? (2) How strong does selection have to be to generate clines in the face of high migration, which is expected to homogenize HCN frequencies? (3) How does colonization history influence the formation of urban-rural clines in cyanogenesis? While I am modeling these simulations on the cyanogenesis system in clover, these results will be more broadly applicable to multi-locus phenotypic clines in other systems.


## How to use the code

The Python code in SEC_Code/Python/ contains multiple scripts: (1) Parameters.py, (2) Functions.py, (3) Simulate.py, and (4) Merge_csv.py. **Only Simulate.py is required to actually run the simulation as it will import the other scripts automatically**.

1. Open Parameters.py and set parameter values. Some of these take arguments from the command line (e.g. sys.argv), although this is not necessary
2. Run the simulate function by typing python Simulate.py into the command line. **Note the simulation use Python 2.7**
3. See dataset exported to specified path.

While the above approach is useful for generating small amounts of data and testing code functionality, it does not provide sufficient data to address the above questions. For this we need to run multiple iteration (e.g. 1000 simulations) of varying parameter combinations. This is a computationally intensive process that almost certainly require the use of an external server. To do this, I recommend the following approach be taken. Note this assumes you have access to a server with multiple cores of computing power.

1. Ensure you have GNU Parallel installed. Installation instructions can be found [here](https://www.gnu.org/software/parallel/).
2. Ensure the parameters you would like to vary accept arguments from the command line (e.g. using sys.argv[1], sys.argv[2], etc.).
3. Run the simulations in parallel using the following syntax
	`parallel python Simulations.py ::: x1 x2 x3 ::: y1 y2 y3 ::: export_path`
	where x1, x2, x3 refer to different values of the first parameter of interest and y1, y2, y3 refer to the sencond parameter of interest. More than 2 parameters can be varied. Hoever beware, the above syntax will perform all pariwise combinations of specified parameter values.
4. See all exported dataset in specified path (one for each run of the script).
5. Run the Merge_csv.py script to merge all of the dataset exported from above. Make sure to specify the path to the datasets from the command line.
