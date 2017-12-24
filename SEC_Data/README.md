This directory contains the summary datasets created by running the executable R-scripts in ../dataset-creation. These are the datasets used in all analyses. The summary datasets can be found in subfolders that are named based the simulations and datasets that generated them, as follows:
* _allFill_Kvary: Fully colonized_. Drift and migration only. Varying carrying capacity across matrix, with largest population at rural end and smallest at urban end.
    * _allFill_Kvary_KminVary_: Varies the strength of the gradient in drift by varying the minimum carrying capacity at the urban end.
    * _allFill_Kvary_Kmin10_AllMig_: Strongest gradient in drift. Explores many migration rates fo in-depth understanding of the effects of migration.
* _allFill_Selection_: Fully colonized. Selection and migration only. Maximum strength of selection varied across matrix.
* _allFill_Kvary_AlleleFreq_: Fully colonized. Drift, migration and allele frequency variation.
* _allFill_Kvary_U-R_Selection_: Fully colonized, Drift, migration and selection. Carrying capacity varied, with largest population at urban end and smallest at rural end. Maximum strength of selection varied across matrix, running in opposite direction to gradient in drift.
* _oneFill_Bottlenecks_: Rural-most population at carrying capacity. Serial founder events. Drift and migration only.
* _oneFill_Bottlenecks_AlleleFreq_: Rural-most population at carrying capacity. Serial founder events. Drift, migration and allele frequency variation.
* _oneFill_Bottlenecks_Selection_: Rural-most population at carrying capacity. Drift, migration and selection. Selection in same direction as gradient in drift.
* _oneFill_Bottlenecks_U-R_Selection_: Urban-most population at carrying capacity. Serial founder events colonizing rural environment. Selection gradient opposite to drift gradient.

The above directories contain the following datasetes:
1. _RegSummary\*.csv_: Slope and P-value of a regression of the within-population frequency of HCN against distance from the urban-most population. Performed every generation for every simulation in all parameter combinations.
2. _MeansProps\*.csv_: Mean slope across 1000 simulations, proportion of significantly (_P_ < 0.05) positive (i.e. less HCN in urban environment) and negative (i.e. more HCN in urban environment) clines. Calculated for all parameter combinations.

_oneFill\*_ datasets contain an additional summary dataset:
1. _FreqFirstGen\*.csv_: Frequency of HCN in the first generation in which the population was founded. Useful for assessing the dynamics of HCN loss or fixation during serial founder events.
