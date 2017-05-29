## Short-term goals

1. ~~Create goals~~
2. ~~Add error caching to simulate function using `try`, `except`. Required for cases where rounding results in negative probability of sampling alleles. Currently, this results in termination of the script and loss of all simulations. Error caching will result in only the current iteration being terminated, allowing the simulations to continue following the error.~~ *I opted for a different solution here: I inserted a series of `if` statements into the function that calculates the sampling probabilities and forced the probabilities to 0 or 1 if they were outside this range due to rounding.*
3. ~~Consider exporting data to [Feather format](http://blog.cloudera.com/blog/2016/03/feather-a-fast-on-disk-format-for-data-frames-for-r-and-python-powered-by-apache-arrow/) instead of currently used `.csv` format. Feather [appears to provide faster read times into R](https://blog.dominodatalab.com/the-r-data-i-o-shootout/).~~ *While the feather format works as instructed, it is not currently suitable for long-term data storage, which is prefered for this project.*
4. ~~Add function to test if the matrix is full and add this as binary variable to dataframe. Will be used to test for transient dynamics and effects of bottlenecks (only present when matrix is not yet full).~~
5. Add function to calculate Fst, which will be used to estmiate realized migration rate to compare this to parameter value for migration rate.
6. Add selection to simulations and allow the strength of selection on different alleles to vary across the landscape.
7. Explore the combined effects of adaptive (i.e. selection) and non-adaptive (e.g. gene flow, drift) processes in generating multi-locus cline.
8. Explore the effects of colonization history (i.e. starting urban and colozing rural, vice-versa) on the formation of clines.

## Other goals

1. Upload R script that can be used to analyze data coming from simulations.
2. Generate figures showing example outputs of simulations.
3. Make code object-oriented and build a GUI to visualize simulations "in real time"


