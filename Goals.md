## Short-term goals

~~1. Create goals~~
2. Add error caching to simulate function using `try`, `except`. Required for cases where rounding results in negative probability of sampling alleles. Currently, this results in termination of the script and loss of all simulations. Error caching will result in only the current iteration being terminated, allowing the simulations to continue following the error.
3. Consider exporting data to [Feather format](http://blog.cloudera.com/blog/2016/03/feather-a-fast-on-disk-format-for-data-frames-for-r-and-python-powered-by-apache-arrow/) instead of currently used `.csv` format. Feather [appears to provide faster read times into R](https://blog.dominodatalab.com/the-r-data-i-o-shootout/).
4. Add selection to simulations and allow the strength of selection on different alleles to vary across the landscape.
5. Explore the combined effects of adaptive (i.e. selection) and non-adaptive (e.g. gene flow, drift) processes in generating multi-locus cline.
6. Explore the effects of colonization history (i.e. starting urban and colozing rural, vice-versa) on the formation of clines.

## Other goals

1. Upload R script that can be used to analyze data coming from simulations.
2. Generate figures showing example outputs of simulations.

## Long-term goals

1. Build a GUI for simulations that would allow parameter values to be input by user directly without having to hard-code them.

