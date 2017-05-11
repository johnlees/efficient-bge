# efficient-bge
scripts and methods for faster bacterial genomic epidemiology. Currently
clustering

## Installation
You will need available on the command line

* [mash](http://mash.readthedocs.io/en/latest/)
* [bcftools](http://www.htslib.org/) >= v1.1
* python3

You will also need the following python libraries

* [scipy](https://www.scipy.org/scipylib/index.html)
* [scikit-learn](http://scikit-learn.org/) >= v0.18
* [numpy](http://www.numpy.org/)
* [matplotlib](http://matplotlib.org/)

these can be installed with

    pip install scipy scikit-learn numpy matplotlib

See <http://stackoverflow.com/questions/21784641/installation-issue-with-matplotlib-python> if you are having trouble with matplotlib.

## create_distances.py
Create a pair-wise distance matrix for use by either clustering
algorithm. For `fast_cluster.py` it will determine cluster membership and
number of clusters; for `slow_cluster.py` it will determine number of clusters only.

    python create_distances.py -v listeria_snps.vcf --assembly_file assemblies.txt

The vcf is needed to ensure consistency of order between the two files.
If you don't have one use `fast_cluster.py` directly.

Currently uses mash by default. This estimates distances using the
minHash distances with k-mers of length 21. Alternatively, a distance
matrix produced from `seaview` can be used.

Additional distance methods are being added. Please raise an issue if
there's one you'd like.

## fast_cluster.py
Creates clusters of genetically related isolates, giving results comparable with [BAPS](http://www.helsinki.fi/bsg/software/BAPS/). Uses t-SNE to embed the samples in 2D, then DBSCAN to cluster in this space. It can run `create_distances.py` for you, or you can use the result yourself using the `--dist_mat` option.

This is designed for very quick and dirty clustering when you don't have an alignment.

### Usage
Briefly

    python fast_cluster.py --assembly assemblies.txt --output clusters --epsilon 0.2 --baps figInfo.txt

Will run clustering, and compare it to the output from BAPS.

`assemblies.txt` should be a tab delimited file with the first column sample names, and the second column the FASTA file with assembled contigs.

The most important parameter to change is epsilon. After the first run, look for the point of inflexion on the y-axis of clusters.k_distances.pdf (e.g. 0.1) and run again with

    python fast_cluster.py assemblies.txt --embedding clusters.embedding.csv --epsilon 0.1 -b figInfo.txt

You can also use heirarchical clustering, which requires you to choose
a number of clusters

    python ../efficient-bge/fast_cluster.py --assembly assemblies.txt --output clusters --hier --clusters 22

### Output
The output .csv can be dragged-and-dropped into [phandango](http://phandango.net/) to display the clusters against a phylogeny.

The embedding and inferred clusters will be written as a pdf.

## slow_cluster.py
Slightly slower, but much more detailed cluster analysis. Uses non-negative
matrix factorisation to assign samples into a given number of clusters.
More fully featured than `fast_cluster.py`, and includes:

* Built in multithreading
* Hierarchical clustering (find clusters within clusters)
* Allows use of a bin cluster, for mixed ancestry samples

You need an alignment as input, but the clusters will be more accurate and you get much more information about cluster assignment.

### Usage
An example command line would be

    python slow_cluster.py -v snps.vcf -o clusters --max_clusters 10 --min_clusters 2 --mac 3 --distances clusters.npy --entropy 0.05 --hier 2 --threads 4 --write_all --structure

* The input is `snps.vcf`, which can be produced by [snp-sites](https://github.com/sanger-pathogens/snp-sites) from an alignment or directly from [bcftools](http://www.htslib.org/).
* `--min_clusters` and `--max_clusters` set the range of cluster numbers
  – all within this range will be tried.
* `--mac` sets the minimum minor allele count which contributes to the analysis. The default is 2.
* `--distances` is from `create_distances.py`, and is used to choose the optimal number of clusters.
* `--entropy` sets the level of mixing before a sample is put in the bin. The minimum is 0 (no mixing) and the reported maximum is for equal weighting from each cluster. Higher values bin fewer samples.
* `--hier` sets the levels of clustering to perform. For each level, the algorithm will be run recursively on the clusters found from the  previous level.
* `--threads` sets the number of CPUs to use.
* `--write_all` writes all clusterings, not just the optimal, to the output.
* `--structure` draws structure-like stacked bar plots of each clustering assignment, which are a useful extra diagnostic.

### Output
As with `fast_cluster.py`, a .csv file of the cluster assignments will beproduced.

`*.divergences.pdf` will show the fit to the distances of the cluster
assignments for each cluster number tried, and the best fit chosen.

`*.structure.pdf` show cluster weightings. Samples are along the x-axis
ordered by cluster assignment, and the stacked bar shows the weight of assignment to each cluster (each of
which is its own colour). Boundaries between cluster assignment are shown by vertical white dashed
lines. If a bin cluster was used this will be on the left,
separated by a vertical red dashed line.


