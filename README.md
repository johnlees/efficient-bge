# efficient-bge
scripts and methods for fast bacterial genomic epidemiology

## fast_cluster.py
Creates clusters of genetically related isolates, giving results comparable with [BAPS](http://www.helsinki.fi/bsg/software/BAPS/).

### Installation
You will need available on the command line

* [mash](http://mash.readthedocs.io/en/latest/)
* python3

You will also need the following python libraries
* [scikit-learn](http://scikit-learn.org/) > v0.18
* [numpy](http://www.numpy.org/)
* [matplotlib](http://matplotlib.org/)

these can be installed with

    pip install scikit-learn numpy matplotlib
    
### Usage
Briefly

    python fast_cluster.py assemblies.txt --epsilon 0.2 --baps figInfo.txt
    
Will run clustering, and compare it to the output from BAPS.

assemblies.txt should be a tab delimited file with the first column sample names, and the second column the FASTA file with assembled contigs.

The most important parameter to change is epsilon. After the first run, look for the point of inflexion on the y-axis of clusters.k_distances.pdf (e.g. 0.1) and run again with

    python fast_cluster.py assemblies.txt --embedding clusters.embedding.csv --epsilon 0.1 -b figInfo.txt
