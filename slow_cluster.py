#!/usr/bin/env python
# -*- coding: ASCII -*-
#

# python libs
import os,sys
import argparse
import subprocess
import itertools
import csv
import tempfile

# multithreading
from multiprocessing.pool import ThreadPool
import threading
import time

# external libraries
import numpy as np
from scipy import stats
from sklearn.decomposition import NMF
import matplotlib.pyplot as plt
import matplotlib.colors

# internal libraries
import ef_cluster

separator = "\t"
sleep_time = 0.01

# Functions
def run_nmf(alignment, num_clusters, alpha, mixing):
    model = NMF(n_components = num_clusters, init = 'nndsvd', alpha = alpha, l1_ratio = mixing, verbose = 0)
    return(model.fit_transform(alignment))

def write_csv(file_name, clusters, samples):
    csv_out = open(file_name, 'w')
    csv_sep = ','

    csv_out.write(csv_sep.join(("name","cluster:o")) + "\n")
    for sample in range(0, len(clusters)):
        csv_out.write(csv_sep.join((samples[sample], str(clusters[sample]))) + "\n")

# Get options
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
io = parser.add_argument_group('Input/output')
efficiency = parser.add_argument_group('Efficiency options')
options = parser.add_argument_group('Method options')
plot = parser.add_argument_group('Plot options')

io.add_argument("-v","--vcf", dest="vcf", help="vcf file, readable by bcftools")
io.add_argument("-o","--output", dest="output_prefix", help="Output prefix", default="clusters")
io.add_argument("-d","--distances", dest="dist_mat", help="Pairwise distance matrix", default=None)
io.add_argument("--verbose", dest="verbose", help="Print verbose messages to stderr", action='store_true', default=False)
io.add_argument("--write_all", dest="write_all", help="Write all clusterings computed", action='store_true', default=False)
efficiency.add_argument("-t", "--threads", dest="threads", help="Number of threads to use", type=int, default=1)
options.add_argument("--hier", dest="hier", help="Depth of clustering to perform",default=1, type=int)
options.add_argument("--mac", dest="min_mac", help="Minimum allele count",default=2, type=int)
options.add_argument("--max_clusters", dest="max_clusters", help="Maximum number of clusters", default=10, type=int)
options.add_argument("--min_clusters", dest="min_clusters", help="Minimum number of clusters", default=2, type=int)
options.add_argument("--entropy", dest="max_entropy", help="Samples with an entropy higher than this will be assigned to a bin cluster", default=None, type=float)
options.add_argument("--regularisation", dest="alpha", help="Regularisation constant", default=0, type=float)
options.add_argument("--reg_mixing", dest="mixing", help="Regularisation L1/L2 mixing", default=0, type=float)
plot.add_argument("--structure", dest="structure", help="Create structure-like plots for all clusterings", action='store_true', default=False)
args = parser.parse_args()

# Read in header
samples = ef_cluster.read_header(args.vcf)

# Read in distances
if os.path.isfile(args.dist_mat):
    sys.stderr.write("Reading distance matrix\n")
    distances = np.load(args.dist_mat)
    # check size matches data (other checks in make_distance)
    if (distances.shape[0] != len(samples) or distances.shape[1] != len(samples)):
        sys.stderr.write("Sample size of distance matrix does not match vcf\n")
        sys.exit(1)

# Read in alignment
sys.stderr.write("Reading alignment\n")
tmp_csv = tempfile.NamedTemporaryFile()
bcftools_command = "bcftools norm -m - " + args.vcf + " | bcftools view -c " + str(args.min_mac) + ":minor - | bcftools query -f '[,%GT]\\n' - | sed 's/^,//' > " + tmp_csv.name
retcode = subprocess.call(bcftools_command, shell=True)
if retcode < 0:
    sys.stderr.write("bcftools conversion failed with code " + str(retcode) + "\n")
    tmp_csv.close()
    sys.exit(1)

alignment = np.loadtxt(tmp_csv.name, delimiter=',')
alignment = np.transpose(alignment)
tmp_csv.close()

best_clusters = 1
upper_clusters = [[0] * len(samples)]
if not args.max_entropy == None:
    upper_bin_clusters = upper_clusters[:] # with entropy bin

# Iterate through cluster levels requested
sys.stderr.write("Beginning clustering of " + str(len(samples)) + " samples\n")
for depth in range(0, args.hier):
    new_clusters = [0] * len(samples)
    new_bin_clusters = np.copy(new_clusters)
    for upper_cluster in set(upper_clusters[depth]):
        # Slice the alignment taking the cluster of the upper level, and enforcing the MAC cutoff again
        # This is the numpy version. Bit awkward as mixing lists (which deal in values) and numpy (which deals in indices)
        # sub_samples_idx = np.where(upper_clusters[depth] == upper_cluster)[0]
        # sub_samples = [samples[i] for i in sub_samples_idx]

        sub_samples_idx = [idx for idx in range(len(upper_clusters[depth])) if upper_clusters[depth][idx] == upper_cluster]
        sub_samples = [samples[i] for i in sub_samples_idx]
        if depth > 0:
            sub_alignment = alignment[sub_samples_idx,:]
            new_mac = np.sum(sub_alignment, axis = 0)
            sub_alignment = sub_alignment[:,(new_mac >= args.min_mac) & (new_mac <= len(sub_samples) - args.min_mac)]
        else:
            sub_alignment = alignment

        # For all cluster sizes, run NMF in parallel
        if args.verbose:
            sys.stderr.write("Running NMF at depth " + str(depth) + ", cluster " + str(upper_cluster) + " with " + str(args.max_clusters - args.min_clusters + 1) + " cluster sizes\n")
            sys.stderr.write("Sub-alignment size after frequency filtering " + str(sub_alignment.shape) + "\n")

        pool = ThreadPool(processes=args.threads)
        nmf_results = []
        for num_clusters in range(args.min_clusters, args.max_clusters + 1):
            nmf_results.append(pool.apply_async(run_nmf, (sub_alignment, num_clusters, args.alpha, args.mixing)))

        # Wait for NMF to complete
        while not(all(a_thread.ready() for a_thread in nmf_results)):
            time.sleep(sleep_time) # or 'pass' to not wait

        # Collect results
        divergences = []
        cluster_results = []
        cluster_bin_results = []
        for clustering_results, num_clusters in zip(nmf_results, range(args.min_clusters, args.max_clusters + 1)):
            decomposition = clustering_results.get()

            # assign to max val cluster, normalise each row by dividing by its sum
            clusters = np.argmax(decomposition, axis = 1)
            normalisation = decomposition.sum(axis=1, keepdims = True)
            normalisation = np.where(normalisation > 0, normalisation, 1) # Catch div/0 errors
            decomposition = decomposition/normalisation

            # evaluate cluster distances
            if os.path.isfile(args.dist_mat):
                ind_mat = np.zeros((len(sub_samples), num_clusters))
                for cluster in range(0, num_clusters):
                    ind_mat[clusters == cluster, cluster] = 1
                cluster_mat = np.dot(ind_mat, np.transpose(ind_mat))
                cluster_mat = cluster_mat/cluster_mat.sum()

                if depth > 0:
                    sub_distances = distances[sub_samples_idx, sub_samples_idx]
                else:
                    sub_distances = distances
                divergence = np.linalg.norm(cluster_mat - sub_distances)
                if args.verbose:
                    sys.stderr.write("Divergence between " + str(num_clusters) + " clusters and distances is " + str(divergence) + "\n")
                divergences.append(divergence)

            # give clusters unique values
            clusters += args.max_clusters * upper_cluster
            cluster_results.append(clusters)

            # bin high entropy samples
            # note stats.entropy([1/num_clusters] * num_clusters = -log_2(1/N)
            if not args.max_entropy == None:
                if args.verbose:
                    sys.stderr.write("Possible entropy range 0-" + "{0:.2f}".format(stats.entropy([1/num_clusters] * num_clusters))
                        + "; binning over " + str(args.max_entropy) + "\n")
                sample_entropy = np.apply_along_axis(stats.entropy, 1, decomposition)
                binned_clusters = np.copy(clusters)
                binned_clusters[sample_entropy > args.max_entropy] = -1
                cluster_bin_results.append(binned_clusters)

            # Draw structure plot
            output_prefix = args.output_prefix + ".level" + str(depth) + "_cluster" + str(upper_cluster) + "." + str(num_clusters) + "clusters"
            if args.structure == True:
                # first sort by assigned cluster, find breaks
                sort_order = np.argsort(binned_clusters)

                breaks = []
                prev = None
                for sample, sample_idx in zip(sort_order, range(0, binned_clusters.size)):
                    if not binned_clusters[sample] == prev:
                        breaks.append(sample_idx)
                        prev = binned_clusters[sample]

                bars = [] # In case I add a legend at some point
                ind = np.arange(len(sub_samples))
                colours = plt.cm.Spectral(np.linspace(0, 1, num_clusters))
                width = 1

                # draw stacked bar
                for cluster in range(0, num_clusters):
                    if cluster == 0:
                        bars.append(plt.bar(ind, decomposition[sort_order,cluster], width, color=colours[cluster],))
                        bottoms = decomposition[sort_order,cluster]
                    else:
                        bars.append(plt.bar(ind, decomposition[sort_order,cluster], width, color=colours[cluster],
                                bottom=bottoms))
                        bottoms += decomposition[sort_order,cluster]

                # draw breaks between clusters
                for line in breaks:
                    if line == breaks[1] and binned_clusters[sort_order[0]] == -1: # bin cluster
                        plt.axvline(line, linewidth=2, color='r', linestyle='dashed')
                    elif not line == 0:
                        plt.axvline(line, linewidth=2, color='w', linestyle='dashed')

                plt.title('Structure plot for %d clusters' % num_clusters)
                plt.ylabel('Assigned weight')
                plt.xlabel('Samples')
                plt.ylim([0,1])
                plt.xlim([0,len(sub_samples)])
                plt.tick_params(
                    axis='y',
                    which='both',
                    left='off',
                    right='off',
                    labelbottom='off')
                plt.savefig(output_prefix + ".structure.pdf")
                plt.close()

            # Write output (file for phandango)
            if args.write_all:
                write_csv(output_prefix + ".csv", clusters, sub_samples)

        # choose best clustering
        best_clusters = args.min_clusters + divergences.index(min(divergences))
        for sample_idx, cluster_idx in zip(sub_samples_idx, range(0, len(sub_samples_idx))):
            new_clusters[sample_idx] = cluster_results[best_clusters - args.min_clusters][cluster_idx]
            new_bin_clusters[sample_idx] = cluster_bin_results[best_clusters - args.min_clusters][cluster_idx]

        # draw distances for each cluster
        output_prefix = args.output_prefix + ".level" + str(depth) + "_cluster" + str(upper_cluster)
        if len(divergences) > 1:
            if args.verbose:
                sys.stderr.write("Minimum divergence " + str(min(divergences)) + " at " + str(best_clusters) + " clusters\n")

            colours = ['red', 'blue']
            levels = [0, 1]
            min_colours = np.where(divergences == min(divergences), 0, 1)
            cmap, norm = matplotlib.colors.from_levels_and_colors(levels=levels, colors=colours, extend='max')

            cluster_vals = np.arange(args.min_clusters, args.max_clusters + 1)
            plt.plot(cluster_vals, divergences, 'k')
            plt.scatter(cluster_vals, divergences, s=40, c=min_colours, cmap=cmap, norm=norm)

            plt.title('Lowest divergence at %d clusters' % best_clusters)
            plt.xlabel('Number of clusters')
            plt.ylabel('Divergence')
            plt.xticks(cluster_vals, cluster_vals)
            plt.savefig(output_prefix + ".divergences.pdf")
            plt.close()

    upper_clusters.append(new_clusters)
    upper_bin_clusters.append(new_bin_clusters)
    output_file = args.output_prefix + ".level" + str(depth) + "_best_clusters.csv"
    if args.max_entropy == None:
        write_csv(output_file, new_clusters, samples)
    else:
        write_csv(output_file, new_bin_clusters, samples)


# Done
sys.stderr.write("Done\n")

