#!/usr/bin/env python
# -*- coding: ASCII -*-
#

import os,sys
import argparse
import subprocess
import itertools
import csv
import tempfile

import numpy as np
from scipy import stats
from sklearn.decomposition import NMF
import matplotlib.pyplot as plt
import matplotlib.colors

import ef_cluster

separator = "\t"

# Get options
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
io = parser.add_argument_group('Input/output')
options = parser.add_argument_group('Method options')
plot = parser.add_argument_group('Plot options')

io.add_argument("-v","--vcf", dest="vcf", help="vcf file, readable by bcftools")
io.add_argument("-o","--output", dest="output_prefix", help="Output prefix", default="clusters")
io.add_argument("-d","--distances", dest="dist_mat", help="Pairwise distance matrix", default=None)
io.add_argument("-b", "--baps", dest="baps_file", help="BAPS clusters, for comparison", default=None)
options.add_argument("-m", "--mac", dest="min_mac", help="Minimum allele count",default=2, type=int)
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

# For all cluster sizes
divergences = []
for num_clusters in range(args.min_clusters, args.max_clusters + 1):
    # Run NMF
    sys.stderr.write("Running NMF on " + str(len(samples)) + " samples with " + str(num_clusters) + " clusters\n")
    model = NMF(n_components = num_clusters, init = 'nndsvd', alpha = args.alpha, l1_ratio = args.mixing, verbose = 0)
    decomposition = model.fit_transform(alignment)

    # assign to max val cluster, normalise each row by dividing by its sum
    clusters = np.argmax(decomposition, axis = 1)
    decomposition = decomposition/decomposition.sum(axis=1, keepdims = True)

    # evaluate cluster distances
    if os.path.isfile(args.dist_mat):
        ind_mat = np.zeros((len(samples), num_clusters))
        for cluster in range(0, num_clusters):
            ind_mat[clusters == cluster, cluster] = 1
        cluster_mat = np.dot(ind_mat, np.transpose(ind_mat))
        cluster_mat = cluster_mat/cluster_mat.sum()

        divergence = np.linalg.norm(cluster_mat - distances)
        sys.stderr.write("Divergence between clusters and distances is " + str(divergence) + "\n")
        divergences.append(divergence)

    # bin high entropy samples
    # note stats.entropy([1/num_clusters] * num_clusters = -log_2(1/N)
    if not args.max_entropy == None:
        sys.stderr.write("Max possible entropy " + "{0:.2f}".format(stats.entropy([1/num_clusters] * num_clusters))
                + "; binning over " + str(args.max_entropy) + "\n")
        sample_entropy = np.apply_along_axis(stats.entropy, 1, decomposition)
        clusters[sample_entropy > args.max_entropy] = -1

    # Write output (file for phandango)
    csv_out = open(args.output_prefix + "." + str(num_clusters) + "clusters.csv", 'w')
    csv_sep = ','

    csv_out.write(csv_sep.join(("name","cluster:o")) + "\n")
    for sample in range(0, clusters.size):
        csv_out.write(csv_sep.join((samples[sample], str(clusters[sample]))) + "\n")

    csv_out.close()

    # Draw structure plot
    if args.structure == True:
        # first sort by assigned cluster, find breaks
        sort_order = np.argsort(clusters)

        breaks = []
        prev = None
        for sample, sample_idx in zip(sort_order, range(0, clusters.size)):
            if not clusters[sample] == prev:
                breaks.append(sample_idx)
                prev = clusters[sample]

        bars = [] # In case I add a legend at some point
        ind = np.arange(len(samples))
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
            if line == breaks[1] and clusters[sort_order[0]] == -1: # bin cluster
                plt.axvline(line, linewidth=2, color='r', linestyle='dashed')
            elif not line == 0:
                plt.axvline(line, linewidth=2, color='w', linestyle='dashed')

        plt.title('Structure plot for %d clusters' % num_clusters)
        plt.ylabel('Assigned weight')
        plt.xlabel('Samples')
        plt.ylim([0,1])
        plt.xlim([0,len(samples)])
        plt.tick_params(
            axis='y',
            which='both',
            left='off',
            right='off',
            labelbottom='off')
        plt.savefig(args.output_prefix + "." + str(num_clusters) + "clusters.structure.pdf")
        plt.close()

# draw distances for each cluster
if len(divergences) > 1:
    best_clusters = args.min_clusters + divergences.index(min(divergences))
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
    plt.savefig(args.output_prefix + ".divergences.pdf")
    plt.close()

#TODO
# Compare with BAPS, if provided
if os.path.isfile(str(args.baps_file)):
    baps = {}
    baps_clusters = set()

    # Read in file
    with(open(args.baps_file)) as f:
        for _ in range(2):
            next(f)
        for line in f:
            line = line.rstrip()
            (row,cluster,lane) = line.split(separator)
            baps[lane] = int(cluster)
            baps_clusters.add(int(cluster))

    # Build array of arrays
    baps_sets = [[] for _ in range(len(baps_clusters))]
    for lane in baps.keys():
        baps_sets[baps[lane]-1].append(lane)

    # Compare with fast cluster results
    confusion_out = open(args.output_prefix + ".baps_confusion.txt", 'w')
    confusion_out.write(separator.join(["BAPS cluster","Total in BAPS cluster"] + [str(x) for x in sorted(set(found_clusters))]) + "\n")

    score = 0
    for cluster in baps_clusters:
        num_in_cluster = [0] * len(set(found_clusters))
        for lane in baps_sets[cluster-1]:
            fast_cluster = found_clusters[file_index[lane]]
            if not args.hier: # DBSCAN uses cluster = -1 for noise
                num_in_cluster[fast_cluster+1] += 1
            else: # others are 0-index based
                num_in_cluster[fast_cluster] += 1

        score += max(num_in_cluster)
        confusion_out.write(separator.join([str(cluster), str(len(baps_sets[cluster-1]))] + [str(x) for x in num_in_cluster]) + "\n")

    confusion_out.close()
    sys.stderr.write("Compared to " + str(len(baps_clusters)) + " BAPS clusters, " + str(score) + " samples of " + str(len(file_num)) + " are in the same clusters\n")

