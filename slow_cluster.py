#!/usr/bin/env python
# -*- coding: ASCII -*-
#

import os,sys
import argparse
import subprocess
import itertools
import csv

import numpy as np
from sklearn.decomposition import NMF
from sklearn.cluster import DBSCAN,AgglomerativeClustering
import matplotlib.pyplot as plt
import matplotlib.markers

separator = "\t"

# Get options
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
io = parser.add_argument_group('Input/output')
options = parser.add_argument_group('Method options')

# TODO add in sample names
io.add_argument("--alignment", dest="alignment", help="csv of alignment (MxN), produced by bcftools")
io.add_argument("--samples", dest="sample_names", help="sample names, in the same order as the alignment")
io.add_argument("-o","--output", dest="output_prefix", help="Output prefix", default="clusters")
io.add_argument("-b", "--baps", dest="baps_file", help="BAPS clusters, for comparison", default=None)
options.add_argument("-m", "--mac", dest="min_mac", help="Minimum allele count",default=2, type=int)
options.add_argument("--max_clusters", dest="max_clusters", help="Maximum number of clusters", default=100, type=int)
args = parser.parse_args()

# Read in files
#file_index = {}
#file_name_index = {}
#file_num = []
#file_names = []
#i = 0
#
#if os.path.isfile(str(args.assembly_file)):
#    with open(str(args.assembly_file)) as f:
#        for line in f:
#            line = line.rstrip()
#            (name, file_name) = line.split(separator)
#            file_index[name] = i
#            file_name_index[file_name] = i
#            file_num.append(name)
#            file_names.append(file_name)
#            i += 1

samples = []
with open(str(args.sample_names)) as f:
    for line in f:
        line = line.rstrip()
        samples.append(line)

# Read in alignment
sys.stderr.write("Reading alignment\n")
alignment = np.loadtxt(args.alignment, delimiter=',')
alignment = np.transpose(alignment)

# Run NMF
num_clusters = 2
sys.stderr.write("Running NMF on " + str(num_clusters) + " clusters\n")
model = NMF(n_components = num_clusters, verbose = 3)
decomposition = model.fit_transform(alignment)

# normalise each row by dividing by its sum
# TODO assign to max val, unless close to others in which case bin
clusters = np.argmax(decomposition, axis = 1)
decomposition = decomposition/decomposition.sum(axis=1, keepdims = True)

print(clusters)
print(decomposition)

# Write output (file for phandango)
csv_out = open(args.output_prefix + str(".csv"), 'w')
csv_sep = ','

csv_out.write(csv_sep.join(("name","cluster:o")) + "\n")
for sample in range(0, clusters.size):
    csv_out.write(csv_sep.join((samples[sample], str(clusters[sample]))) + "\n")

csv_out.close()

#TODO
# Draw stacked bar plot
#unique_labels = set(found_clusters)
#colours = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
#markers = itertools.cycle(matplotlib.markers.MarkerStyle.filled_markers)
#for k, col in zip(unique_labels, colours):
#    class_member_mask = (found_clusters == k)
#    xy = embedding[class_member_mask]
#
#    if k == -1:
#        # Black used for noise.
#        col = 'k'
#        plt.plot(xy[:, 0], xy[:, 1], marker='o', markerfacecolor=col,
#             markeredgecolor='k', markersize=6, linestyle='None')
#    else:
#        plt.plot(xy[:, 0], xy[:, 1], marker=next(markers), markerfacecolor=col,
#             markeredgecolor='k', markersize=10, linestyle='None')
#
#plt.title('Estimated number of clusters: %d' % len(set(found_clusters)))
#plt.savefig(args.output_prefix + '.pdf')
#plt.close()

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

