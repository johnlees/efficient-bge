#!/usr/bin/env python
# -*- coding: ASCII -*-
#

import os,sys
from optparse import OptionParser
import subprocess

import numpy as np
from sklearn.manifold import TSNE
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt

import pdb

separator = "\t"

# Get options
parser = OptionParser()
parser.add_option("-a","--assembly_file", dest="assembly_file", help="Tab separated file with sample name and assembly location on each line")
parser.add_option("-o","--output", dest="output_prefix", help="Output prefix", default="clusters")
parser.add_option("-d","--dimensions",dest="dimensions",help="Number of t-SNE dimensions to embed to",type='int',default=2)
parser.add_option("-b", "--baps", dest="baps_file", help="BAPS clusters, for comparison",default=None)
parser.add_option("-m", "--mash", dest="mash_exec", help="Location of mash executable",default='mash')
(options, args) = parser.parse_args()

# Read in files
file_index = {}
file_name_index = {}
file_num = []
file_names = []
i = 0

with open(options.assembly_file) as f:
    for line in f:
        line = line.rstrip()
        (name, file_name) = line.split(separator)
        file_index[name] = i
        file_name_index[file_name] = i
        file_num.append(name)
        file_names.append(file_name)
        i += 1

# Run mash
sys.stderr.write("Calculating distances with mash\n")

distances = np.zeros((len(file_num), len(file_num)))
try:
    if not os.path.isfile("reference.msh"):
        mash_command = str(options.mash_exec) + " sketch -o reference " + separator.join(file_names)
        retcode = subprocess.call(mash_command, shell=True)
        if retcode < 0:
            print >>sys.stderr, "Mash failed with signal", -retcode

    pdb.set_trace()
    p = subprocess.Popen([str(options.mash_exec) + ' dist reference.msh reference.msh'], stdout=subprocess.PIPE, shell=True)
    for line in iter(p.stdout.readline, ''):
        line = line.decode('utf-8')
        line = line.rstrip()
        (name1, name2, dist, p, matches) = line.split(separator)
        if distances[file_name_index[name1], file_name_index[name2]] == 0:
            distances[file_name_index[name1], file_name_index[name2]] = dist

except OSError as e:
    print >>sys.stderr, "Mash Execution failed:", e

# Run t-SNE
sys.stderr.write("Embedding samples into " + str(options.dimensions) + " dimensions with t-SNE\n")
model = TSNE(n_components = options.dimensions, metric = "precomputed")
embedding = model.fit_transform(distances)

# Run DBSCAN
sys.stderr.write("Clustering samples with DBSCAN\n")
embedding = StandardScaler().fit_transform(embedding)
dbscan_clusters = DBSCAN().fit_predict(embedding)

# Write output (file for phandango)
sys.stderr.write("Found " + int(len(set(dbscan_clusters))) + " clusters\n")
csv_out = open(options.output_prefix, 'w')
csv_sep = ','

csv_out.write(csv_sep.join("name","cluster:o") + "\n")
for sample in dbscan_clusters:
    csv_out.write(csv_sep.join(file_num[sample], dbscan_clusters[sample]) +"\n")

# Draw plot (see http://scikit-learn.org/stable/auto_examples/cluster/plot_dbscan.html)
unique_labels = set(dbscan_clusters)
colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
for k, col in zip(unique_labels, colors):
    if k == -1:
        # Black used for noise.
        col = 'k'

    plt.plot(embedding, 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=10)

    xy = X[class_member_mask & ~core_samples_mask]
    plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=6)

plt.savefig(options.output_prefix + '.pdf')

# Compare with BAPS, if provided
if options.baps_file == True:
    baps = {}
    baps_clusters = set()

    # Read in file
    with(open(options.baps_file)) as f:
        for _ in range(2):
            next(f)
        for line in f:
            line = line.rstrip()
            (row,cluster,lane) = line.split(separator)
            baps[lane] = cluster
            baps_clusters.add(cluster)

    # Build array of arrays
    baps_sets = [[] for _ in range(len(baps_clusters))]
    for lane in baps.keys():
        baps_sets[baps[lane]].append(lane)

    # Compare with fast cluster results
    score = 0
    for cluster in baps_clusters:
        num_in_cluster = []
        for lane in baps_sets[cluster]:
            fast_cluster = dbscan_clusters[file_index[lane]]
            if fast_cluster == -1:
                num_in_cluster[0] += 1
            else:
                num_in_cluster[fast_cluster] += 1

        score += max(num_in_cluster)

    sys.stderr.write("Compared to " + str(len(baps_clusters)) + " BAPS clusters, " + str(score) + " samples are in the same clusters\n")


