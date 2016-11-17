#!/usr/bin/env python
# -*- coding: ASCII -*-
#

import os,sys
import argparse
import subprocess

import numpy as np
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt

# Remove when done
import pdb

separator = "\t"

# Get options
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("assembly_file", help="Tab separated file with sample name and assembly location on each line")
parser.add_argument("-o","--output", dest="output_prefix", help="Output prefix", default="clusters")
parser.add_argument("-d","--dimensions",dest="dimensions",help="Number of t-SNE dimensions to embed to",type=int,default=2)
parser.add_argument("--dist_mat",dest="distmat", help="Pre-computed distances.csv matrix", default=None)
parser.add_argument("--embedding",dest="embedding_file", help="Pre-computed t-SNE embedding", default=None)
parser.add_argument("-b", "--baps", dest="baps_file", help="BAPS clusters, for comparison",default=None)
parser.add_argument("-m", "--mash", dest="mash_exec", help="Location of mash executable",default='mash')
parser.add_argument("--min_pts", dest="min_pts", help="Minimum number of samples in each cluster",default=5, type=int)
parser.add_argument("--epsilon", dest="epsilon", help="Distance between DBSCAN clusters (pick with knn_plot)",default=0.1,type=float)
args = parser.parse_args()

# Read in files
file_index = {}
file_name_index = {}
file_num = []
file_names = []
i = 0

with open(args.assembly_file) as f:
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

if os.path.isfile(str(args.distmat)):
    distances = np.loadtxt(str(args.distmat), delimiter=",")
elif not os.path.isfile(str(args.embedding_file)):
    distances = np.zeros((len(file_num), len(file_num)))
    try:
        if not os.path.isfile("reference.msh"):
            mash_command = str(args.mash_exec) + " sketch -o reference " + separator.join(file_names)
            retcode = subprocess.call(mash_command, shell=True)
            if retcode < 0:
                print >>sys.stderr, "Mash failed with signal", -retcode

        p = subprocess.Popen([str(args.mash_exec) + ' dist reference.msh reference.msh'], stdout=subprocess.PIPE, shell=True)
        for line in iter(p.stdout.readline, ''):
            line = line.decode('utf-8')
            line = line.rstrip()
            if line != '':
                (name1, name2, dist, p, matches) = line.split(separator)
                if distances[file_name_index[name1], file_name_index[name2]] == 0:
                    distances[file_name_index[name1], file_name_index[name2]] = dist
            else:
                break

    except OSError as e:
        print >>sys.stderr, "Mash Execution failed:", e

    np.savetxt(str(args.output_prefix) + ".distances.csv", distances, delimiter=",")

# Run t-SNE
if os.path.isfile(str(args.embedding_file)):
    embedding = np.loadtxt(str(args.embedding_file), delimiter=",")
else:
    sys.stderr.write("Embedding samples into " + str(args.dimensions) + " dimensions with t-SNE\n")
    model = TSNE(n_components = args.dimensions, metric = "precomputed", n_iter=1000, verbose = 2, method = 'exact', learning_rate=1000, early_exaggeration = 4, n_iter_without_progress = 500, perplexity = 35)
    embedding = model.fit_transform(distances)

    np.savetxt(str(args.output_prefix) + ".embedding.csv", embedding, delimiter=",")

# Draw a k-distances plot
sys.stderr.write("Clustering samples with DBSCAN\n")
embedding = StandardScaler().fit_transform(embedding)
nbrs = NearestNeighbors(n_neighbors=args.min_pts, algorithm='kd_tree').fit(embedding)
nbr_distances = nbrs.kneighbors(embedding)

dist_list = nbr_distances[0].flatten()
dist_list.sort()

plt.plot(dist_list)
plt.title('k-distances plot')
plt.xlabel('Ordered distances')
plt.ylabel(str(args.min_pts) + '-NN distance')
plt.savefig(args.output_prefix + '.k_distances.pdf')
plt.close()

# Run DBSCAN
dbscan_clusters = DBSCAN(eps = args.epsilon, min_samples = args.min_pts).fit_predict(embedding)

# Write output (file for phandango)
sys.stderr.write("Found " + str(len(set(dbscan_clusters))) + " clusters\n")
csv_out = open(args.output_prefix + str(".csv"), 'w')
csv_sep = ','

csv_out.write(csv_sep.join(("name","cluster:o")) + "\n")
for sample in range(0, dbscan_clusters.size):
    csv_out.write(csv_sep.join((file_num[sample], str(dbscan_clusters[sample]))) + "\n")

csv_out.close()

# Draw plot (see http://scikit-learn.org/stable/auto_examples/cluster/plot_dbscan.html)
unique_labels = set(dbscan_clusters)
colours = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
for k, col in zip(unique_labels, colours):
    class_member_mask = (dbscan_clusters == k)
    xy = embedding[class_member_mask]

    if k == -1:
        # Black used for noise.
        col = 'k'
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=6)
    else:
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
             markeredgecolor='k', markersize=10)

plt.title('Estimated number of clusters: %d' % len(set(dbscan_clusters)))
plt.savefig(args.output_prefix + '.pdf')
plt.close()

# Compare with BAPS, if provided
if os.path.isfile(args.baps_file):
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
    confusion_out.write(separator.join(("BAPS cluster","# in BAPS cluster","# in fast cluster")) + "\n")

    score = 0
    for cluster in baps_clusters:
        num_in_cluster = [0] * len(set(dbscan_clusters))
        for lane in baps_sets[cluster-1]:
            fast_cluster = dbscan_clusters[file_index[lane]]
            if fast_cluster == -1:
                num_in_cluster[0] += 1
            else:
                num_in_cluster[fast_cluster] += 1

        score += max(num_in_cluster)
        confusion_out.write(separator.join((str(cluster), str(len(baps_sets[cluster-1])), str(max(num_in_cluster)))) + "\n")

    confusion_out.close()
    sys.stderr.write("Compared to " + str(len(baps_clusters)) + " BAPS clusters, " + str(score) + " samples of " + str(len(file_num)) + " are in the same clusters\n")

