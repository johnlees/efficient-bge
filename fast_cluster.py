#!/usr/bin/env python
# -*- coding: ASCII -*-
#

import os,sys
import argparse
import subprocess
import itertools
import csv

import numpy as np
from sklearn.manifold import TSNE,MDS
from sklearn.decomposition import NMF
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import DBSCAN,AgglomerativeClustering
import matplotlib.pyplot as plt
import matplotlib.markers

import ef_cluster

separator = "\t"

# Get options
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
io = parser.add_argument_group('Input/output')
methods = parser.add_argument_group('Method')
options = parser.add_argument_group('Method options')

io.add_argument("--assembly_file", help="Tab separated file with sample name and assembly location on each line")
io.add_argument("-o","--output", dest="output_prefix", help="Output prefix", default="clusters")
io.add_argument("--dist_mat",dest="distmat", help="Pre-computed distances.csv.npy matrix", default=None)
io.add_argument("--seaview_mat",dest="seaview_mat", help="Pre-computed distances matrix from seaview", default=None)
io.add_argument("--embedding",dest="embedding_file", help="Pre-computed t-SNE embedding", default=None)
io.add_argument("-b", "--baps", dest="baps_file", help="BAPS clusters, for comparison", default=None)
methods.add_argument("--mds", dest="mds", action='store_true', default=False, help="Use MDS instead of t-SNE")
methods.add_argument("--hier", dest="hier", action='store_true', default=False, help="Use heirarchical clustering instead of embedding")
options.add_argument("-d","--dimensions",dest="dimensions",help="Number of t-SNE dimensions to embed to",type=int,default=2)
options.add_argument("-m", "--mash", dest="mash_exec", help="Location of mash executable",default='mash')
options.add_argument("--kmer_size", dest="kmer_size", help="K-mer size for mash sketches", default="21")
options.add_argument("--sketch_size", dest="sketch_size", help="Size of mash sketch", default="1000")
options.add_argument("--min_pts", dest="min_pts", help="Minimum number of samples in each cluster",default=5, type=int)
options.add_argument("--epsilon", dest="epsilon", help="Distance between DBSCAN clusters (pick with knn_plot)",default=0.1,type=float)
options.add_argument("--clusters", dest="clusters", help="Number of clusters to form if using hierarchical clustering", default=20, type=int)
args = parser.parse_args()

# Need one and only one of these inputs (XOR)
if not bool(os.path.isfile(str(args.assembly_file))) ^ bool(os.path.isfile(str(args.seaview_mat))):
    sys.stderr.write("Need one of --assembly_file or --seaview_mat\n")
    sys.exit(1)

# Read in files
file_index = {}
file_name_index = {}
file_num = []
file_names = []
i = 0

if os.path.isfile(str(args.assembly_file)):
    with open(str(args.assembly_file)) as f:
        for line in f:
            line = line.rstrip()
            (name, file_name) = line.split(separator)
            file_index[name] = i
            file_name_index[file_name] = i
            file_num.append(name)
            file_names.append(file_name)
            i += 1

# Create distance matrix
if os.path.isfile(str(args.distmat)):
    distances = np.load(str(args.distmat))
elif os.path.isfile(str(args.seaview_mat)):
    # Read from seaview
    with open(str(args.seaview_mat), 'r') as f:
        reader = csv.reader(f)
        seaview_mat = list(reader)

    # Create file objects as for assemblies (may be used in BAPS comparison)
    file_num = seaview_mat[0][1:len(seaview_mat[0])]
    for file_name in file_num:
        file_index[file_name] = i

    # Distances remove the row and column names
    distances = np.zeros((len(file_num), len(file_num)))
    line_nr = 0
    for line in seaview_mat:
        if line_nr > 0:
            distances[line_nr-1,:] = np.asarray(line[1:len(seaview_mat[0])])
        line_nr += 1

elif not os.path.isfile(str(args.embedding_file)) and not args.hier:
    # Run mash
    sys.stderr.write("Calculating distances with mash\n")
    distances = np.zeros((len(file_num), len(file_num)))
    try:
        if not os.path.isfile("reference.msh"):
            ef_cluster.run_mash_sketch(args.mash_exec, file_names, args.kmer_size, args.sketch_size)

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
        sys.stderr.write("Mash Execution failed: " + str(e))

    np.save(str(args.output_prefix) + ".distances.csv", distances)

# Run embedding
if os.path.isfile(str(args.embedding_file)):
    embedding = np.loadtxt(str(args.embedding_file), delimiter=",")
else:
    # metric MDS
    if args.mds:
        sys.stderr.write("Embedding samples into " + str(args.dimensions) + " dimensions with MDS\n")
        model = MDS(n_components = args.dimensions, metric=True, dissimilarity = "precomputed", verbose = 1)
    # default is t-SNE
    else:
        sys.stderr.write("Embedding samples into " + str(args.dimensions) + " dimensions with t-SNE\n")
        #model = TSNE(n_components = args.dimensions, metric = "precomputed", n_iter=1000, verbose = 2, method = 'exact', learning_rate=1000, early_exaggeration = 4, n_iter_without_progress = 500, perplexity = 35)
        model = TSNE(n_components = args.dimensions, metric = "precomputed", n_iter=1000, verbose = 2, method = 'barnes_hut', learning_rate=1000, early_exaggeration = 4, n_iter_without_progress = 100, perplexity = 35)

    embedding = model.fit_transform(distances)
    np.savetxt(str(args.output_prefix) + ".embedding.csv", embedding, delimiter=",")

# Draw a k-distances plot
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

# Run clustering
if args.hier:
    sys.stderr.write("Finding " + str(args.clusters) + " clusters with hierarchical clustering\n")
    found_clusters = AgglomerativeClustering(args.clusters, affinity="precomputed", linkage="average").fit_predict(distances)
# Default is DBSCAN
else:
    sys.stderr.write("Clustering samples with DBSCAN\n")
    found_clusters = DBSCAN(eps = args.epsilon, min_samples = args.min_pts).fit_predict(embedding)
    sys.stderr.write("Found " + str(len(set(found_clusters))) + " clusters\n")

# Write output (file for phandango)
csv_out = open(args.output_prefix + str(".csv"), 'w')
csv_sep = ','

csv_out.write(csv_sep.join(("name","cluster:o")) + "\n")
for sample in range(0, found_clusters.size):
    csv_out.write(csv_sep.join((file_num[sample], str(found_clusters[sample]))) + "\n")

csv_out.close()

# Draw plot (see http://scikit-learn.org/stable/auto_examples/cluster/plot_dbscan.html)
unique_labels = set(found_clusters)
colours = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
markers = itertools.cycle(matplotlib.markers.MarkerStyle.filled_markers)
for k, col in zip(unique_labels, colours):
    class_member_mask = (found_clusters == k)
    xy = embedding[class_member_mask]

    if k == -1:
        # Black used for noise.
        col = 'k'
        plt.plot(xy[:, 0], xy[:, 1], marker='o', markerfacecolor=col,
             markeredgecolor='k', markersize=6, linestyle='None')
    else:
        plt.plot(xy[:, 0], xy[:, 1], marker=next(markers), markerfacecolor=col,
             markeredgecolor='k', markersize=10, linestyle='None')

plt.title('Estimated number of clusters: %d' % len(set(found_clusters)))
plt.savefig(args.output_prefix + '.pdf')
plt.close()

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

