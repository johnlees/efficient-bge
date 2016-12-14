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
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt
import matplotlib.markers

separator = "\t"
mash_chunk_size = 500

# Functions
def run_mash_sketch(file_names, kmer_size, sketch_size):
    # Needs to be done in chunks to prevent massive command lines
    mash_sep = " "
    for chunk in range(((len(file_names) - 1) // mash_chunk_size) + 1):
        chunk_start = chunk * mash_chunk_size
        chunk_end = (chunk+1) * mash_chunk_size
        if chunk_end > len(file_names):
            chunk_end = len(file_names)

        mash_command = str(args.mash_exec) + " sketch -k " kmer_size " -s " sketch_size " -o reference" + str(chunk) + " " + mash_sep.join(file_names[chunk_start:chunk_end])
        retcode = subprocess.call(mash_command, shell=True)
        if retcode < 0:
            sys.stderr.write("Mash sketch failed with signal " + str(retcode) + "\n")
            sys.exit(1)

    if (len(file_names) // mash_chunk_size) > 0:
        paste_join = ".msh reference"
        mash_paste_command = str(args.mash_exec) + " paste reference reference" + paste_join.join([str(chunk) for chunk in range(((len(file_names) - 1) // mash_chunk_size) + 1)]) + ".msh"
        retcode = subprocess.call(mash_paste_command, shell=True)
        if retcode < 0:
            sys.stderr.write("Mash paste failed with signal " + str(retcode) + "\n")
            sys.exit(1)
        else:
            for chunk in range(((len(file_names) - 1) // mash_chunk_size) + 1):
                os.remove("reference" + str(chunk) + ".msh")
    else:
        os.rename("reference0.msh", "reference.msh")

# Get options
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("--assembly_file", help="Tab separated file with sample name and assembly location on each line")
parser.add_argument("-o","--output", dest="output_prefix", help="Output prefix", default="clusters")
parser.add_argument("-d","--dimensions",dest="dimensions",help="Number of t-SNE dimensions to embed to",type=int,default=2)
parser.add_argument("--dist_mat",dest="distmat", help="Pre-computed distances.csv.npy matrix", default=None)
parser.add_argument("--mds", dest="mds", action='store_true', default=False, help="Use MDS instead of t-SNE")
parser.add_argument("--seaview_mat",dest="seaview_mat", help="Pre-computed distances matrix from seaview", default=None)
parser.add_argument("--embedding",dest="embedding_file", help="Pre-computed t-SNE embedding", default=None)
parser.add_argument("-b", "--baps", dest="baps_file", help="BAPS clusters, for comparison", default=None)
parser.add_argument("-m", "--mash", dest="mash_exec", help="Location of mash executable",default='mash')
parser.add_argument("--kmer_size", dest="kmer_size", help="K-mer size for mash sketches", default="21")
parser.add_argument("--sketch_size", dest="sketch_size", help="Size of mash sketch", default="1000")
parser.add_argument("--min_pts", dest="min_pts", help="Minimum number of samples in each cluster",default=5, type=int)
parser.add_argument("--epsilon", dest="epsilon", help="Distance between DBSCAN clusters (pick with knn_plot)",default=0.1,type=float)
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

elif not os.path.isfile(str(args.embedding_file)):
    # Run mash
    sys.stderr.write("Calculating distances with mash\n")
    distances = np.zeros((len(file_num), len(file_num)))
    try:
        if not os.path.isfile("reference.msh"):
            run_mash_sketch(file_names, args.kmer_size, args.sketch_size)

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

# Run t-SNE
if os.path.isfile(str(args.embedding_file)):
    embedding = np.loadtxt(str(args.embedding_file), delimiter=",")
else:
    # default is t-SNE
    if not args.mds:
        sys.stderr.write("Embedding samples into " + str(args.dimensions) + " dimensions with t-SNE\n")
        model = TSNE(n_components = args.dimensions, metric = "precomputed", n_iter=1000, verbose = 2, method = 'exact', learning_rate=1000, early_exaggeration = 4, n_iter_without_progress = 500, perplexity = 35)
    else:
        sys.stderr.write("Embedding samples into " + str(args.dimensions) + " dimensions with MDS\n")
        model = MDS(n_components = args.dimensions, metric=True, dissimilarity = "precomputed", verbose = 1)

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
markers = itertools.cycle(matplotlib.markers.MarkerStyle.filled_markers)
for k, col in zip(unique_labels, colours):
    class_member_mask = (dbscan_clusters == k)
    xy = embedding[class_member_mask]

    if k == -1:
        # Black used for noise.
        col = 'k'
        plt.plot(xy[:, 0], xy[:, 1], marker='o', markerfacecolor=col,
             markeredgecolor='k', markersize=6, linestyle='None')
    else:
        plt.plot(xy[:, 0], xy[:, 1], marker=next(markers), markerfacecolor=col,
             markeredgecolor='k', markersize=10, linestyle='None')

plt.title('Estimated number of clusters: %d' % len(set(dbscan_clusters)))
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
    confusion_out.write(separator.join(["BAPS cluster","Total in BAPS cluster"] + [str(x) for x in sorted(set(dbscan_clusters))]) + "\n")

    score = 0
    for cluster in baps_clusters:
        num_in_cluster = [0] * len(set(dbscan_clusters))
        for lane in baps_sets[cluster-1]:
            fast_cluster = dbscan_clusters[file_index[lane]]
            num_in_cluster[fast_cluster+1] += 1

        score += max(num_in_cluster)
        confusion_out.write(separator.join([str(cluster), str(len(baps_sets[cluster-1]))] + [str(x) for x in num_in_cluster]) + "\n")

    confusion_out.close()
    sys.stderr.write("Compared to " + str(len(baps_clusters)) + " BAPS clusters, " + str(score) + " samples of " + str(len(file_num)) + " are in the same clusters\n")

