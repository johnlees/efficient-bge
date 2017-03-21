#!/usr/bin/env python
# -*- coding: ASCII -*-
#

import os,sys
import argparse
import subprocess
import itertools
import csv
import re

import numpy as np
from sklearn.decomposition import NMF
import matplotlib.pyplot as plt

separator = "\t"

# Get options
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
io = parser.add_argument_group('Input/output')
options = parser.add_argument_group('Method options')

io.add_argument("--vcf", dest="vcf", help="vcf file, readable by bcftools")
io.add_argument("-o","--output", dest="output_prefix", help="Output prefix", default="clusters")
io.add_argument("-b", "--baps", dest="baps_file", help="BAPS clusters, for comparison", default=None)
options.add_argument("-m", "--mac", dest="min_mac", help="Minimum allele count",default=2, type=int)
options.add_argument("--max_clusters", dest="max_clusters", help="Maximum number of clusters", default=10, type=int)
options.add_argument("--min_clusters", dest="min_clusters", help="Minimum number of clusters", default=2, type=int)
options.add_argument("--regularisation", dest="alpha", help="Regularisation constant", default=0, type=float)
options.add_argument("--reg_mixing", dest="mixing", help="Regularisation L1/L2 mixing", default=0, type=float)
args = parser.parse_args()

samples = []
with open(str(args.sample_names)) as f:
    for line in f:
        line = line.rstrip()
        samples.append(line)

# Read in alignment
sys.stderr.write("Reading alignment\n")
tmp_csv = NamedTemporaryFile()
bcftools_command = "bcftools norm -m - " + args.vcf + " | bcftools view -c " + str(args.min_mac) + ":minor - | bcftools query -f '[,%GT]\\n' - | sed 's/^,//' > " + tmp_csv.name()
retcode = subprocess.call(bcftools_command, shell=True)
if retcode < 0:
    sys.stderr.write("bcftools conversion failed with code " + str(retcode) + "\n")
    tmp_csv.close()
    sys.exit(1)

alignment = np.loadtxt(tmp_csv.name(), delimiter=',')
alignment = np.transpose(alignment)
tmp_csv.close()

# Read in header
p = subprocess.Popen(['bcftools view -h ' + args.vcf], stdout=subprocess.PIPE, shell=True)
for line in iter(p.stdout.readline, ''):
    line = line.decode('utf-8')
    line = line.rstrip()
    header_line = re.match("^#CHROM", line)
    if header_line != None:
        samples = line.split(separator)
        del samples[0:8]
        break

# Run NMF
for num_clusters in range(args.min_clusters, args.max_clusters):
    sys.stderr.write("Running NMF on " + str(len(samples)) + " samples with " + str(num_clusters) + " clusters\n")
    model = NMF(n_components = num_clusters, init = 'nndsvd', alpha = args.alpha, l1_ratio = args.mixing, verbose = 2)
    decomposition = model.fit_transform(alignment)

    # normalise each row by dividing by its sum
    # TODO assign to max val, unless close to others in which case bin
    clusters = np.argmax(decomposition, axis = 1)
    decomposition = decomposition/decomposition.sum(axis=1, keepdims = True)

    print(clusters)
    print(decomposition)

    # Write output (file for phandango)
    csv_out = open(args.output_prefix + "." + str(num_clusters) + "clusters.csv", 'w')
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

