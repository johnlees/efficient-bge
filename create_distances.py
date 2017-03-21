#!/usr/bin/env python
# -*- coding: ASCII -*-
#

import os,sys
import argparse
import subprocess
import itertools
import csv

import numpy as np

separator = "\t"
mash_chunk_size = 500

# Get options
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
io = parser.add_argument_group('Input/output')
options = parser.add_argument_group('Method options')

io.add_argument("-v","--vcf", dest="vcf", help="vcf file, readable by bcftools")
io.add_argument("--assembly_file", help="Tab separated file with sample name and assembly location on each line")
io.add_argument("-o","--output", dest="output_prefix", help="Output prefix", default="clusters")
io.add_argument("--seaview_mat",dest="seaview_mat", help="Pre-computed distances matrix from seaview", default=None)
options.add_argument("-m", "--mash", dest="mash_exec", help="Location of mash executable",default='mash')
options.add_argument("--kmer_size", dest="kmer_size", help="K-mer size for mash sketches", default="21")
options.add_argument("--sketch_size", dest="sketch_size", help="Size of mash sketch", default="1000")
args = parser.parse_args()

# Need one and only one of these inputs (XOR)
if not bool(os.path.isfile(str(args.assembly_file))) ^ bool(os.path.isfile(str(args.seaview_mat))):
    sys.stderr.write("Need one of --assembly_file or --seaview_mat\n")
    sys.exit(1)

if not os.path.isfile(str(args.vcf)):
    sys.stderr.write("Need to supply vcf, to ensure consistency between distances and clusters\n")
    sys.exit(1)
else:
    samples = ef_cluster.read_header(args.vcf)

# Read in files
file_num = []
file_names = []
if os.path.isfile(str(args.assembly_file)):
    with open(str(args.assembly_file)) as f:
        for line in f:
            line = line.rstrip()
            (name, file_name) = line.split(separator)
            file_num.append(name)
            file_names.append(file_name)

# Check matches vcf
ordered = 1
if not set(file_num) == set(samples):
    sys.stderr.write("Sample names in vcf and assembly file do not match\n")
    sys.exit(1)
elif not file_num == samples:
    sys.stderr.write("Sample names in vcf and assembly file are in different orders. This may be slower\n")
    ordered = 0

# Create distance matrix
if os.path.isfile(str(args.seaview_mat)):
    # Read from seaview
    with open(str(args.seaview_mat), 'r') as f:
        reader = csv.reader(f)
        seaview_mat = list(reader)

    # Create file objects as for assemblies
    file_num = seaview_mat[0][1:len(seaview_mat[0])]
    if not file_num == samples:
        sys.stderr.write("Sample names in vcf and seaview file are not matching, and in the same order\n")
        sys.exit(1)

    # Distances remove the row and column names
    distances = np.zeros((len(file_num), len(file_num)))
    line_nr = 0
    for line in seaview_mat:
        if line_nr > 0:
            distances[line_nr-1,:] = np.asarray(line[1:len(seaview_mat[0])])
        line_nr += 1

else:
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
                if ordered == 0:
                    pos1 = samples.index(file_num[file_names.index(name1)])
                    pos2 = samples.index(file_num[file_names.index(name2)])
                else:
                    pos1 = file_num[file_names.index(name1)]
                    pos2 = file_num[file_names.index(name2)]

                if distances[pos1, pos2] == 0:
                    distances[pos1, pos2] = dist
            else:
                break

    except OSError as e:
        sys.stderr.write("Mash Execution failed: " + str(e))

np.save(str(args.output_prefix), distances)

