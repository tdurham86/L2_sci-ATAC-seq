#! /usr/bin/env python

#Author: Timothy Durham (c) 2018

import argparse
import numpy
import os
from plumbum import local
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--cell_name_to_cluster_map')
    parser.add_argument('--cuts_bed_file')
    parser.add_argument('--out_dir')
    args = parser.parse_args()

    #get cell names for each cluster
    cluster_map = numpy.loadtxt(args.cell_name_to_cluster_map, delimiter='\t', dtype=object)
    for clust_num in set(cluster_map[:,1]):
        clust_out_dir = os.path.join(args.out_dir, clust_num)
        if not os.path.isdir(clust_out_dir):
            os.makedirs(clust_out_dir)
        with open(os.path.join(clust_out_dir, '{!s}.indextable.txt'.format(clust_num)), 'w') as out:
            out.write('\n'.join(['\t'.join(cluster_map[idx]) for idx in numpy.where(cluster_map[:,1] == clust_num)[0]]) + '\n')
    cell_names_to_clusters = dict([tuple(cluster_map[idx]) for idx in range(cluster_map.shape[0])])

#    with open(args.cell_name_to_cluster_map) as map_in:
#        cell_names_to_clusters = dict([(elt.strip().split()[0], 
#                                        elt.strip().split()[1]) for elt in map_in])

    #parse bam file to make per-cluster bam files
    total_cuts_bed = args.cuts_bed_file
    out_dir = args.out_dir
    if not os.path.isdir(out_dir):
        os.makedirs(out_dir)

    bed_paths = {elt:os.path.join(out_dir, elt, '{!s}.cuts.bed'.format(elt)) for elt in set(cell_names_to_clusters.values())}
    bed_files = {}
    for elt, path in bed_paths.items():
        if not os.path.isdir(os.path.dirname(path)):
            os.makedirs(os.path.dirname(path))
        bed_files[elt] = open(path, 'w')
#    bed_files = {elt:open(path, 'w') for elt,path in bed_paths.items()}
    not_assignable_path = os.path.join(out_dir, 'unassignable.cuts.bed')
    not_assignable = open(not_assignable_path, 'w')
    with open(total_cuts_bed) as lines_in:
        for idx, line in enumerate(lines_in):
            if idx and not idx % 500000:
                sys.stdout.write('Processed {!s} BED records.\n'.format(idx))
                sys.stdout.flush()
            cell_id = line.split()[3].split(':')[0]
            try:
                bed_files[cell_names_to_clusters[cell_id]].write(line)
            except KeyError:
                not_assignable.write(line)
    sys.stdout.write('Processed {!s} BED records.\n'.format(idx))
    sys.stdout.flush()
    not_assignable.close()
    for bed in bed_files.values(): 
        bed.close()

    sys.stdout.write('Splitting BED file complete.\n')
    sys.stdout.flush()
