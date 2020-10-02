#! /usr/bin/env python

#Author: Timothy Durham (c) 2018

import argparse
import glob
import numba
import numpy
import os
from plumbum import local
from scipy import sparse as sps
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('peak_bed')
    parser.add_argument('cells_indextable')
    parser.add_argument('--out_file')
    parser.add_argument('--barcode_len', type=int, default=37)
    parser.add_argument('-S', '--sum_reads', 
                        action='store_true', default=False)
    args = parser.parse_args()

    #allpeaks = '../5_peaks/all_subclusters_summits.slop50.merged.bed'
    allpeaks = args.peak_bed
    peakcoords = numpy.loadtxt(allpeaks, dtype=str)
#MODIFIED TO ALLOW EXACTLY OVERLAPPING PEAKS
#    peakcoords_dict = {'{!s}:{!s}-{!s}'.format(*elt):idx 
#                       for idx, elt in enumerate(peakcoords)}
#    print(len(peakcoords_dict))
    peakcoords_dict = {}
    num_peaks = None
    for idx, elt in enumerate(peakcoords):
        peakcoords_dict.setdefault('{!s}:{!s}-{!s}'.format(*elt), []).append(idx)
        num_peaks = idx + 1
    print(num_peaks)

    #cellfile = '../../all_L2_sciATAC_merged.readdepth.cells.indextable.noWXbatch.txt'
    cellfile = args.cells_indextable
    with open(cellfile) as cells_in:
        cellnames_dict = {elt.split()[0]:idx for idx, elt in enumerate(cells_in)}
    print(len(cellnames_dict))

    cell_by_peak = sps.dok_matrix((len(cellnames_dict),num_peaks), dtype=int)

    print(args.sum_reads)

    #read in lines from bedtools intersect -wo
#    input_file = open('./test_input.bed')
    for idx, line in enumerate(sys.stdin):
#    for idx, line in enumerate(input_file.readlines()):
        line = line.strip().split()
        try:
            cell_idx = cellnames_dict[line[3][:args.barcode_len]]
        except KeyError:
            continue
        peak_idx = peakcoords_dict['{!s}:{!s}-{!s}'.format(*line[6:9])]
        try:
            if args.sum_reads is True:
                cell_by_peak[cell_idx, peak_idx] += 1
            else:
                cell_by_peak[cell_idx, peak_idx] = 1
        except:
            print(cell_idx, peak_idx, idx, line)
            raise
#    input_file.close()
    print(idx)
    cell_by_peak_binary = cell_by_peak.tocsr().toarray()

    #remove any cells with no peaks or peaks with no cells
    cellnames_array = numpy.array([elt[0] for elt in sorted(cellnames_dict.items(), key=lambda x:x[1])], dtype=object)
    peaknames_array = peakcoords.copy()
    while True:
        if numpy.any(numpy.sum(cell_by_peak_binary, axis=1) < 1):
            idx_to_keep = numpy.where(numpy.sum(cell_by_peak_binary, axis=1) > 0)[0]
            cell_by_peak_binary = cell_by_peak_binary[idx_to_keep,:]
            cellnames_array = cellnames_array[idx_to_keep]
        elif numpy.any(numpy.sum(cell_by_peak_binary, axis=0) < 1):
            idx_to_keep = numpy.where(numpy.sum(cell_by_peak_binary, axis=0) > 0)[0]
            cell_by_peak_binary = cell_by_peak_binary[:,idx_to_keep]
            peaknames_array = peaknames_array[idx_to_keep,:]
        else:
            break

    #with open('../cell_by_peak_binary.bow', 'w') as out:
    with open(args.out_file, 'w') if args.out_file else sys.stdout as out:
        out.write('{!s}\n{!s}\n1000\n'.format(*cell_by_peak_binary.shape))
        for cell_idx, peak_idx in zip(*numpy.where(cell_by_peak_binary)):
#            out.write('{!s}\t{!s}\t{!s}\n'.format(cell_idx + 1, peak_idx + 1, 1))
            out.write('{!s}\t{!s}\t{!s}\n'.format(cell_idx + 1, peak_idx + 1, cell_by_peak_binary[cell_idx, peak_idx]))

    with open(os.path.splitext(args.out_file)[0] + '.zeros_filtered.indextable.txt', 'w') as out:
        for idx in range(cellnames_array.shape[0]):
            out.write('{!s}\thypodermis\n'.format(cellnames_array[idx]))

    with open(os.path.splitext(args.out_file)[0] + '.zeros_filtered.bed', 'w') as out:
        for idx in range(peaknames_array.shape[0]):
            out.write('\t'.join(peaknames_array[idx])+'\n')
