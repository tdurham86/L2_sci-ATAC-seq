#! /usr/bin/env python

#Author: Timothy Durham (c) 2018

import argparse
import numpy
from plumbum import local
import pyBigWig
import sys

def group_peaks(narrowPeak):
    grouped_peaks = {}
    with sys.stdin if narrowPeak == 'stdin' else open(narrowPeak) as lines_in:
        for line in lines_in:
            line = line.strip().split()
            line[1] = int(line[1])
            line[2] = int(line[2])
            grouped_peaks.setdefault(tuple(line[:3]), []).append(line)
    return grouped_peaks

def split_peak(peak_lines, bw_data, min_peak_radius=25):
    '''Takes narrowPeak lines corresponding to the same peak region
    and a pyBigWig object containing the signal for that region.
    Identifies the peak summit(s) and adjusts the peak coordinates for 
    each peak based on the prominence of its
    signal using a changepoint-detecting convolution. Returns the 
    bed+3 lines with the coordinates set to these new boundaries.
    '''
    bw_slice = [peak_lines[0][0], 
                peak_lines[0][1] - 150, 
                peak_lines[0][2] + 150]
    if bw_slice[1] < 0:
        bw_slice[1] = 0
    if bw_slice[2] > bw_data[0].chroms(bw_slice[0]):
        bw_slice[2] = bw_data[0].chroms(bw_slice[0])
    data_vector = numpy.array([elt.values(*bw_slice) for elt in bw_data]).sum(axis=0)
    data_vector[numpy.where(numpy.isnan(data_vector))[0]] = 0
    data_vector[numpy.where(data_vector < 0)[0]] = 0
    smoothed_data_vector = numpy.convolve(data_vector, numpy.ones(100)/100, mode='valid')
    test_step = numpy.convolve(numpy.hstack([numpy.ones(100),
                                             0 - numpy.ones(100)]),
                               numpy.ones(6)/6, mode='valid')
    conv_dist = numpy.convolve(smoothed_data_vector - numpy.mean(smoothed_data_vector), test_step, mode='valid')
    pos_vals = conv_dist >= 0
    neg_vals = conv_dist <= 0
    summit_idx_list = []
    for summit_idx in numpy.where(numpy.logical_and(pos_vals[:-1], ~pos_vals[1:]))[0]:
        if (numpy.sum(pos_vals[summit_idx - min_peak_radius:summit_idx]) >= min_peak_radius and
            numpy.sum(neg_vals[summit_idx:summit_idx + min_peak_radius + 1]) >= min_peak_radius):
            summit_idx_list.append(summit_idx)

    xvals = numpy.arange(conv_dist.shape[0])+97
    new_summits = []
    for idx, summit_loc in enumerate(summit_idx_list):
        lower_end = summit_loc - min_peak_radius
        lower_bound = numpy.where(conv_dist[:lower_end] < 0)[0]
        lower_bound = (numpy.where(conv_dist[lower_bound[-1]:lower_end + 1] > 0)[0][0] + 97 + lower_bound[-1] if len(lower_bound) > 0 else xvals[0])
        upper_start = summit_loc + min_peak_radius
        upper_bound = numpy.where(conv_dist[upper_start:] > 0)[0]
        if len(upper_bound) > 0:
            upper_bound = numpy.where(conv_dist[upper_start:upper_bound[0] + upper_start] < 0)[0][-1] + upper_start + 97
        else:
            upper_bound = xvals[-1]
        peak_coord = peak_lines[0][:3]
        peak_name = peak_lines[0][3].rstrip('abcdefghijklmnopqrstuvwxyz')
        if len(summit_idx_list) > 1:
            peak_name += chr(97 + idx)
        peak_score = data_vector[summit_loc + 50]
        new_summits.append([peak_coord[0], 
                            peak_coord[1] + lower_bound - 100, 
                            peak_coord[1] + upper_bound - 100,
                            peak_name,
                            peak_score])
        if new_summits[-1][1] < peak_lines[0][1]:
            new_summits[-1][1] = peak_lines[0][1]
        if new_summits[-1][2] > peak_lines[0][2]:
            new_summits[-1][2] = peak_lines[0][2]
    return new_summits

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('narrowPeak')
    parser.add_argument('signal_bigWig', nargs='+')
    parser.add_argument('--outfile')
    parser.add_argument('--ignore_single_summits', action='store_true', default=False)
    args = parser.parse_args()

    grouped_peaks = group_peaks(args.narrowPeak)
    bw_data = [pyBigWig.open(elt) for elt in args.signal_bigWig]

    with open(args.outfile, 'w') if args.outfile is not None else sys.stdout as out:
        for coord, lines in sorted(grouped_peaks.items(), key=lambda x:x[0]):
            if args.ignore_single_summits is True and len(lines) == 1:
                out.write('\t'.join([str(elt) for elt in lines[0][:5]]) + '\n')
            else:
                try:
                    to_write = split_peak(lines, bw_data)
                except:
                    sys.stderr.write(str(coord) + '\n')
                    sys.stderr.write('\n'.join(['\t'.join([str(elt) for elt in line]) for line in lines]) + '\n')
                    raise
                for line in to_write:
                    out.write('\t'.join([str(elt) for elt in line]) + '\n')
