#! /usr/bin/env python

#Author: Timothy Durham (c) 2018

import argparse
import os
import sys

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('in_bed')
    parser.add_argument('out_bed')
    parser.add_argument('rgb')
    parser.add_argument('--label', default='peakname')
    parser.add_argument('--scale_scores', action='store_true', default=False)
    args = parser.parse_args()

    if os.path.isfile(args.rgb):
        with open(args.rgb) as map_in:
            color_map = dict([line.strip().split() for line in map_in])
    else:
        color_map = None

    with sys.stdin if args.in_bed == 'stdin' else open(args.in_bed) as lines_in:
        with sys.stdout if args.out_bed == 'stdout' else open(args.out_bed, 'w') as out:
            if args.scale_scores is True:
                all_lines = lines_in.readlines()
                all_scores = [float(line.strip().split()[4]) for line in all_lines]
                min_score = min(all_scores)
                max_score = max(all_scores) - min_score
                lines_in = all_lines
            for line in lines_in:
                line = line.strip().split()
                if len(line) == 3:
                    if color_map is None:
                        line.extend([args.label, '1000', '.', line[1], line[1], args.rgb])
                    else:
                        line.extend([args.label, '1000', '.', line[1], line[1], color_map[args.label]])
                elif len(line) == 5:
                    if args.scale_scores is True:
                        newscore = (((float(line[4]) - min_score)/max_score) * 900) + 100
                        line[4] = str(round(newscore))
                    if color_map is None:
                        line.extend(['.', line[1], line[1], args.rgb])
                    else:
                        line.extend(['.', line[1], line[1], color_map[line[3]]])
                out.write('\t'.join(line) + '\n')
