#!/home/mathieu/miniconda3/bin/python

from sys import argv
from pandas import read_csv, cut, interval_range, DataFrame
from numpy import median, round

s = argv[1]
cn = int(argv[2])

depth = read_csv(f'/mnt/HDD3/lrma/depth/{s}.fil.depth.gz', sep='\t', header=None)
depth.rename({0:'chrom', 1:'pos', 2:'depth'}, axis=1, inplace=True)
med = depth['depth'].median()/cn # divided by cn parameter to account for global ploidy of subg

wdw = 10000
depth['bin'] = cut(depth['pos'], bins=interval_range(start=1, end=12e6, freq=wdw, closed='left'))
depth_median = depth.groupby(['chrom','bin']).apply(lambda x: median(x['depth'])/med).rename('depth').reset_index()
depth_median['Bin'] = depth_median['bin'].apply(lambda x: x.mid)

depth_median[['chrom', 'Bin', 'depth']].to_csv(f'/mnt/HDD3/lrma/depth_bins/{s}.bins.tsv', sep='\t')
