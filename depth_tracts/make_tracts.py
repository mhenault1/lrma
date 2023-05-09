#!/home/mathieu/miniconda3/bin/python

from sys import argv
from pandas import read_csv, cut, interval_range, DataFrame
from numpy import median, round

s = argv[1]
cn = int(argv[2])

depth = read_csv(f'/mnt/HDD3/lrma/depth/{s}.fil.depth.gz', sep='\t', header=None)
med = depth[2].median()/cn # divided by cn parameter to account for global ploidy of subg

wdw = 30000
depth['bin'] = cut(depth[1], bins=interval_range(start=1, end=12e6, freq=wdw, closed='left'))
depth_median = depth.groupby([0,'bin'])[2].apply(lambda x: round(median(x)/med)).rename('median').reset_index()
depth_median['Bin'] = depth_median['bin'].apply(lambda x: x.mid)

tracts = []
for chrom, df in depth_median.loc[~depth_median['median'].isna()].groupby(0):

    latest = -1
    idx = 0

    new_tract = [[], latest]

    while idx <= df.shape[0]-1:
        i = df.index[idx]
        call = depth_median.loc[i, 'median']
        if call != latest:
            #if there is an active new tract, dump it first
            if len(new_tract[0]) > 0:
                tracts.append(new_tract)
            #define new tract with the first index in
            new_tract = [[i], call]
            #update last 
            latest = call
            idx +=1
        else:
            new_tract[0].append(i)
            latest = call
            idx +=1
    #final dump
    tracts.append(new_tract)

Tracts = []
for (t, call) in tracts:
    start = depth_median.loc[t[0], 'bin'].left
    end = depth_median.loc[t[-1], 'bin'].right
    chrom = depth_median.loc[t[0], 0]
    Tracts.append([s, chrom, start, end, call])

Tracts = DataFrame(Tracts, columns=['s_subg','chrom','start','end','call'])
Tracts.to_csv(f'/mnt/HDD3/lrma/depth_tracts/{s}.tracts.csv')
