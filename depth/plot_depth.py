#!/home/mathieu/miniconda3/bin/python

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from sys import argv

# parse args
s = str(argv[1])
wdw = int(argv[2])

depth = pd.read_csv(f'/mnt/HDD3/lrma/depth/{s}.fil.depth.gz', sep='\t', header=None)
median = depth[2].median()

depth['bin'] = pd.cut(depth[1], bins=pd.interval_range(start=1, end=12e6, freq=wdw, closed='left'))
depth['bin_left'] = depth['bin'].apply(lambda x: x.left)

fig, axes = plt.subplots(nrows=16, figsize=[7,22])

chrom_ax_dict = dict(zip(sorted(set(depth[0])), range(16)))

hist_bin = 0.25
hist_by_bin = depth.groupby([0,'bin_left'])[2].apply(lambda x: pd.Series(np.histogram(x/median, np.arange(0,5,hist_bin))[0])).reset_index()


for chrom, df in hist_by_bin.groupby(0):

    dat = df.pivot_table(index='level_2', columns='bin_left', values=2)
    add_bins = [b for b in np.arange(0,df['bin_left'].max()+wdw, wdw)+1 if b not in dat.columns]
    for b in add_bins:
        dat[b] = 0
    dat = dat.loc[:, sorted(dat.columns)]

    ax = axes[chrom_ax_dict[chrom]]
    ax.imshow(dat, cmap='viridis', origin='lower', aspect='auto', interpolation='none')

depth_median = depth.groupby([0,'bin_left'])[2].apply(lambda x: np.round(np.median(x)/median)).rename('median').reset_index()

for chrom, df1 in depth_median.groupby(0):

    ax = axes[chrom_ax_dict[chrom]]
    ax.axhline(1/hist_bin, lw=1, color='white')
    ax.plot(df1['bin_left'].apply(lambda x: x/wdw), df1['median']/hist_bin, c='red')
    ax.set_title(chrom)

    ax.set_yticks(np.arange(5)/hist_bin)
    ax.set_yticklabels(np.arange(5))

plt.tight_layout()
plt.savefig(f'/mnt/HDD3/lrma/depth/fig/depth_{s}_{wdw}.png', dpi=300)
plt.close()
