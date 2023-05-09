from sys import argv
import pandas as pd
import numpy as np
import pysam

#define sorting function
def sort_counts(x, s1, s2):
    c1, c2 = x[[s1, s2]]
    if np.max([c1 ,c2]) >= 2:
        if c1 >= 2*c2:
            return s1
        if c2 >= 2*c1:
            return s2
        else:
            return -1
    else:
        return -1

#import tables of strain identities and cross parents
nano_strains = pd.read_csv('/mnt/HDD3/lrma/script/nano_strains.csv', index_col=0)
cross_parents = pd.read_csv('/mnt/HDD3/lrma/private_variants/cross_parents.txt', sep=';', header=None, index_col=0, squeeze=True)

#parse strain from command line
s = argv[1]
cross = nano_strains.loc[s, 'cross']

#all crosses except parental strains
if cross != 'P':
    s1, s2 = cross_parents.loc[cross].split(',')

    #import table of private snps between parents
    tab = pd.read_csv(f'/mnt/HDD3/lrma/private_variants/private_variants_{cross}.tab', sep='\t', header=None,
     dtype={0:str,1:np.int32,2:str,3:str,4:str,5:str,6:str})
    for gt_col in (5,6):
        tab[gt_col] = np.where(tab[gt_col].isin(['1|1','1/1']), tab[4], tab[3])
    tab.index = tab[1].values-1

    #import bam
    bam = pysam.AlignmentFile(f'/mnt/HDD3/lrma/minimap_sort/{s}.sort.minimap.bam', 'rb')

    #parse reads and score snps
    reads = {}
    for tig, tab_sub in tab.groupby(0):
        # get list of positions on the contig
        loci = bam.pileup(tig)
        #iterate pileups per position
        for pileup in loci:
            i = pileup.reference_pos
            # if position is a snp that allows to discriminate parents, continue
            if i in tab_sub.index:
                gt1, gt2 = tab_sub.loc[i, [5,6]]
                for (rid, nt) in zip(pileup.get_query_names(), pileup.get_query_sequences()):
                    # init reads dict
                    if rid not in reads:
                        reads[rid] = {s1:0, s2:0}
                    nt = nt.upper()
                    if nt == gt1:
                        reads[rid][s1] += 1
                    elif nt == gt2:
                        reads[rid][s2] += 1

    # sort reads
    sort = pd.DataFrame(reads).T
    sort['sort'] = sort.apply(lambda x : sort_counts(x, s1, s2), axis=1)

    #export read names
    for gt, df in sort.groupby('sort'):
        with open(f'/mnt/HDD3/lrma/sort/{s}.{gt}.reads', 'w') as handle:
            handle.write('\n'.join(df.index))
