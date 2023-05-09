#!/bin/bash

wtdbg2_launch() { # input: 1-base_name
    if [ ! -d $1 ]
    then
        mkdir $1
    fi

    cd $1
    wtdbg2 -p 0 -k 15 -AS 2 -s 0.05 -g 12m -t 1 -i /mnt/HDD3/lrma/sort/${1}.reads.fastq -o $1
    wtpoa-cns -t 1 -i ${1}.ctg.lay.gz -fo ${1}.cns.fa

}

export -f wtdbg2_launch

parallel -j15 wtdbg2_launch {/} ::: $(ls -1 /mnt/HDD3/lrma/sort/*.reads.fastq | sed 's/.reads.fastq//')
