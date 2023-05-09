#!/bin/bash
CROSS=$(cut -d ';' -f1 /mnt/HDD3/lrma/private_variants/cross_parents.txt)
PARENTS=$(cut -d ';' -f2 /mnt/HDD3/lrma/private_variants/cross_parents.txt)

parallel -j13 --link 'bcftools annotate -x INFO/AN /mnt/HDD2/seqdata/freebayes/freebayes.parents.combined.filtered.vcf.gz | \
bcftools view -s {2} | vcfallelicprimitives | bcftools view -g^het -g^miss -c2 -C2 | bcftools query -f "%CHROM\t%POS\t%TYPE\t%REF\t%ALT[\t%GT]\n" > \
 private_variants_{1}.tab' ::: $CROSS ::: $PARENTS


