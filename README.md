# lrma

Scripts and data for the analysis of transposable elements in hybrid yeast mutation accumulation (MA) lines with long reads

## Main scripts 


### sort_reads.ipynb
Sorting of long-read libraries into subgenomes
### sv.ipynb
Classification of TE orthogroups into SV types
### te_annot.ipynb
Processing of TE annotations

## Subfolders description

### depth/
Compute depth of coverage from sorted libraries mapped on rDNA-masked parental assemblies
### depth_bins/
Binned depth of coverage data for manual visualization of polyploidy/aneuploidy/LOH loci
### depth_tracts/
Depth of coverage data merged by tracts
### liftover/
Conversion of subgenome assembly coordinates into parental coordinates
### medaka/
Polished subgenome-level assemblies
### minimap_aln/
Reciprocal alignment of subgenome assemblies and rDNA-masked parental assemblies for liftover and reverse-liftover
### minimap_misassembly/ (deprecated)
Mapping of sorted libraries on subgenome assemblies
### minimap_sort/
Mapping of long-read libraries for sorting
### mummer_crosses/
Alignment of parental genomes for Ty loci orthology
### mummer_polished/
Alignment of polished subgenome assemblies against parental assemblies
### private_variants/
Extraction of genomic variants that discriminate the parents of each cross
### REannotate_output/
Processing of RepeatMasker TE annotations
### RepeatMasker/
TE annotation in subgenome assemblies
### reverse_liftover/
Conversion of parental coordinates into subgenome assembly coordinates
### script/
Main scripts
### seqkit_stats/
Stats on raw basecalled read libraries
### sort/
Sorting of long-read libraries into subgenomes
### wtdbg2/
Assembly of subgenomes
