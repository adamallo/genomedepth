# genomedepth
Small C program to efficiently generate summary statistics of sequencing depth and gap size. It uses the result of *bedtools genomecov -bga -ibam* as input.

##Memory usage
This program has been optimized for speed, not RAM usage. This grows linearly with the number of genomic positions and with the number of gaps. The main components that use ram should not use more than sizeof(int)*npos*1.5.

##Example
Using a bed file to subset a bam file efficiently (removing duplicates and secondary alignments) using samtools, and generating the input coverage bedgraph with bedtools. The bedgraph needs to be intersected again with the bed file, since samtools includes all overlapping read pairs, even if they don't align with the desired regions directly.

```
samtoolsview -F 1280 -u -h -b -L test.bed -M test.bam | bedtools genomecov -bga -ibam - | bedtools intersect -a - -b test.bed | genomedepth -g gap.histo -d depth.histo -r 100 -o test_genomedepth.out 
```
