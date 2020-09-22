# ChIP-seq pipeline. 
*derived from ENCODE3 ChIP-seq pipeline*

## dependencies
bowtie2, samtools, picard, bedtools, macs2, epic2, phantompeakqualtools  

run  
```
snakemake --cores 20 -p -s fq2bam_pe.smk
```
