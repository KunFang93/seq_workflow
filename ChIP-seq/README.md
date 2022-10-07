# ChIP-seq pipeline. 
*derived from ENCODE3 ChIP-seq pipeline*

## dependencies
check env.yml, which is exported from centos7.0 conda env

## Workflow
<img src="https://github.com/KunFang93/seq_workflow/blob/master/ChIP-seq/ChIP-seq_pipeline_design_100622.pdf" width="900">

## Quick Start  
```
snakemake --cores 20 -p -s Snakefile.smk --configfile config.yaml
```

## TBD
1. Add conservation methods for broad peaks.  
2. Integrate MPSC into pipeline.  
