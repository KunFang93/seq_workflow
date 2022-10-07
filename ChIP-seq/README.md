# ChIP-seq pipeline. 
* Customized ChIP-seq pipeline, partially derived from ENCODE ChIP-seq pipeline*

## dependencies
check env.yml, which is exported from centos7.0 conda env

## Workflow
<img src="https://github.com/KunFang93/seq_workflow/blob/master/ChIP-seq/workflow.tiff" width="900">

## Quick Start  
```
snakemake --cores 20 -p -s Snakefile.smk --configfile config.yaml
```

## TBD
1. Add pipeline for **Transcription Factors**
2. Add conservation methods for broad peaks.  
3. Integrate MPSC into pipeline.  
