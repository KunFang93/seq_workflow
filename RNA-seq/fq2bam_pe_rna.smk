import sys
import os

fqdir = '/data/kun/QB_Wang/Fuwen/raw_data/merged_fq'
outdir = '/data/kun/QB_Wang/Fuwen/raw_data/bam'
log_dir = '/data/kun/QB_Wang/Fuwen/raw_data/log'
hisat2_idx = '/data/kun/Align_Index/hisat2_ucsc_hg19/hg19/genome'

if not os.path.exists(outdir):
    os.mkdir(outdir)
if not os.path.exists(log_dir):
    os.mkdir(log_dir)

samples = glob_wildcards(os.path.join(fqdir,"{S}_R1.fastq.gz")).S
threads = 50
mapq_cutoff = 30

rule all:
    input:
        expand("{out_dir}/{sample}.filt.srt.bam", out_dir=outdir,sample=samples)

print("Running fastq file quality control")
rule trimgalore_filt:
    input:
         R1 = lambda sample: os.path.join(fqdir,"{sample}_R1.fastq.gz"),
         R2 = lambda sample: os.path.join(fqdir,"{sample}_R2.fastq.gz")
    params: thread = threads
    output:
          R1_trim = fqdir + "/{sample}_R1_val_1.fq.gz",
          R2_trim = fqdir + "/{sample}_R2_val_2.fq.gz"
    shell:
         "trim_galore --trim-n --paired -j {params.thread} -o {fqdir} {input.R1} {input.R2}"

print("Running mapping process")
rule hisat2_map:
    input:
         R1 = lambda sample: os.path.join(fqdir, "{sample}_R1_val_1.fq.gz"),
         R2 = lambda sample: os.path.join(fqdir, "{sample}_R2_val_2.fq.gz")
    output:
          bam_tmp1 = protected(outdir + "/{sample}.bam")
    params:
           thread = threads
    log:
        log_dir+ "/{sample}_hisat2_map.log"
    shell:
         "hisat2 -p {params.thread} --rna-strandness RF -x {hisat2_idx} -1 {input.R1} -2 {input.R2}| "
         "samtools sort --threads {params.thread} -O BAM -o {output.bam_tmp1} - 2>{log}"


print("Running post mapping fitering")
rule post_map_filt:
    input:
         raw_bam = lambda sample: os.path.join(outdir, "{sample}.bam")
    params: thread = threads, mapq_cut = mapq_cutoff
    output:
          filt_bam = outdir + "/{sample}.filt.srt.bam",
          tmp_file1 = temp(outdir + "/{sample}_tmp_filt.bam"),
          tmp_file2 = temp(outdir + "/{sample}_tmp_filt.fixmate.tmp")
    run:
        # Remove  unmapped, mate unmapped
        # not primary alignment, reads failing platform
        # Remove low MAPQ reads
        # Only keep properly paired reads
        # Obtain name sorted BAM file
        shell("samtools view --threads {params.thread} -F 1804 -f 2 -q {params.mapq_cut} -u {input.raw_bam} "
              "| samtools sort -n --threads {params.thread} -o {output.tmp_file1} - ")
        # Remove orphan reads (pair was removed)
        # and read pairs mapping to different chromosomes
        # Obtain position sorted BAM
        shell("samtools fixmate --threads {params.thread} -r {output.tmp_file1} {output.tmp_file2}")
        shell("samtools view --threads {params.thread} -F 1804 -f 2 -u {output.tmp_file2} | "
              "samtools sort --threads {params.thread} -o {output.filt_bam} -")
