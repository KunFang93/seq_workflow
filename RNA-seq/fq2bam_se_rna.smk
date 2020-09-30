import os

fqdir = '/data/kun/Jingwei/data/RNA-seq/Jin_V_03192019'
outdir = '/data/kun/Jingwei/data/RNA-seq/bam'
log_dir = '/data/kun/Jingwei/data/RNA-seq/log'
hisat2_idx = '/data/kun/Align_Index/hisat2_ucsc_hg19/hg19/genome'
if not os.path.exists(outdir):
    os.mkdir(outdir)
if not os.path.exists(log_dir):
    os.mkdir(log_dir)
# samples = glob_wildcards(os.path.join(fqdir,"{S}_SE.fastq.gz")).S
samples = glob_wildcards(os.path.join(fqdir,"{S}.fastq.gz")).S
print(samples)
threads = 40
mapq_cutoff = 30

rule all:
    input:
        expand("{out_dir}/{sample}.filt.bam", out_dir=outdir,sample=samples)

print("Running fastq file quality control")
rule trimgalore_filt:
    input:
         R1 = lambda sample: os.path.join(fqdir,"{sample}.fastq.gz"),
    params: thread = threads
    output:
          R1_trim = fqdir + "/{sample}_trimmed.fq.gz",
    shell:
         "trim_galore --trim-n -j {params.thread} -o {fqdir} {input.R1}"

print("Running mapping process")
rule hisat2_map:
    input:
         R1 = lambda sample: os.path.join(fqdir, "{sample}_trimmed.fq.gz"),
    output:
          bam_tmp1 = protected(outdir + "/{sample}.bam")
    params:
           thread = threads
    log:
        log_dir+ "/{sample}_hisat2_map.log"
    shell:
         "hisat2 -p {params.thread} --rna-strandness RF -x {hisat2_idx} -U {input.R1} | "
         "samtools sort --threads {params.thread} -O BAM -o {output.bam_tmp1} - 2>{log}"


print("Running post mapping fitering")
rule post_map_filt:
    input:
         raw_bam = lambda sample: os.path.join(outdir, "{sample}.bam")
    params: thread = threads, mapq_cut = mapq_cutoff
    output:
          filt_bam = outdir + "/{sample}.filt.bam",
    run:
        # Remove  unmapped, mate unmapped
        # not primary alignment, reads failing platform
        # Remove low MAPQ reads
        # Only keep properly paired reads
        # Obtain name sorted BAM file
        shell("samtools view --threads {params.thread} -F 1804 -q {params.mapq_cut} -b {input.raw_bam} -o {output.filt_bam}")
