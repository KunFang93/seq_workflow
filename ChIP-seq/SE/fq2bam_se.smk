import sys
import os
# try:
#     fqdir = sys.argv[1]
#     outdir = sys.argv[2]
# except IndexError:
#     print("Usage: python fq2bam.rule fq_directory output_directory")
#     exit(1)

fqdir = '/data/kun/prostate_cancer_project/raw_data/abl/Zhong_chen'
outdir = '/data/kun/prostate_cancer_project/raw_data/abl/bam'
log_dir = '/data/kun/prostate_cancer_project/raw_data/abl/log'
bwt2_idx = '/data/kun/Align_Index/main_chromosome_hg19/index/hg19'
if not os.path.exists(outdir):
    os.mkdir(outdir)
if not os.path.exists(log_dir):
    os.mkdir(log_dir)
# samples = glob_wildcards(os.path.join(fqdir,"{S}_SE.fastq.gz")).S
samples = glob_wildcards(os.path.join(fqdir,"{S}.fastq")).S
print(samples)
threads = 10
mapq_cutoff = 30

rule all:
    input:
        expand("{out_dir}/{sample}_SE.filt.nodup.bam", out_dir=outdir,sample=samples)

print("Running fastq file quality control")
rule trimgalore_filt:
    input:
         R1 = lambda sample: os.path.join(fqdir,"{sample}_SE.fastq.gz"),
    params: thread = threads
    output:
          R1_trim = fqdir + "/{sample}_SE_trimmed.fq.gz",
    shell:
         "trim_galore --trim-n -j {params.thread} -o {fqdir} {input.R1}"

print("Running mapping process")
rule bowtie2_map:
    input:
         R1 = lambda sample: os.path.join(fqdir, "{sample}_SE_trimmed.fq.gz"),
    output:
          bam_tmp1 = protected(outdir + "/{sample}_SE.bam")
    params:
           thread = threads
    log:
        log_dir+ "/{sample}_bowtie2_map.log"
    shell:
         "bowtie2 -p {params.thread} -x {bwt2_idx} -U {input.R1} | "
         "samtools sort --threads {params.thread} -O BAM -o {output.bam_tmp1} - 2>{log}"


print("Running post mapping fitering")
rule post_map_filt:
    input:
         raw_bam = lambda sample: os.path.join(outdir, "{sample}_SE.bam")
    params: thread = threads, mapq_cut = mapq_cutoff
    output:
          filt_bam = outdir + "/{sample}_SE.filt.bam",
    run:
        # Remove  unmapped, mate unmapped
        # not primary alignment, reads failing platform
        # Remove low MAPQ reads
        # Only keep properly paired reads
        # Obtain name sorted BAM file
        shell("samtools view --threads {params.thread} -F 1804 -q {params.mapq_cut} -b {input.raw_bam} -o {output.filt_bam}")

print("Running deduplicate process")
rule remove_duplicates:
    input:
        filt_bam = lambda sample: os.path.join(outdir, "{sample}_SE.filt.bam")
    params: thread = threads
    output:
          final_bam = outdir + "/{sample}_SE.filt.nodup.bam",
          dup_file_qc = outdir + "/{sample}_SE.dup.qc",
          tmp_out = temp(outdir + "/{sample}_tmp_filt_dedup.bam")
    run:
        # Mark duplicates
        shell("picard MarkDuplicates INPUT={input.filt_bam} OUTPUT={output.tmp_out} METRICS_FILE={output.dup_file_qc}"
              " VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=false")
        # Remove duplicates
        # Index final position sorted BAM
        # Create final name sorted BAM
        shell("samtools view --threads {params.thread} -F 1804 -b {output.tmp_out} > {output.final_bam}")
        # Index Final BAM file
        shell("samtools index {output.final_bam}")
        # remove input filt.srt file
        shell("rm {input.filt_bam}")
