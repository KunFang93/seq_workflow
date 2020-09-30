import os

# only accept 2 replicates currently
t_bamdir = '/data/kun/prostate_cancer_project/raw_data/abl/bam/broad_histone'
peakdir = '/data/kun/prostate_cancer_project/raw_data/abl/broadpeak_histone'
blacklist = '/data/kun/Align_Index/Blacklist/lists/hg19-blacklist.v2.bed.gz'

if len(os.listdir(t_bamdir)) == 0:
    print("Empty directory or directory not exists")
    exit(1)
if not os.path.exists(peakdir):
    os.mkdir(peakdir)

samples_prefix = glob_wildcards(os.path.join(t_bamdir,"{S}_rep1.filt.srt.nodup.bam")).S
threads = 20

rule all:
    input:
         expand("{out_dir}/{sample}.PooledInRep1AndRep2.filt.broadPeak", out_dir=peakdir,sample=samples_prefix)

rule PooltRep:
    input:
         t_rep1 = lambda sample: os.path.join(t_bamdir,"{sample}_rep1.filt.srt.nodup.bam"),
         t_rep2 = lambda sample: os.path.join(t_bamdir,"{sample}_rep2.filt.srt.nodup.bam"),
    params: thread = threads
    output:
          pool_t = t_bamdir + "/{sample}_pool.filt.srt.nodup.bam"
    run:
        shell("samtools merge --threads {params.thread} -O BAM {output.pool_t} {input.t_rep1} {input.t_rep2}")

rule bampetobedpe:
    input:
         t_rep1 = t_bamdir +"/{sample}_rep1.filt.srt.nodup.bam",
         t_rep2 = t_bamdir +"/{sample}_rep2.filt.srt.nodup.bam",
         t_pool = t_bamdir +"/{sample}_pool.filt.srt.nodup.bam"
    params: thread = threads
    output:
           out_rep1 = t_bamdir + "/{sample}_rep1.filt.srt.nodup.bedpe.gz",
           out_rep2 = t_bamdir + "/{sample}_rep2.filt.srt.nodup.bedpe.gz",
           out_pool = t_bamdir + "/{sample}_pool.filt.srt.nodup.bedpe.gz"
    run:
        shell("samtools sort --threads {params.thread} -n {input.t_rep1} | bedtools bamtobed -bedpe -mate1 -i stdin | gzip -nc > {output.out_rep1}")
        shell("samtools sort --threads {params.thread} -n {input.t_rep2} | bedtools bamtobed -bedpe -mate1 -i stdin | gzip -nc > {output.out_rep2}")
        shell("samtools sort --threads {params.thread} -n {input.t_pool} | bedtools bamtobed -bedpe -mate1 -i stdin | gzip -nc > {output.out_pool}")

rule epic2callpeak:
    input:
         t_rep1 = t_bamdir +"/{sample}_rep1.filt.srt.nodup.bedpe.gz",
         t_rep2 = t_bamdir +"/{sample}_rep2.filt.srt.nodup.bedpe.gz",
         t_pool = t_bamdir +"/{sample}_pool.filt.srt.nodup.bedpe.gz"
    output:
          p_rep1 = peakdir + '/{sample}_rep1_peaks.broadPeak',
          p_rep2 = peakdir + '/{sample}_rep2_peaks.broadPeak',
          p_pool = peakdir + '/{sample}_pool_peaks.broadPeak'
    run:
        shell("epic2 -t {input.t_rep1} -bin 150 -g 2 -fs 200 -o {output.p_rep1} ")
        shell("epic2 -t {input.t_rep2} -bin 150 -g 2 -fs 200 -o {output.p_rep2} ")
        shell("epic2 -t {input.t_pool} -bin 150 -g 2 -fs 200 -o {output.p_pool} ")

rule overlappeak:
    input:
         p_rep1 = lambda sample: os.path.join(peakdir,"{sample}_rep1_peaks.broadPeak"),
         p_rep2 = lambda sample: os.path.join(peakdir,"{sample}_rep2_peaks.broadPeak"),
         p_pool = lambda sample: os.path.join(peakdir,"{sample}_pool_peaks.broadPeak")
    output:
          overlap_out = peakdir + '/{sample}.PooledInRep1AndRep2.broadPeak'
    run:
        shell("""intersectBed -wo -a {input.p_pool} -b {input.p_rep1} | awk 'BEGIN{{FS="\t";OFS="\t"}}{{s1=$3-$2; s2=$9-$8; if (($13/s1 >= 0.5) || ($13/s2 >= 0.5)) {{print $0}}}}' | cut -f 1-6 | sort | uniq | intersectBed -wo -a stdin -b {input.p_rep2} | awk 'BEGIN{{FS="\t";OFS="\t"}}{{s1=$3-$2; s2=$9-$8; if (($13/s1 >= 0.5) || ($13/s2 >= 0.5)) {{print $0}}}}' | cut -f 1-6 | sort | uniq > {output.overlap_out}""")


rule filtblacklist:
    input:
         raw_input = peakdir + '/{sample}.PooledInRep1AndRep2.broadPeak'
    output:
          final_out = peakdir + '/{sample}.PooledInRep1AndRep2.filt.broadPeak'
    run:
        shell("""bedtools intersect -v -a {input.raw_input} -b {blacklist} | sort | uniq | grep -P 'chr[\dXY]+[ \t]' > {output.final_out}""")
