import os

# only accept 2 replicates currently
t_bamdir = '/data/kun/prostate_cancer_project/raw_data/abl/bam/narrow_histone'
peakdir = '/data/kun/prostate_cancer_project/raw_data/abl/narrowpeak_histone'
blacklist = '/data/kun/Align_Index/Blacklist/lists/hg19-blacklist.v2.bed.gz'
macs_g_params = 'hs'

if len(os.listdir(t_bamdir)) == 0:
    print("Empty directory or directory not exists")
    exit(1)
if not os.path.exists(peakdir):
    os.mkdir(peakdir)
samples_prefix = glob_wildcards(os.path.join(t_bamdir,"{S}_rep1.filt.srt.nodup.bam")).S
threads = 20

rule all:
    input:
         expand("{out_dir}/{sample}.PooledInRep1AndRep2.filt.narrowPeak", out_dir=peakdir,sample=samples_prefix)

rule PooltRep:
    input:
         t_rep1 = lambda sample: os.path.join(t_bamdir,"{sample}_rep1.filt.srt.nodup.bam"),
         t_rep2 = lambda sample: os.path.join(t_bamdir,"{sample}_rep2.filt.srt.nodup.bam"),
    params: thread = threads
    output:
          pool_t = t_bamdir + "/{sample}_pool.filt.srt.nodup.bam"
    run:
        shell("samtools merge --threads {params.thread} -O BAM {output.pool_t} {input.t_rep1} {input.t_rep2}")


rule macs2callpeak:
    input:
         # t_rep1 = lambda sample: os.path.join(t_bamdir,"{sample}_rep1.filt.srt.nodup.bam"),
         # t_rep2 = lambda sample: os.path.join(t_bamdir,"{sample}_rep2.filt.srt.nodup.bam"),
         # t_pool = lambda sample: os.path.join(t_bamdir,"{sample}_pool.filt.srt.nodup.bam"),
         t_rep1 = t_bamdir +"/{sample}_rep1.filt.srt.nodup.bam",
         t_rep2 = t_bamdir +"/{sample}_rep2.filt.srt.nodup.bam",
         t_pool = t_bamdir +"/{sample}_pool.filt.srt.nodup.bam",
    output:
          p_rep1 = peakdir + '/{sample}_rep1_peaks.narrowPeak',
          p_rep2 = peakdir + '/{sample}_rep2_peaks.narrowPeak',
          p_pool = peakdir + '/{sample}_pool_peaks.narrowPeak'
    run:
        shell("macs2 callpeak -t {input.t_rep1} -f BAMPE -n {wildcards.sample}_rep1 -g {macs_g_params} --outdir {peakdir}")
        shell("macs2 callpeak -t {input.t_rep2} -f BAMPE -n {wildcards.sample}_rep2 -g {macs_g_params} --outdir {peakdir}")
        shell("macs2 callpeak -t {input.t_pool} -f BAMPE -n {wildcards.sample}_pool -g {macs_g_params} --outdir {peakdir}")

rule overlappeak:
    input:
         p_rep1 = lambda sample: os.path.join(peakdir,"{sample}_rep1_peaks.narrowPeak"),
         p_rep2 = lambda sample: os.path.join(peakdir,"{sample}_rep2_peaks.narrowPeak"),
         p_pool = lambda sample: os.path.join(peakdir,"{sample}_pool_peaks.narrowPeak")
    output:
          overlap_out = peakdir + '/{sample}.PooledInRep1AndRep2.narrowPeak'
    run:
        shell("""intersectBed -wo -a {input.p_pool} -b {input.p_rep1} | awk 'BEGIN{{FS="\t";OFS="\t"}}{{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {{print $0}}}}' | cut -f 1-10 | sort | uniq | intersectBed -wo -a stdin -b {input.p_rep2} | awk 'BEGIN{{FS="\t";OFS="\t"}}{{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {{print $0}}}}' | cut -f 1-10 | sort | uniq > {output.overlap_out}""")


rule filtblacklist:
    input:
         raw_input = peakdir + '/{sample}.PooledInRep1AndRep2.narrowPeak'
    output:
          final_out = peakdir + '/{sample}.PooledInRep1AndRep2.filt.narrowPeak'
    run:
        shell("""cat {input.raw_input} | sort | uniq | awk 'BEGIN{{OFS="\t"}}{{if ($5>1000) $5=1000; print $0}}' | grep -P 'chr[\dXY]+[ \t]' > {output.final_out}""")
