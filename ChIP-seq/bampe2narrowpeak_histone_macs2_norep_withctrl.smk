import os

# only accept 2 replicates currently
t_bamdir = '/data/kun/modes_SFs_SEs/raw_data/tmp/GSE63109/bam'
peakdir = '/data/kun/modes_SFs_SEs/raw_data/tmp/GSE63109/peaks'
samples = 'T47D_H3K4me2'
blacklist = '/data/kun/Align_Index/Blacklist/lists/hg19-blacklist.v2.bed.gz'
ctrl_f = 'T47D_INPUT.filt.srt.nodup.bam'
macs_g_params = 'hs'

if len(os.listdir(t_bamdir)) == 0:
    print("Empty directory or directory not exists")
    exit(1)
if not os.path.exists(peakdir):
    os.mkdir(peakdir)

threads = 20

rule all:
    input:
         expand("{out_dir}/{sample}.filt.narrowPeak", out_dir=peakdir,sample=samples.split(' '))

rule macs2callpeak:
    input:
         t_rep1 = t_bamdir + "/{sample}.filt.srt.nodup.bam",
         t_ctrl = t_bamdir + "/" + ctrl_f
    output:
          p_rep1 = peakdir + '/{sample}.narrowPeak',
    run:
        shell("macs2 callpeak -t {input.t_rep1} -c {input.t_ctrl} -f BAMPE -n {wildcards.sample} -g {macs_g_params} --outdir {peakdir} -B")

rule filtblacklist:
    input:
         raw_input = peakdir + '/{sample}.narrowPeak'
    output:
          final_out = peakdir + '/{sample}.filt.narrowPeak'
    run:
        shell("""cat {input.raw_input} | sort | uniq | awk 'BEGIN{{OFS="\t"}}{{if ($5>1000) $5=1000; print $0}}' | grep -P 'chr[\dY]+[ \t]' > {output.final_out}""")