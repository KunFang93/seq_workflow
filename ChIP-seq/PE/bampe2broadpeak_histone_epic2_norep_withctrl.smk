import os

# only accept 2 replicates currently
t_bamdir = '/data/kun/modes_SFs_SEs/raw_data/tmp/GSE63109/bam'
peakdir = '/data/kun/modes_SFs_SEs/raw_data/tmp/GSE63109/peaks'
samples = 'T47D_H3K36me3'
blacklist = '/data/kun/Align_Index/Blacklist/lists/hg19-blacklist.v2.bed.gz'
ctrl_f = 'T47D_INPUT.filt.srt.nodup.bam'

if len(os.listdir(t_bamdir)) == 0:
    print("Empty directory or directory not exists")
    exit(1)
if not os.path.exists(peakdir):
    os.mkdir(peakdir)

samples_prefix = glob_wildcards(os.path.join(t_bamdir,"{S}_rep1.filt.srt.nodup.bam")).S
threads = 20

rule all:
    input:
         expand("{out_dir}/{sample}.b150.g2.fs200.blfilt.broadPeak", out_dir=peakdir,sample=samples.split(' '))

rule bampetobedpe_ctrl:
    input:
        t_ctrl = t_bamdir + "/{}".format(ctrl_f)
    params: thread = threads
    output:
           out_ctrl = t_bamdir + "/{}".format(ctrl_f)
    run:
        shell("samtools sort --threads {params.thread} -n {input.t_ctrl} | bedtools bamtobed -bedpe -mate1 -i stdin | gzip -nc > {output.out_ctrl}")
    

rule bampetobedpe:
    input:
         t_rep1 = t_bamdir + "/{sample}.filt.srt.nodup.bam"
    params: thread = threads
    output:
           out_rep1 = t_bamdir + "/{sample}.filt.srt.nodup.bedpe.gz"
    run:
        shell("samtools sort --threads {params.thread} -n {input.t_rep1} | bedtools bamtobed -bedpe -mate1 -i stdin | gzip -nc > {output.out_rep1}")
        
rule epic2callpeak:
    input:
         t_rep1 = t_bamdir +"/{sample}.filt.srt.nodup.bedpe.gz",
         t_ctrl = t_bamdir +"/{}".format(ctrl_f)
    output:
          p_out1 = peakdir + '/{sample}.b150.g2.fs200.broadPeak',

    run:
        shell("epic2 -t {input.t_rep1} -c {input.t_ctrl} -bin 150 -g 2 -fs 200 -o {output.p_out1} ")


rule filtblacklist:
    input:
         raw_input = peakdir + '/{sample}.b150.g2.fs200.broadPeak'
    output:
          final_out = peakdir + '/{sample}.b150.g2.fs200.blfilt.broadPeak'
    run:
        shell("""bedtools intersect -v -a {input.raw_input} -b {blacklist} | sort | uniq | grep -P 'chr[\dY]+[ \t]' > {output.final_out}""")
