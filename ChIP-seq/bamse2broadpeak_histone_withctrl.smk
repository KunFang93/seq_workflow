import os

# only accept 2 replicates currently
t_bamdir = '/data/kun/Cao_R/Bmal1_032822/mapped/chip'
peakdir = '/data/kun/Cao_R/Bmal1_032822/results/chip/pool_method'
if not os.path.exists(peakdir):
    os.mkdir(peakdir)
samples_prefix = glob_wildcards(os.path.join(t_bamdir,"{S}_H3K4_rep1.filt.nodup.bam")).S
macs_g_params = 'mm'
exclusion_range_min = -500
threads = 40
def cal_readlen(input):
    readlen = int(subprocess.check_output("samtools view " + str(input) +
                                                    " | head -n 1000 | gawk '{print length($10)}' | sort | uniq -c | "
                                                    "perl -ane '$_ =~ s/^[ ]+//g;print $_' | sort -k 1nr,1nr | head -1 "
                                                    "| cut -f2 -d \" \" ",shell=True).split()[0])
    return readlen

def extract_fraglen(inputfile):
    final_fraglen = 250
    with open(inputfile, 'r') as input_file:
        for line in input_file:
            line_info = line.strip().split()
            fraglen_list = line_info[2].split(',')
            for fraglen in fraglen_list:
                if int(fraglen) >= 100 and int(fraglen) <= 500:
                    final_fraglen = int(fraglen)
                    break
                else:
                    continue
    input_file.close()
    return final_fraglen

rule all:
    input:
         expand("{out_dir}/{sample}_H3K4.PooledInRep1AndRep2.filt.broadPeak", out_dir=peakdir,sample=samples_prefix)


rule PooltRep:
    input:
         t_rep1 = lambda sample: os.path.join(t_bamdir,"{sample}_H3K4_rep1.filt.nodup.bam"),
         t_rep2 = lambda sample: os.path.join(t_bamdir,"{sample}_H3K4_rep2.filt.nodup.bam"),
    params: thread = threads
    output:
          pool_t = t_bamdir + "/{sample}_H3K4_pool.filt.nodup.bam"
    run:
        shell("samtools merge --threads {params.thread} -O BAM {output.pool_t} {input.t_rep1} {input.t_rep2}")

rule estimateFragLen:
    input:
         t_rep1 = lambda sample: os.path.join(t_bamdir,"{sample}_H3K4_rep1.filt.nodup.bam"),
         t_rep2 = lambda sample: os.path.join(t_bamdir,"{sample}_H3K4_rep2.filt.nodup.bam"),
         t_pool = lambda sample: os.path.join(t_bamdir,"{sample}_H3K4_pool.filt.nodup.bam")
    params: thread = threads, ermin = exclusion_range_min
    output:
          cc_scores_rep1 = peakdir + '/{sample}_H3K4_rep1.cc.qc',
          cc_plot_rep1 = peakdir + '/{sample}_H3K4_rep1.cc.plot.pdf',
          cc_scores_rep2 = peakdir + '/{sample}_H3K4_rep2.cc.qc',
          cc_plot_rep2 = peakdir + '/{sample}_H3K4_rep2.cc.plot.pdf',
          cc_scores_pool = peakdir + '/{sample}_H3K4_pool.cc.qc',
          cc_plot_pool = peakdir + '/{sample}_H3K4_pool.cc.plot.pdf'
    run:
        rep1_read_len = cal_readlen(input.t_rep1)
        rep2_read_len = cal_readlen(input.t_rep2)
        pool_read_len = cal_readlen(input.t_pool)
        print(rep1_read_len, type(rep1_read_len))
        # calculate exclusion range maximum
        # TF ChIP-seq:  max(read_len + 10, 50)
        # Histone ChIP-seq:  max(read_len + 10, 100)
        ermax_rep1 = max(rep1_read_len + 10, 100)
        ermax_rep2 = max(rep2_read_len + 10, 100)
        ermax_pool = max(pool_read_len + 10, 100)

        # estimate fragment length
        shell("run_spp.R -c={input.t_rep1} -p={params.thread} -filtchr=chrM -savp={output.cc_plot_rep1} "
              "-out={output.cc_scores_rep1} -x={params.ermin}:"+str(ermax_rep1))
        shell("run_spp.R -c={input.t_rep2} -p={params.thread} -filtchr=chrM -savp={output.cc_plot_rep2} "
              "-out={output.cc_scores_rep2} -x={params.ermin}:"+str(ermax_rep2))
        shell("run_spp.R -c={input.t_pool} -p={params.thread} -filtchr=chrM -savp={output.cc_plot_pool} "
              "-out={output.cc_scores_pool} -x={params.ermin}:"+str(ermax_pool))

rule macs2callpeak:
    input:
         t_rep1 = lambda sample: os.path.join(t_bamdir,"{sample}_H3K4_rep1.filt.nodup.bam"),
         t_rep2 = lambda sample: os.path.join(t_bamdir,"{sample}_H3K4_rep2.filt.nodup.bam"),
         t_pool = lambda sample: os.path.join(t_bamdir,"{sample}_H3K4_pool.filt.nodup.bam"),
         t_ctrl = lambda sample: os.path.join(t_bamdir,"{sample}_Input.filt.nodup.bam"),
         cc_score_rep1 = peakdir + '/{sample}_H3K4_rep1.cc.qc',
         cc_score_rep2 = peakdir + '/{sample}_H3K4_rep2.cc.qc',
         cc_score_pool = peakdir + '/{sample}_H3K4_pool.cc.qc',
    output:
          p_rep1 = peakdir + '/{sample}_H3K4_rep1_peaks.broadPeak',
          p_rep2 = peakdir + '/{sample}_H3K4_rep2_peaks.broadPeak',
          p_pool = peakdir + '/{sample}_H3K4_pool_peaks.broadPeak'
    run:
        current_sample = wildcards.sample
        shell("macs2 callpeak -t {input.t_rep1} -c {input.t_ctrl} -f BAM -n {wildcards.sample}_H3K4_rep1 --outdir {peakdir} -g {macs_g_params} -p 1e-2"
              " --nomodel --shift 0 --extsize " + str(extract_fraglen(input.cc_score_rep1)) + " -B --broad")
        shell("macs2 callpeak -t {input.t_rep2} -c {input.t_ctrl} -f BAM -n {wildcards.sample}_H3K4_rep2 --outdir {peakdir} -g {macs_g_params} -p 1e-2"
              " --nomodel --shift 0 --extsize " + str(extract_fraglen(input.cc_score_rep2)) + " -B --broad")
        shell("macs2 callpeak -t {input.t_pool} -c {input.t_ctrl} -f BAM -n {wildcards.sample}_H3K4_pool --outdir {peakdir} -g {macs_g_params} -p 1e-2"
              " --nomodel --shift 0 --extsize " + str(extract_fraglen(input.cc_score_pool)) + " -B --broad")

rule overlappeak:
    input:
         p_rep1 = lambda sample: os.path.join(peakdir,"{sample}_H3K4_rep1_peaks.broadPeak"),
         p_rep2 = lambda sample: os.path.join(peakdir,"{sample}_H3K4_rep2_peaks.broadPeak"),
         p_pool = lambda sample: os.path.join(peakdir,"{sample}_H3K4_pool_peaks.broadPeak")
    output:
          overlap_out = peakdir + '/{sample}_H3K4.PooledInRep1AndRep2.broadPeak'
    run:
        shell("""intersectBed -wo -a {input.p_pool} -b {input.p_rep1} | awk 'BEGIN{{FS="\t";OFS="\t"}}{{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {{print $0}}}}' | cut -f 1-10 | sort | uniq | intersectBed -wo -a stdin -b {input.p_rep2} | awk 'BEGIN{{FS="\t";OFS="\t"}}{{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {{print $0}}}}' | cut -f 1-10 | sort | uniq > {output.overlap_out}""")


rule filtblacklist:
    input:
         raw_input = peakdir + '/{sample}_H3K4.PooledInRep1AndRep2.broadPeak',
         pool_input = peakdir + '/{sample}_H3K4_pool_peaks.broadPeak'
    output:
          final_out = peakdir + '/{sample}_H3K4.PooledInRep1AndRep2.filt.broadPeak',
          pool_out = peakdir + '/{sample}_H3K4_pool_peaks.filt.broadPeak'
    run:
        shell("""cat {input.raw_input} | sort | uniq | awk 'BEGIN{{OFS="\t"}}{{if ($5>1000) $5=1000; print $0}}' \
              | grep -P 'chr[\dXY]+[ \t]' > {output.final_out}""")
        shell("""cat {input.pool_input} | sort | uniq | awk 'BEGIN{{OFS="\t"}}{{if ($5>1000) $5=1000; print $0}}' \
              | grep -P 'chr[\dXY]+[ \t]' > {output.pool_out}""")
