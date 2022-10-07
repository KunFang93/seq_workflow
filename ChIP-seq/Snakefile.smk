import os
import subprocess
from snakemake.io import glob_wildcards, expand

# useful functions for later fragment esitmation
def cal_readlen(input):
    readlen = int(subprocess.check_output("samtools view " + str(input) +
                                                    " | head -n 1000 | gawk '{print length($10)}' | sort | uniq -c | "
                                                    "perl -ane '$_ =~ s/^[ ]+//g;print $_' | sort -k 1nr,1nr | head -1 "
                                                    "| cut -f2 -d \" \" ",shell=True).split()[0])
    return readlen

def extract_fraglen(inputfile):
    final_fraglen = 180
    with open(inputfile, 'r') as input_file:
        for line in input_file:
            line_info = line.strip().split()
            fraglen_list = line_info[2].split(',')
            for fraglen in fraglen_list:
                if int(fraglen) >= 100 and int(fraglen) <= 250:
                    final_fraglen = int(fraglen)
                    break
                else:
                    continue
    input_file.close()
    return final_fraglen

def estimateFragLen_core(current_input,cc_scores,cc_plot,threads):
    read_len = cal_readlen(current_input)
    # calculate exclusion range maximum
    # TF ChIP-seq:  max(read_len + 10, 50)
    # Histone ChIP-seq:  max(read_len + 10, 100)
    ermax= max(read_len + 10, 100)
    # estimate fragment length
    if not os.path.exists(cc_plot):
        subprocess.check_output("run_spp.R -c={} -p={} -filtchr=chrM -savp={} -out={} -x=-500:{}".format(current_input,threads,cc_plot,cc_scores,ermax),shell=True)
    else:
        print("Plot file already exist, skipped.")

def estimateFragLen(wcs):
    sample_name = wcs.sample
    if config['biorep_method'] == 'pool':
        current_input = trt_dir + '/' + sample_name + "_pool" + samples_suffix
        cc_scores = out_dir + '/' + sample_name + '_pool.cc.qc'
        cc_plot = out_dir + '/' + sample_name +'_pool.cc.plot.pdf'
        estimateFragLen_core(current_input,cc_scores,cc_plot,threads)
        return extract_fraglen(cc_scores)
    elif config['biorep_method'] == 'intersect':
        reps = wcs.rep
        current_input = trt_dir + '/' + sample_name + "_" + reps + samples_suffix
        cc_scores = out_dir + '/' + sample_name + '_'+ reps + '.cc.qc'
        cc_plot = out_dir + '/' + sample_name +'_' + reps + '.cc.plot.pdf'
        estimateFragLen_core(current_input,cc_scores,cc_plot,threads)
        return extract_fraglen(cc_scores)

def all_in(wcs):
    if config['biorep_method'] == 'intersect':
        return expand("{out_dir}/{sample}.PooledInRep1AndRep2.filt.narrowPeak", out_dir = out_dir, sample = samples_prefix)
    elif config['biorep_method'] == 'pool':
        print('wcs',wcs)
        return expand("{out_dir}/{sample}_pool.filt.{peaktype}Peak", out_dir = out_dir, sample = samples_prefix, peaktype = config['peak_type'])

def peakcall_in(wcs):
    sample_name = wcs.sample
    print('peakcall_in',sample_name)
    if peak_caller == 'macs2':
        if config['biorep_method'] == 'pool':
            return expand("{trt_dir}/{sample}_pool{suffix}", trt_dir = trt_dir, sample = sample_name, suffix = samples_suffix)
    elif peak_caller == 'epic2':
        if not paired:
            return expand("{trt_dir}/{sample}_pool.{suffix}bed.gz", trt_dir = trt_dir, sample = sample_name, suffix = samples_suffix_sub)
        else:
            return expand("{trt_dir}/{sample}_pool.{suffix}bedpe.gz", trt_dir = trt_dir, sample = sample_name, suffix = samples_suffix_sub)
    else:
        print("Only support Macs2 and Epic2 in current version")
        exit(1)

def PoolRep_in(wcs):
    sample_name = wcs.sample
    print('PooRep_in',sample_name)
    if len(reps) > 1:
        return expand("{out_dir}/{sample}_{rep}{samples_suffix}",out_dir=trt_dir,sample=sample_name,rep=reps,samples_suffix=samples_suffix)
    else:
        return expand("{out_dir}/{sample}{samples_suffix}",out_dir=trt_dir,sample=sample_name,samples_suffix=samples_suffix)

def PoolRep_ctrl_in(wcs):
    sample_name = wcs.sample
    if len(reps) > 1:
        return expand("{out_dir}/{sample}_{rep}{samples_suffix}",out_dir=ctrl_dir,sample=sample_name,rep=reps,samples_suffix=samples_suffix)
    else:
        return expand("{out_dir}/{sample}{samples_suffix}",out_dir=ctrl_dir,sample=sample_name,samples_suffix=samples_suffix)

def ctrl_in(wcs):
    if peak_caller == 'epic2':
        if paired:
            return "-c " + ctrl_dir + "/Pool_Input.bedpe.gz" if len(ctrl_samples) > 0 else ""
        else:
            return "-c " + ctrl_dir + "/Pool_Input.bed.gz" if len(ctrl_samples) > 0 else ""
    elif peak_caller == 'macs2':
        return "-c " + ctrl_dir + "/Pool_Input.bam" if len(ctrl_samples) > 0 else ""


def ctrl_file(ctrl_dir, ctrl_samples, threads, paired):
    # print("Detect ctrl file..Processing")
    if not os.path.exists("{}/Pool_Input.bam".format(ctrl_dir)):
        outbam = ctrl_dir + "/Pool_Input.bam"
        inputbam = " ".join(["{}/{}".format(ctrl_dir,ctrl_s) for ctrl_s in ctrl_samples])
        if len(ctrl_samples) > 1:
            subprocess.check_output("samtools merge --threads {} -O BAM {} {}".format(threads,outbam,inputbam),shell=True)
        else:
            subprocess.check_output("ln -s {} {}".format(inputbam,outbam),shell=True)
    else:
        print("{}/Pool_Input.bam existed, skipped".format(ctrl_dir))
    if peak_caller == 'epic2':
        if paired:
            if not os.path.exists("{}/Pool_Input.bedpe.gz".format(ctrl_dir)):
                subprocess.check_output("samtools sort --threads {} -n {}/Pool_Input.bam "
                                        "| bedtools bamtobed -bedpe -mate1 -i stdin "
                                        "| pigz -p {} -nc > {}/Pool_Input.bedpe.gz".format(threads,ctrl_dir,threads,ctrl_dir), shell=True)
            else:
                print("{}/Pool_Input.bedpe.gz existed, skipped".format(ctrl_dir))
        else:
            if not os.path.exists("{}/Pool_Input.bed.gz".format(ctrl_dir)):
                subprocess.check_output("samtools sort --threads {} -n {}/Pool_Input.bam "
                                        "| bedtools bamtobed -i stdin "
                                        "| pigz -p {} -nc > {}/Pool_Input.bed.gz".format(threads,ctrl_dir,threads,ctrl_dir), shell=True)
                print("{}/Pool_Input.bed.gz existed, skipped".format(ctrl_dir))

  

trt_dir = config['samples_dir']
ctrl_dir = config['ctrl_dir']
out_dir = config['result_dir']
log_dir = config['log_dir']

# double check if trt_dir exists
if len(os.listdir(trt_dir)) == 0:
    # print("Empty directory or directory not exists")
    exit(1)

# create result and log directory if they are not exists
if not os.path.exists(out_dir):
    os.mkdir(out_dir)
if not os.path.exists(log_dir):
    os.mkdir(log_dir)

# get samples_prefix from configurefile
samples_prefix = config['samples_prefix'].split(' ')
samples_suffix = config['samples_suffix']
samples_suffix_sub = '.'.join(samples_suffix.split('.')[:-1])
blacklist = config['remove_blacklist']
paired = config['paired_end']


if len(samples_suffix_sub) > 1:
    samples_suffix_sub = samples_suffix_sub[1:] + '.'

ctrl_samples = config['ctrl_samples'].split(' ')
if ctrl_samples == ['']:
    ctrl_samples = []

# get some params
nrep = config['biorep']
if nrep > 1:
    reps = ["rep{}".format(i) for i in range(1,nrep+1)]
else:
    reps = []
threads = config['threads']
peak_caller = list(config['peak_caller'].keys())[0]
peak_params = config['peak_caller'].get(peak_caller,'')

if len(ctrl_samples) > 0:
    ctrl_file(ctrl_dir,ctrl_samples,threads,paired)

rule all:
    input:
        all_in

# print("Pool Reps")
rule PoolRep:
    output:
        trt_dir + "/{sample}_pool" + samples_suffix
    input:
        PoolRep_in
    params:
        thread = threads
    run:
        if len(reps) > 1:
            shell("parallel  samtools index ::: {wildcards.sample}*bam")
            shell("multiBamSummary bins -p {params.thread} --bamfiles {input} -o {wildcards.sample}.results.npz")
            shell("plotCorrelation -in {wildcards.sample}.results.npz --corMethod pearson --skipZeros --whatToPlot scatterplot -o {wildcards.sample}_PearsonCorr.png --outFileCorMatrix {wildcards.sample}_PearsonCorr.tab")                 
            shell("samtools merge --threads {params.thread} -O BAM {output} {input}")             
        else:
            shell("ln -s {input} {output}")

if config['biorep_method'] == 'pool':
    if peak_caller == 'macs2':
        # print("Macs2")
        rule macs2_peakcalling:
            output:
                out_dir + "/{sample}_pool_peaks.{peaktype}Peak"
            input:
                trt_bam = peakcall_in,
            params:
                ctrl = ctrl_in,
                param = config['peak_caller'].get('macs2',''),
                cc_score = estimateFragLen
            run:
                shell("macs2 callpeak -t {input} {params.ctrl} -f BAM -n {wildcards.sample}_pool --outdir {out_dir} {params.param}"
                      " --nomodel --shift 0 --extsize {params.cc_score} -B")
    elif peak_caller == 'epic2':
        # print("Epic2")
        rule bam2bed:
            output:
                trt_dir + "/{sample}_pool.{suffix}.gz"
            input:
                trt_dir + "/{sample}_pool" + samples_suffix
            params:
                thread = threads,
                paired_p = "-bedpe -mate1" if paired else ""
            run:
                shell("samtools sort --threads {params.thread} -n {input} | bedtools bamtobed {params.paired_p} -i stdin | pigz -p {params.thread} -nc > {output}")           
         
        rule epic2_peakcalling:
            output:
                out_dir + "/{sample}_pool." + config['peak_type'] + "Peak"
            input:
                trt_bam = peakcall_in,
            params:
                ctrl = ctrl_in,
                param = config['peak_caller'].get('epic2',''),
            run:
                shell("epic2 -t {input} {params.ctrl} {params.param} -o {output}")

    rule filtBlacklist:
        input:
            out_dir + "/{sample}_pool_peaks.{peaktype}Peak"
        output:
            out_dir + "/{sample}_pool.filt.{peaktype}Peak"
        run:
            shell("""bedtools intersect -v -a {input} -b {blacklist} | sort | uniq | grep -P 'chr[\dM]+[ \t]' > {output}""")

elif config['biorep_method'] == 'intersect':
    rule macs2_peakcalling:
        output:
            out_dir + "/{sample}_{rep}_peaks.narrowPeak"
        input:
            trt_dir + "/{sample}_{rep}" + samples_suffix
        params:
            ctrl = ctrl_in,
            param = config['peak_caller'].get('macs2',''),
            cc_score = estimateFragLen
        run:
            shell("macs2 callpeak -t {input} {params.ctrl} -f BAM -n {wildcards.sample}_{wildcards.rep} --outdir {out_dir} {params.param}"
                  " --nomodel --shift 0 --extsize {params.cc_score} -B")

    rule overlappeak:
        input:
             p_rep1 = out_dir + "/{sample}_rep1_peaks.narrowPeak",
             p_rep2 = out_dir + "/{sample}_rep2_peaks.narrowPeak",
             p_pool = out_dir + "/{sample}_pool_peaks.narrowPeak"
        output:
            overlap_out = out_dir + '/{sample}.PooledInRep1AndRep2.narrowPeak'
        run:
            shell("""intersectBed -wo -a {input.p_pool} -b {input.p_rep1} | awk 'BEGIN{{FS="\t";OFS="\t"}}{{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {{print $0}}}}' | cut -f 1-10 | sort | uniq | intersectBed -wo -a stdin -b {input.p_rep2} | awk 'BEGIN{{FS="\t";OFS="\t"}}{{s1=$3-$2; s2=$13-$12; if (($21/s1 >= 0.5) || ($21/s2 >= 0.5)) {{print $0}}}}' | cut -f 1-10 | sort | uniq > {output.overlap_out}""")


    rule filtBlacklist:
        input:
            out_dir + '/{sample}.PooledInRep1AndRep2.narrowPeak'
        output:
            out_dir + '/{sample}.PooledInRep1AndRep2.filt.narrowPeak'
        run:
            shell("""bedtools intersect -v -a {input} -b {blacklist} | sort | uniq | grep -P 'chr[\dM]+[ \t]' > {output}""")







