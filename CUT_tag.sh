# CUT
# Jiulin Chan
# WT_ChIP-seq_H3K27ac_48h.R1.fastq.gz
# WT_ChIP-seq_H3K27ac_48h.R2.fastq.gz

# ------------------------
dir_work="/storage/public/jiulin/lxl/ATAC_CUT_RNA/"
dir_raw="${dir_work}raw/"
dir_clean="${dir_work}clean/"
dir_bam="${dir_work}bam/"
dir_peak="${dir_work}peak/"
ref_fa="/storage/public/xiaolong/hulab/hu/ref_genome_110/Danio_rerio.GRCz11.110.fa"
# ------------------------

name="WT_ChIP-seq_H3K27ac_48h"
fastp \
    -i ${dir_raw}${name}.R1.fastq.gz \
    -I ${dir_raw}${name}.R2.fastq.gz \
    -o ${dir_clean}${name}.R1.fastp.gz \
    -O ${dir_clean}${name}.R2.fastp.gz \
    -h ${dir_clean}${name}.fastp.html \
    -j ${dir_clean}${name}.fastp.json \
    --thread 16 \
    --length_required 8


bowtie2 -p 60 -q -I 10 -X 1000 --dovetail --no-unal --very-sensitive-local --no-mixed --no-discordant -x ${ref_fa%.fa} -1 ${dir_clean}${name}.R1.fastp.gz -2 ${dir_clean}${name}.R2.fastp.gz 2>${dir_bam}log.${name}.bowtie2.log | samtools view -F 4 -u - | samtools sort -@ 5 -m 20G -o ${dir_bam}${name}.bam -  
# filt alignment
samtools view -b -f 2 -q 30 -o ${dir_bam}${name}.f2.q30.bam ${dir_bam}${name}.bam
# preseq
# /public/home/aczhzv5pmn/software/Anaconda/anaconda3-3d/envs/preseq/bin/preseq lc_extrap -e 1e+8 -P -B -D -v -o A1.dat A1.f2.q30.bam 2>A1.log 
# /public/home/aczhzv5pmn/software/Anaconda/anaconda3-3d/envs/R/bin/Rscript preseq.R A1.dat A1.log A1
# picard
sambamba markdup -r -t 20 ${dir_bam}${name}.f2.q30.bam ${dir_bam}${name}.f2.q30.sort.dedup.bam >${dir_bam}log.${name}.sambamba.log 2>&1
samtools index -@ 15 ${dir_bam}${name}.f2.q30.sort.dedup.bam

# # insertion size
# java -jar /public/home/aczhzv5pmn/software/picard/picard-2.25.6/picard.jar CollectInsertSizeMetrics\
    # I=A1.f2.q30.bam\
    # O=cut_tag.insert_size.txt\
    # H=cut_tag.insert_size_histogram.png\                                                                    
    # M=0.5
# call peak
GSIZE=$(bowtie2-inspect ${ref_fa%.fa} | perl -ne 'BEGIN{$n=0} next if(/^>/);s/[^ATGC]//gi; $n+=length($_); END{print int($n*0.85);}') # 计算基因组大小
/home/jiulin/miniforge3/envs/macs2/bin/macs2 callpeak\
    -t ${dir_bam}${name}.f2.q30.sort.dedup.bam\
    -f BAMPE\
    -n ${name}\
    -g 1368780147\
    -B\
    --SPMR\
    --keep-dup all\
    --outdir ${dir_bam}\
    >${dir_bam}log.${name}.macs2.log 2>&1
    