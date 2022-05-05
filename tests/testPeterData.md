```bash
cd /data/home/ywang120/myData/PeterData/bbsrcPipelineOnPeterData
ln -s ../rawReads/12* .

REF=/data/home/ywang120/myData/PeterData/reference_genome/dmel_r5.44/dmel-all-chromosome-r5.44.fasta


for r1 in $(ls *read1.fq.gz);do
    echo ${r1}
    r2=${r1//read1/read2}
    tt=${r1##*-}
    sid=${tt:0:2}
    echo ${sid}
    r1_new=pet_dmel_01_S${sid}_R1.fastq.gz
    r2_new=pet_dmel_01_S${sid}_R2.fastq.gz
    mv ${r1} ${r1_new}
    mv ${r2} ${r2_new}
done
```



```
X       1       22422827        M       1
X       1       22422827        F       2
*       *       *       M       2
*       *       *       F       2

```

```
S32	F
S43	M
S62	M
S74	M
S75	F
S79	F
S84	F
S87	F
S88	F
S89	F
S90	F
S93	F
S94	F
S95	F
```



```python

import os


READ_PATH = "/data/home/ywang120/myData/PeterData/bbsrcPipelineOnPeterData/"
PREFIX = "dmel_01"



REF   = '/data/home/ywang120/myData/PeterData/reference_genome/dmel_r5.44/dmel-all-chromosome-r5.44.fasta'
CHROM = ['2L', '2R', '3L', '3R', 'X', '4']

########################################################

def ind_id(sample_id):
    return sample_id.split('_')[3]

SAMPLE = [x.replace("_R1.fastq.gz","") for x in os.listdir(READ_PATH) if x.endswith('R1.fastq.gz')]
#######################

rule final:
    input:
        PREFIX + "_allsample.vcf.gz.tbi"

rule bwa_map:
    input:
        fa = REF,
        r1 = lambda w:"{t1}{t2}_R1.fastq.gz".format(t1=READ_PATH, t2=w.sample),
        r2 = lambda w:"{t1}{t2}_R2.fastq.gz".format(t1=READ_PATH, t2=w.sample)
    output:    
        ss = temp("{sample}.bam")
    params:
        rg = lambda  w: r"@RG\tID:{tt}\tSM:{tt}".format(tt=ind_id(w.sample))
    threads: 8
    shell:
        "bwa mem -R '{params.rg}' -t {threads} {input.fa} {input.r1} {input.r2} | samtools view -Sb -@ {threads} - > {output.ss}"

rule samtools_sort:
    input:
        ii = "{sample}.bam"
    output:
        # DO NOT delete the file, picard MarkDuplicates remove duplicate not in pairs which may leave a single read end
        # which may not perform simulation that convert bam  to fastq
        oo = "{sample}.sort.bam" 
    threads: 8
    shell:
        "samtools sort -t {threads} {input.ii} > {output.oo}"

rule samtools_index:
    input:
        ii = "{sample}.sort.bam" 
    output:
        oo = "{sample}.sort.bam.bai"
    shell:
        "samtools index {input.ii}"

####
rule markRemoveDuplicate:
    input:
        ii = "{sample}.sort.bam",
        hh = "{sample}.sort.bam.bai"
    output:
        oo = "{sample}.sort.noDup.bam",
        qq = temp("{sample}.sort.noDup.metrics.txt")
    params:
        pp = "REMOVE_DUPLICATES=true ASSUME_SORTED=true VALIDATION_STRINGENCY=LENIENT"
    shell:
        "picard MarkDuplicates INPUT={input.ii} OUTPUT={output.oo} METRICS_FILE={output.qq} "
        "{params.pp}"
            
           
rule indexNodupBam:
    input:
        ii = "{sample}.sort.noDup.bam"
    output:
        oo = "{sample}.sort.noDup.bam.bai"
    shell:
        "samtools index {input.ii}"
        
       
rule makeBamList:
    input: 
        ii = expand("{sample}.sort.noDup.bam.bai", sample=SAMPLE)
    output:
        oo = "bam.list"
    shell:
        "find . -regex '.*sort.noDup.bam$' > {output.oo}"

rule bcfcall:
    input:
        ii = "bam.list",
        fa = REF,
        jj = "ploidy.config",
        kk = "samples.config"
    output:
        oo = temp("temp_{chr}.bcf")
    params:
        pp = "-a DP,AD,ADF,ADR --max-depth 250 --min-MQ 20 --regions {chr}"
    shell:
        "bcftools mpileup -Ou "
        "--fasta-ref {input.fa} "
        "{params.pp} "
        "--bam-list {input.ii} | "
        "bcftools call -mv -Ob -f GQ "
        "--ploidy-file  {input.jj} "
        "--samples-file {input.kk} "
        "-o {output.oo}"
        
rule bcf2vcf:
    input:
        ii = "temp_{chr}.bcf"
    output:
        oo = "temp_{chr}.vcf.gz"
    shell:
        "bcftools view {input.ii} | bgzip > {output.oo}"

rule gatherChrVCF:
    input:
        ii = expand("temp_{chr}.vcf.gz", chr=CHROM)
    output:
        oo = PREFIX + "_allsample.vcf.gz"
    shell:
        "bcftools concat -Ou {input.ii} | bcftools sort -Oz -o {output.oo}"

rule indexVcf:
    input:
        ii = PREFIX + "_allsample.vcf.gz"
    output:
        oo = PREFIX + "_allsample.vcf.gz.tbi"
    shell:
        "tabix -p vcf {input.ii}"
```




```bash
SPECIES='dmel'
FAM='34'
PED="/data/home/ywang120/myData/bbsrc/idata/ped/${SPECIES}_${FAM}.ped"
CHR="chromosome.config"
REF="/data/home/ywang120/myData/PeterData/reference_genome/dmel_r5.44/dmel-all-chromosome-r5.44.fasta"
##########################


IN1=dmel_01_allsample.vcf.gz
OUT1=${SPECIES}_${FAM}_allsample.snp.vcf.gz
gatk SelectVariants \
    -R ${REF} \
    -V ${IN1} \
    -select-type-to-include SNP \
    -O ${OUT1}

```