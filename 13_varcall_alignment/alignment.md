# Variant Calling workflow - Alignments


## Before you start

```{bash}
git clone https://github.com/lescai-teaching/datasets_class.git
```


## Prepare your folders


```{bash}
mkdir -p variant_calling
cd variant_calling
mkdir -p raw_data
cd raw_data
```



```{bash}
ln -s /config/workspace/datasets_class/germline_calling/reads/*.gz .
```


```{bash}
cd ..
mkdir -p alignment
cd alignment
```


## Align the reads to the reference


```{bash}
bwa mem \
-t 2 \
-R "@RG\tID:sim\tSM:normal\tPL:illumina\tLB:sim" \
/config/workspace/datasets_class/reference/sequence/Homo_sapiens_assembly38_chr21.fasta \
/config/workspace/variant_calling/raw_data/normal_1.000+disease_0.000_1.fq.gz \
/config/workspace/variant_calling/raw_data/normal_1.000+disease_0.000_2.fq.gz \
| samtools view -@ 8 -bhS -o normal.bam -
```



```{bash}
bwa mem \
-t 2 \
-R "@RG\tID:sim\tSM:disease\tPL:illumina\tLB:sim" \
/config/workspace/datasets_class/reference/sequence/Homo_sapiens_assembly38_chr21.fasta \
/config/workspace/variant_calling/raw_data/normal_0.000+disease_1.000_1.fq.gz \
/config/workspace/variant_calling/raw_data/normal_0.000+disease_1.000_2.fq.gz \
| samtools view -@ 8 -bhS -o disease.bam -
```


## Sort and index BAM files

```{bash}
samtools sort -o normal_sorted.bam normal.bam
samtools sort -o disease_sorted.bam disease.bam
```


```{bash}
samtools index normal_sorted.bam
samtools index disease_sorted.bam
```


## Mark duplicates

```{bash}
gatk MarkDuplicates \
-I normal_sorted.bam \
-M normal_metrics.txt \
-O normal_md.bam

gatk MarkDuplicates \
-I disease_sorted.bam \
-M disease_metrics.txt \
-O disease_md.bam
```


## Base Quality Score Recalibration


### Calculate recalibration


```{bash}
gatk BaseRecalibrator \
   -I normal_md.bam \
   -R /config/workspace/datasets_class/reference/sequence/Homo_sapiens_assembly38_chr21.fasta \
   --known-sites /config/workspace/datasets_class/reference/gatkbundle/dbsnp_144.hg38_chr21.vcf.gz \
   --known-sites /config/workspace/datasets_class/reference/gatkbundle/Mills_and_1000G_gold_standard.indels.hg38_chr21.vcf.gz \
   -O normal_recal_data.table

gatk BaseRecalibrator \
   -I disease_md.bam \
   -R /config/workspace/datasets_class/reference/sequence/Homo_sapiens_assembly38_chr21.fasta \
   --known-sites /config/workspace/datasets_class/reference/gatkbundle/dbsnp_144.hg38_chr21.vcf.gz \
   --known-sites /config/workspace/datasets_class/reference/gatkbundle/Mills_and_1000G_gold_standard.indels.hg38_chr21.vcf.gz \
   -O disease_recal_data.table
```

### Apply recalibration to alignments

```{bash}
gatk ApplyBQSR \
   -R /config/workspace/datasets_class/reference/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -I normal_md.bam \
   --bqsr-recal-file normal_recal_data.table \
   -O normal_recal.bam

gatk ApplyBQSR \
   -R /config/workspace/datasets_class/reference/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -I disease_md.bam \
   --bqsr-recal-file disease_recal_data.table \
   -O disease_recal.bam
```