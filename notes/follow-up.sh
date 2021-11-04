## clone repository

git clone https://github.com/lescai-teaching/datasets_class.git

## sym link so we do not change the repository itself

mkdir -p variant_calling
cd variant_calling
mkdir -p raw_data
cd raw_data
ln -s /config/workspace/datasets_class/germline_calling/reads/*.gz .
cd ..
mkdir -p alignment
cd alignment

## now we can perform the alignment with BWA

bwa mem \
-t 2 \
-R "@RG\tID:sim\tSM:normal\tPL:illumina\tLB:sim" \
/config/workspace/datasets_class/reference/sequence/Homo_sapiens_assembly38_chr21.fasta \
/config/workspace/variant_calling/raw_data/normal_1.000+disease_0.000_1.fq.gz \
/config/workspace/variant_calling/raw_data/normal_1.000+disease_0.000_2.fq.gz \
| samtools view -@ 8 -bhS -o normal.bam -

## Real time: 152.349 sec; CPU: 257.270 sec

bwa mem \
-t 2 \
-R "@RG\tID:sim\tSM:disease\tPL:illumina\tLB:sim" \
/config/workspace/datasets_class/reference/sequence/Homo_sapiens_assembly38_chr21.fasta \
/config/workspace/variant_calling/raw_data/normal_0.000+disease_1.000_1.fq.gz \
/config/workspace/variant_calling/raw_data/normal_0.000+disease_1.000_2.fq.gz \
| samtools view -@ 8 -bhS -o disease.bam -

## Real time: 157.120 sec; CPU: 273.662 sec


# sort the bam file
samtools sort -o normal_sorted.bam normal.bam
samtools sort -o disease_sorted.bam disease.bam

# index the bam file
samtools index normal_sorted.bam
samtools index disease_sorted.bam


# Marking duplicates

gatk MarkDuplicates \
-I normal_sorted.bam \
-M normal_metrics.txt \
-O normal_md.bam

gatk MarkDuplicates \
-I disease_sorted.bam \
-M disease_metrics.txt \
-O disease_md.bam


### recalibrating

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


#### Apply recalibration

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


### variant calling

cd /config/workspace/variant_calling
mkdir -p variants
cd variants

## first single sample discovery

gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R /config/workspace/datasets_class/reference/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -I /config/workspace/variant_calling/alignment/normal_recal.bam \
   -O normal.g.vcf.gz \
   -ERC GVCF

gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R /config/workspace/datasets_class/reference/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -I /config/workspace/variant_calling/alignment/disease_recal.bam \
   -O disease.g.vcf.gz \
   -ERC GVCF

## then consolidate the 2 files

mkdir -p tmp

### on AMD64 this code ######
## combine the files into one
gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
      -V normal.g.vcf.gz \
      -V disease.g.vcf.gz \
      --genomicsdb-workspace-path compared_db \
      --tmp-dir /config/workspace/variant_calling/variants/tmp \
      -L chr21

### on ARM64 (Mac M1 chip) this code
## combine the files into one
 gatk CombineGVCFs \
   -R /config/workspace/datasets_class/reference/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -V normal.g.vcf.gz \
   -V disease.g.vcf.gz \
   -O cohort.g.vcf.gz

### on AMD64 this code ######
### finally we can call the genotypes jointly
gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R /config/workspace/datasets_class/reference/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -V gendb://compared_db \
   --dbsnp /config/workspace/datasets_class/reference/gatkbundle/dbsnp_146.hg38_chr21.vcf.gz \
   -O results.vcf.gz



### on ARM64 (Mac M1 chip) this code
### finally we can call the genotypes jointly
gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R /config/workspace/datasets_class/reference/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -V cohort.g.vcf.gz \
   --dbsnp /config/workspace/datasets_class/reference/gatkbundle/dbsnp_146.hg38_chr21.vcf.gz \
   -O results.vcf.gz


#### ANNOTATE THE SAMPLE

## download hg38 (UCSC) version of database
snpEff download -v hg38

### to execute snpeff we need to contain the memory
snpEff -Xmx4g ann -v hg38 results.vcf.gz >results_ann.vcf
