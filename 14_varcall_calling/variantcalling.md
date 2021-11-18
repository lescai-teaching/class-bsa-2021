# Variant Calling Workflow - Calling the variants


## Prepare your folders


```{bash}
cd /config/workspace/variant_calling
mkdir -p variants
cd variants
```


## Identify variant sites

```{bash}
gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R /config/workspace/datasets_class/reference/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -I /config/workspace/variant_calling/alignment/normal_recal.bam \
   -O normal.g.vcf.gz \
   -ERC GVCF
```


```{bash}
gatk --java-options "-Xmx4g" HaplotypeCaller  \
   -R /config/workspace/datasets_class/reference/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -I /config/workspace/variant_calling/alignment/disease_recal.bam \
   -O disease.g.vcf.gz \
   -ERC GVCF
```


## Combine data

```{bash}
mkdir -p tmp
```



```{bash}
gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport \
    -V normal.g.vcf.gz \
    -V disease.g.vcf.gz \
    --genomicsdb-workspace-path compared_db \
    --tmp-dir /config/workspace/variant_calling/variants/tmp \
    -L chr21
```


## Calculate Genotypes


```{bash}
gatk --java-options "-Xmx4g" GenotypeGVCFs \
   -R /config/workspace/datasets_class/reference/sequence/Homo_sapiens_assembly38_chr21.fasta \
   -V gendb://compared_db \
   --dbsnp /config/workspace/datasets_class/reference/gatkbundle/dbsnp_146.hg38_chr21.vcf.gz \
   -O results.vcf.gz
```