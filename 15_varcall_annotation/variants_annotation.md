# Variant Calling Workflow - Annotating the variants

## Download the database

```{bash}
snpEff download -v hg38
```


## Annotate the VCF file

```{bash}
snpEff -Xmx4g ann -v hg38 results.vcf.gz >results_ann.vcf
```