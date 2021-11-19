# Resequencing Workflow - Annotating the variants

In this third and last part of our workflow, we will attach some meaning to each variants we have called in the VCF file generated in the previous class.
## Download the database

First of all we need to download the cache database we need, which is not included in the container to keep its size to a minimum:

```{bash}
snpEff download -v hg38
```


## Annotate the VCF file

Once the above is complete, we can perform the annotation using this local cache:

```{bash}
snpEff -Xmx4g ann -v hg38 results.vcf.gz >results_ann.vcf
```

## Finding the causal variant(s)


### Inspecting the VCF file

We can start searching for the following criteria:

- variants with *HIGH* impact
- variants with a described phenotype, if any
- status of *LoF* 


### Looking into other databases

Depending on our search into the VCF file, we can then search for the position(s) we found, or their *variant identifier* if available (rs number or cosmic ID) in the following databases:

- [gnomad](https://gnomad.broadinstitute.org)
- [clinvar](https://www.ncbi.nlm.nih.gov/clinvar/)

