# Biopython

## Introduction

The goal of Biopython is to make it as easy as possible to use Python for bioinformatics by creating high-quality, reusable modules and classes. Biopython features include parsers for various Bioinformatics file formats (BLAST, Clustalw, FASTA, Genbank), access to online services (NCBI, Expasy)

For more information see the [cookbook](http://biopython.org/DIST/docs/tutorial/Tutorial.html)

The first thing we need to do is importing the libraries into our environment:

```
import Bio
from Bio import Entrez
from Bio import SeqIO
import textwrap
print(Bio.__version__)
Entrez.email = "A.N.Other@example.com"
```

## Simple Downloads

We can now use *Entrez* methods, to download GenBank records based on their identifier.

We use the method *efetch* and we need to specify a few parameters:

- the database we want to search into (db)
- the type of results we would like (rettype)
- the format of the download (retmode)
- the accession number of the record (id)

The code appears as follows:

```
result = Entrez.efetch(
    db="nucleotide", rettype="fasta", retmode="text", id="MW586689"
)
```
Once we have the result, we can *parse* (or *read*) our record: *parsing* means transforming this information into something usable.

Sometimes parsing means converting a format into a readable one. In other cases, like the following, it means *extracting* some information from a more complex structure.

```
seq_record = SeqIO.read(result, "fasta")
```

We can print some of the information we have *parsed* with the following code:

```
print("sequence identifier %s with fasta sequence of length %i" % (seq_record.id, len(seq_record.seq)))
print("A taste of sequence is %s" % repr(seq_record.seq))
```
Follow your teacher, for a more complete explanation of the code above.


## GenBank records



```
gbresult = Entrez.efetch(db="nucleotide", rettype="gb", retmode="test", id="MW586689")
gb_record = SeqIO.read(gbresult, "genbank")
print(str(gb_record.features))
```


```
protein_list = Entrez.efetch(db="protein", rettype="fasta", retmode="text", id="QRK24690,QRO03507,QRU93410,QRI43434,QRX39425,QRD95445,QRC42505,QRF69711")
```


```
records = SeqIO.parse(protein_list, "fasta")
```

```
for protein in records:
  print("%s" % protein.id)
  print("%s" % protein.seq)
```

```
protein_list = Entrez.efetch(db="protein", rettype="fasta", retmode="text", id="QRK24690,QRO03507,QRU93410,QRI43434,QRX39425,QRD95445,QRC42505,QRF69711")
records = SeqIO.parse(protein_list, "fasta")
handle = open("input.fasta", mode='a')
for record in records:
  handle.write(">%s\n%s\n" % (record.id, record.seq))
  print(">%s\n%s\n" % (record.id, repr(record.seq)))
```




```
def get_spike_from_gb(record):
    for feature in record.features:
        if feature.type == "CDS" and 'S' in feature.qualifiers.get("gene"):
            identifier = feature.qualifiers.get("protein_id")[0]
            sequence = feature.qualifiers.get("translation")[0]
            results = [record.id, identifier, sequence]
            #print("for genome %s we got protein %s with sequence\n%s\n" % (record.id, identifier, sequence))
            return results
```


```
gbresultList = Entrez.efetch(db="nucleotide", rettype="gb", retmode="test", id="MW662150,MW662159,MW642248,MW621433,MW580576")
recordsList = SeqIO.parse(gbresultList, "genbank")
```


```
for record in recordsList:
    spike = get_spike_from_gb(record)
    fileName = "genome-" + spike[0] + "_spike-" + spike[1] + ".fasta"
    handle = open(fileName, mode="w")
    handle.write(">%s_%s\n%s" % tuple(spike))
```