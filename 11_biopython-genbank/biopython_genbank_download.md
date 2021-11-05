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

### Inspecting a single record

In previous examples we have extracted the records in *fasta* format.
Now we're going to extract them in *GenBank* format, i.e. containing all the features you have seen already on the website.

For this we choose:

- download from the *nucleotide* database
- record type (rettype) equals *gb* (i.e. GenBank)
- we select a text download (retmode)

```
gbresult = Entrez.efetch(db="nucleotide", rettype="gb", retmode="text", id="MW586689")
gb_record = SeqIO.read(gbresult, "genbank")
```

Now that we have fetched the record, and parsed the data, we can inspect the structure of the record we have downloaded, by printing the features:

```
print(str(gb_record.features))
```
Now, this is a bit difficult to read.

So let's have a look at the first one:

```
print(str(gb_record.features[0]))
```

We understand a little better the structure now. This feature has three objects, i.e. identified by the names printed without indentation:

- type
- location
- qualifiers

Then, the object qualifier is a dictionary: I can understand this because of the *key* - *value* pair structure.

Therefore, if I wanted to grab the host organism I'd print:

```
print(str(gb_record.features[0].qualifiers["host"]))
```

### Fetching multiple records

If we pass a list to the *id* parameters, this will fetch records corresponding to each of the accession numbers we have passed.

In order to simplify the inspection of the data, we're fetching some proteins from the protein database instead of nucleotined.

```
protein_list = Entrez.efetch(db="protein", rettype="fasta", retmode="text", id="QRK24690,QRO03507,QRU93410,QRI43434,QRX39425,QRD95445,QRC42505,QRF69711")
```

When we have multiple records to parse, we don't use the method *read* but we need to use the method *parse*, which loops through each of them:

```
records = SeqIO.parse(protein_list, "fasta")
```

Now our *records* is an iterable object, which we can loop through in a for loop:

```
for protein in records:
  print("%s" % protein.id)
  print("%s" % protein.seq)
```

If we wanted to download the fasta sequence of all these records, we'd simply write the following code.
Follow your teacher for a detailed explanation:

```
protein_list = Entrez.efetch(db="protein", rettype="fasta", retmode="text", id="QRK24690,QRO03507,QRU93410,QRI43434,QRX39425,QRD95445,QRC42505,QRF69711")
records = SeqIO.parse(protein_list, "fasta")
handle = open("input.fasta", mode='a')
for record in records:
  handle.write(">%s\n%s\n" % (record.id, record.seq))
  print(">%s\n%s\n" % (record.id, repr(record.seq)))
```

In the example above, we have downloaded protein sequences.

Let's inspect again a GenBank record though, at [this link](https://www.ncbi.nlm.nih.gov/nuccore/MW662150)

Actually, the translation of the spike protein is already present in the whole genome record.

We just need to find a way to extract it.

Let's write a function for it.
Follow your teacher for a detailed explanation of the choices in the code:

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

Now we just have to apply it.

First we get some genome records.

```
gbresultList = Entrez.efetch(db="nucleotide", rettype="gb", retmode="test", id="MW662150,MW662159,MW642248,MW621433,MW580576")
recordsList = SeqIO.parse(gbresultList, "genbank")
```

Then we loop through them, we extract the S protein and we write them in a file.

```
for record in recordsList:
    spike = get_spike_from_gb(record)
    fileName = "genome-" + spike[0] + "_spike-" + spike[1] + ".fasta"
    handle = open(fileName, mode="w")
    handle.write(">%s_%s\n%s" % tuple(spike))
```