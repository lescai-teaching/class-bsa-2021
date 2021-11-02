



```
import Bio
from Bio import Entrez
from Bio import SeqIO
import textwrap
print(Bio.__version__)
Entrez.email = "A.N.Other@example.com"
```


```
result = Entrez.efetch(
    db="nucleotide", rettype="fasta", retmode="text", id="MW586689"
)
```


```
seq_record = SeqIO.read(result, "fasta")
```

```
print("sequence identifier %s with fasta sequence of length %i" % (seq_record.id, len(seq_record.seq)))
print("A taste of sequence is %s" % repr(seq_record.seq))
```




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