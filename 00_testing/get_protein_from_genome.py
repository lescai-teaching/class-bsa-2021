import Bio
from Bio import Entrez
from Bio import SeqIO
import textwrap
Entrez.email = "A.N.Other@example.com"


gbresult = Entrez.efetch(db="nucleotide", rettype="gb", retmode="test", id="MW662150")
gb_record = SeqIO.read(gbresult, "genbank")


def get_spike_from_gb(record):
    for feature in record.features:
        if feature.type == "CDS" and 'S' in feature.qualifiers.get("gene"):
            identifier = feature.qualifiers.get("protein_id")[0]
            sequence = feature.qualifiers.get("translation")[0]
            results = [record.id, identifier, sequence]
            #print("for genome %s we got protein %s with sequence\n%s\n" % (record.id, identifier, sequence))
            return results

get_spike_from_gb(gb_record)


gbresultList = Entrez.efetch(db="nucleotide", rettype="gb", retmode="test", id="MW662150,MW662159,MW642248,MW621433,MW580576")
recordsList = SeqIO.parse(gbresultList, "genbank")

for record in recordsList:
    spike = get_spike_from_gb(record)
    fileName = "genome-" + spike[0] + "_spike-" + spike[1] + ".fasta"
    handle = open(fileName, mode="w")
    handle.write(">%s_%s\n%s" % tuple(spike))

