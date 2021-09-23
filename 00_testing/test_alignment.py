from Bio import pairwise2
from Bio.pairwise2 import format_alignment

seqA="ACTACTAGATTACTTACGGATCAGGTACTTTAGAGGCTTGCAACCA"
seqB="TACTCACGGATGAGGTACTTTAGAGGC"


for a in pairwise2.align.localxx(seqA, seqB):
  print(format_alignment(*a, full_sequences=True))



for a in pairwise2.align.globalxx(seqA, seqB):
  print(format_alignment(*a, full_sequences=True))