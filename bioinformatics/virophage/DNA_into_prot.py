#!/usr/bin/env python3

seq = "TCCACCGACGTGGACATGGTCAACGCGCACCCCGTGATCCTCAGCCAGCTCTTTCTGAAGCTCGGACTGGCGTGTCCCGCCCTGGAGAGGTACGTGGACGAGCGGGAGGCGGTGTTGGCTGAG"

# table de complément
complement = {
    "A": "T",
    "T": "A",
    "C": "G",
    "G": "C"
}

genetic_code = {
'TTT':'F','TTC':'F','TTA':'L','TTG':'L',
'TCT':'S','TCC':'S','TCA':'S','TCG':'S',
'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*',
'TGT':'C','TGC':'C','TGA':'*','TGG':'W',
'CTT':'L','CTC':'L','CTA':'L','CTG':'L',
'CCT':'P','CCC':'P','CCA':'P','CCG':'P',
'CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
'CGT':'R','CGC':'R','CGA':'R','CGG':'R',
'ATT':'I','ATC':'I','ATA':'I','ATG':'M',
'ACT':'T','ACC':'T','ACA':'T','ACG':'T',
'AAT':'N','AAC':'N','AAA':'K','AAG':'K',
'AGT':'S','AGC':'S','AGA':'R','AGG':'R',
'GTT':'V','GTC':'V','GTA':'V','GTG':'V',
'GCT':'A','GCC':'A','GCA':'A','GCG':'A',
'GAT':'D','GAC':'D','GAA':'E','GAG':'E',
'GGT':'G','GGC':'G','GGA':'G','GGG':'G'
}

# reverse complement
# seq = "".join(complement[b] for b in seq[::-1])

protein = ""
for i in range(0, len(seq)-2, 3):
    codon = seq[i:i+3]
    protein += genetic_code.get(codon, "X")

print("Séquence inversée :")
print(seq)

print("\nProtéine :")
print(protein)