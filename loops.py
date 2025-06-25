# ----Iterative Statements-----

# for loop----

# Counting Nucleotides:
dna = "ATGCATGC"
counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
for base in dna:
 if base in counts:
  counts[base] += 1
print("Counts:", counts)

# Processing Genes:
genes = ["BRCA1", "TP53", "EGFR"]
for gene in genes:
 print(f"Analyzing gene: {gene}")

 # while loop----

# Reversing a Sequence with While:
 dna = "ATGC"
 index = len(dna) - 1
 reversed_dna = ""
 while index >= 0:
  reversed_dna += dna[index]
  index -= 1
 print("Reversed:", reversed_dna)