#-------Input Operation----------------

#Inputting a DNA Sequence:
dna = input("Enter a DNA sequence (e.g., ATGCATGC): ")
print("Received DNA:", dna)


#Inputting a Gene ID:
gene_id = input("Enter a gene ID (e.g., BRCA1): ")
print(f"Processing gene: {gene_id}")

#Inputting Expression Level:
expression = float(input("Enter gene expression level: "))
print(f"Expression level: {expression}")

#------Output Operation----------

#Displaying Sequence Length:
dna = "ATGCATGC"
print("Sequence Length:", len(dna))

#Showing GC Content:
dna = "ATGGCC"
gc_count = dna.count('G') + dna.count('C')
gc_content = gc_count / len(dna)
print(f"GC Content: {gc_content:.2f}")

#Reporting Validation Result:
dna = "ATGX"
is_valid = all(base in 'ATGC' for base in dna)
print("Is Valid DNA:", is_valid)


#-------VARIABLES------

#Storing a DNA Sequence:
dna_sequence = "ATGCATGC"
print("DNA Sequence:", dna_sequence)

# Storing Gene Expression:
gene_expression = 7.2
print("Expression Level:", gene_expression)

# Storing a List of Genes:
genes = ["BRCA1", "TP53", "EGFR"]
print("Genes:", genes)