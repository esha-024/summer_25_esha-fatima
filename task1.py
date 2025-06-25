import sys

dna = input("Enter DNA sequence: ")

    #  Validate that it contains only A, T, G, C
valid_nucleotides = {'A', 'T', 'G', 'C'}
if not all(base in valid_nucleotides for base in dna):
        print("Invalid DNA sequence.")
        sys.exit(1)

        
if  all(base in valid_nucleotides for base in dna):
        print("Valid DNA sequence. ")
    #  Calculate length and GC content
length = len(dna)
gc_count = dna.count('G') + dna.count('C')
gc_content = gc_count / length

    #  Count each nucleotide
count_A = dna.count('A')
count_T = dna.count('T')
count_G = dna.count('G')
count_C = dna.count('C')

print("Sequence Length:" , length)
print(f"GC Content: {gc_content:.2f}")
print(f"A count: {count_A}")
print(f"T count: {count_T} ")
print(f"G count: {count_G}")
print(f"C count: {count_C}")

    #  Check if GC content > 0.4
if gc_content > 0.4:
        print("High GC content sequence")
if gc_content< 0.4:
        print("Low GC content sequence")    

    #  Reverse using a loop
reversed_seq = ""
for base in dna:
        reversed_seq = base + reversed_seq
print(f"Reversed Sequence: {reversed_seq}")

