# this script will read an input FASTA file and after analyzing the sequence store the result in a csv file

import csv

def calculate_gc_content(sequence):
    #Calculate GC content as a percentage.
    gc_count = sequence.count('G') + sequence.count('C')
    return int((gc_count / len(sequence)) * 100) 

def validate_sequence(sequence):
    #Check if the sequence contains only valid DNA nucleotides.

    return all(base in ('A','T','C','G') for base in sequence)

def read_fasta(file):
    #Read a FASTA file and return a dictionary of ID: sequence.
    sequences = {}
    with open(file, 'r') as file:
        current_id = None
        for line in file:
            line = line.strip()
            if line.startswith('>'):
                current_id = line[1:]
                sequences[current_id] = ''
            elif current_id:
                sequences[current_id] += line.upper()
    return sequences

def analyze_sequences(sequences):
    #Analyze sequences and return results with unique nucleotides.
    results = []
    unique_nucleotides = set()

    for seq_id, sequence in sequences.items():
        valid = validate_sequence(sequence)
        gc = calculate_gc_content(sequence)
        length = len(sequence)
        results.append({
            'ID': seq_id,
            'Length': length,
            'GC_Content(%)': gc,
            'Valid': valid
        })
        unique_nucleotides.update(sequence)

    return results, unique_nucleotides

def save_to_csv(results, output_file):
    #Save sequence analysis to CSV.
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['ID', 'Length', 'GC_Content(%)', 'Valid']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

#  Main Execution 
fasta_file = 'input.fasta'  # here the fasta file name is "input"
output_csv = 'sequence_analysis.csv'  # rsults will store in a csv file named sequence analysis

sequences = read_fasta(fasta_file)
results, unique_nucleotides = analyze_sequences(sequences)
save_to_csv(results, output_csv)

print("Unique nucleotides found across all sequences:", unique_nucleotides)
print("Analysis complete. Results saved to:", output_csv)
