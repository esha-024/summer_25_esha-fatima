import sys

def read_fasta(file_path):
    """Reads sequences from a FASTA file and returns a dictionary of {header: sequence}."""
    sequences = {}
    try:
        with open(file_path, 'r') as f:
            header = None
            seq_lines = []
            for line in f:
                line = line.strip()
                if not line:
                    continue  # skip empty lines
                if line.startswith(">"):
                    if header and seq_lines:
                        sequences[header] = ''.join(seq_lines)
                    header = line
                    seq_lines = []
                elif header:
                    if all(base in "ATGCatgc" for base in line):  # basic sequence validation
                        seq_lines.append(line)
                   
        if header and seq_lines:
            sequences[header] = ''.join(seq_lines)
    except FileNotFoundError:
        print(f"Error: File not found - {file_path}")
        return None
    
    return sequences

def filter_sequences(sequences, min_length):   #Filters the dictionary to only include sequences >= min_length.
    return {h: s for h, s in sequences.items() if len(s) >= min_length}

def write_fasta(sequences, output_file):       #Writes the filtered sequences to a new FASTA file.
    try:
        with open(output_file, 'w') as f:
            for header, seq in sequences.items():
                f.write(f"{header}\n")
                # wrap lines to 60 chars
                for i in range(0, len(seq), 60):
                    f.write(seq[i:i+60] + '\n')
    except IOError as e:
        print(f"Error writing to file: {e}")
        return False
    return True

if __name__ == "__main__":
 
   input_file = input("Enter the path to the input FASTA file: ").strip()
   sequences = read_fasta(input_file)
    
   print(f"\nTotal sequences read: {len(sequences)}")
    
    # get min length from user
   try:
        min_length = int(input("Enter the minimum sequence length: "))
        if min_length < 0:
            raise ValueError
   except ValueError:
        print("Invalid input. Please enter a non-negative integer for length.")
        sys.exit()

   filtered = filter_sequences(sequences, min_length)
   print(f"Sequences ≥ {min_length} bp: {len(filtered)}")

   output_file = input("Enter the output file path for filtered sequences: ").strip()
   output_written = write_fasta(filtered, output_file)
