import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import gc_fraction
from Bio.Restriction import EcoRI, Analysis
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.express as px

fasta_path = "Saccharomyces_cerevisiae.fasta"
output_csv = "sequence_analysis.csv"

#  Sequence Input and Validation
def read_and_validate_fasta(file_path):
    try:
        records = list(SeqIO.parse(file_path, 'fasta'))
        if not records:
            print("ERROR: No sequences found in FASTA file.")
            return []

        valid_records = []
        valid_nucleotides = {'A', 'T', 'G', 'C', 'a', 't', 'g', 'c', 'N', 'n'}

        for record in records:
            seq_sequence = str(record.seq)
            if all(base in valid_nucleotides for base in seq_sequence):
                valid_records.append(record)
            else:
                print(f"WARNING: Sequence {record.id} contains invalid characters and will be skipped.")

        return valid_records

    except FileNotFoundError:
        print(f"ERROR: File '{file_path}' not found.")
        return []
    except Exception as e:
        print(f"ERROR reading FASTA file: {e}")
        return []

#  2. GC Content Calculation 
def calculate_gc(seq):
    try:
        return round(gc_fraction(seq) * 100, 2)
    except Exception:
        return 0.0


# 4. Transcription and Translation 
def transcribe_and_translate(dna_seq):
    try:
        rna = dna_seq.transcribe()
        protein = rna.translate(to_stop=True)
        return str(rna)[:50], str(protein)[:50]
    except Exception as e:
        print(f"ERROR in transcription/translation: {e}")
        return "", ""



# 5. EcoRI Site Counting 
def count_ecori_sites(seq):
    try:
        ana = Analysis(RestrictionBatch=[EcoRI])
        result = ana.full_analysis(Seq(seq))
        return len(result.get(EcoRI, []))
    except Exception:
        return 0

#  6. Analyze All Sequences 
def analyze_sequences(records):
    result_list = []

    for rec in records:
        seq = str(rec.seq).upper()
        gc = calculate_gc(seq)
        #codon = codon_usage(seq)
        rna, protein = transcribe_and_translate(Seq(seq))
        eco = count_ecori_sites(seq)
        length = len(seq)

        result_list.append({
            "Sequence_ID": rec.id,
            "Length": length,
            "GC_Content": gc,
            "EcoRI_Sites": eco,
            "RNA": rna,
            "Protein": protein,
            
        })

    return result_list
def plot_gc_bar(df):
    plt.figure(figsize=(10, 5))
    plt.bar(df["Sequence_ID"], df["GC_Content"], color='green')
    plt.xlabel("Sequence ID")
    plt.ylabel("GC Content (%)")
    plt.title("GC Content per Sequence")
    plt.xticks(rotation=90)
    plt.tight_layout()
    plt.savefig("gc_content_bar.png")
    plt.close()

def plot_histogram(df):
    plt.figure(figsize=(8, 5))
    plt.hist(df["Length"], bins=20, color="skyblue", edgecolor="black")
    plt.title("Distribution of Sequence Lengths")
    plt.xlabel("Sequence Length")
    plt.ylabel("Frequency")
    plt.savefig("length_histogram.png")
    plt.close()

def Scatter_Plot(df):
    fig = px.scatter(
        df, x="GC_Content", y="EcoRI_Sites", text="Sequence_ID",
        size="Length", color="Length",
        title="GC Content vs. EcoRI Sites"
    )
    fig.write_html("restriction_scatter.html")

if __name__ == "__main__":
    records = read_and_validate_fasta(fasta_path)

    if not records:
        print("No valid records to analyze.")
    else:
        print(f"Analyzing {len(records)} valid sequences...")
        analysis_results = analyze_sequences(records)

        # Convert to DataFrame
        df = pd.DataFrame(analysis_results)

        # Save to CSV with error check
        try:
            df.to_csv(output_csv, index=False)
            print(f"Results saved to {output_csv}")
        except Exception as e:
            print(f"ERROR writing to CSV: {e}")

        # Generate and save plots
        plot_gc_bar(df)
        plot_histogram(df)
        Scatter_Plot(df)