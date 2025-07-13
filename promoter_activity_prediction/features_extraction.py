
import pandas as pd
import numpy as np
from math import log2
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.preprocessing import StandardScaler
from scipy.sparse import hstack, csr_matrix

# 1. Load and Clean Data
def load_and_prepare(filepath):
    df = pd.read_csv(filepath)
    df = df.dropna(subset=["6)pmSequence", "15)confidenceLevel"])
    df = df[df["15)confidenceLevel"].isin(["S", "W"])]
    df = df.rename(columns={"6)pmSequence": "sequence", "15)confidenceLevel": "confidence"})
    df["label"] = df["confidence"].map({"S": 1, "W": 0})
    return df

# 2. Basic Sequence Features
def compute_gc_content(seq):
    return (seq.count('G') + seq.count('C')) / len(seq) if seq else 0

def compute_at_content(seq):
    return (seq.count('A') + seq.count('T')) / len(seq) if seq else 0

def dinucleotide_repeat_score(seq):
    return sum(1 for i in range(len(seq) - 3) if seq[i:i+2] == seq[i+2:i+4]) / len(seq) if len(seq) >= 4 else 0

def has_tata_box(seq):
    return int("TATA" in seq[:50])

def shannon_entropy(seq):
    if not seq:
        return 0
    base_freq = [seq.count(nuc) / len(seq) for nuc in "ATGC"]
    return -sum(p * log2(p) for p in base_freq if p > 0)

def get_kmers(seq, k=3):
    return ' '.join([seq[i:i+k] for i in range(len(seq) - k + 1)])

def generate_kmers_column(df, k=3):
    kmers_list = []
    for seq in df["sequence"]:
        kmers_list.append(get_kmers(seq, k))
    return kmers_list
def one_hot_encode(seq, max_len=100):
    mapping = {'A': [1, 0, 0, 0],
               'T': [0, 1, 0, 0],
               'G': [0, 0, 1, 0],
               'C': [0, 0, 0, 1]}
    encoded = []
    for base in seq[:max_len]:  # truncate/pad to fixed length
        encoded.extend(mapping.get(base, [0, 0, 0, 0]))  # unknown bases → zero vector
    # pad if sequence shorter than max_len
    if len(seq) < max_len:
        padding = (max_len - len(seq)) * 4
        encoded.extend([0]*padding)
    return encoded

# 3. Feature Extraction
def extract_all_features(df, k=3):
    df = df.copy()
    df["gc_content"] = df["sequence"].apply(compute_gc_content)
    df["at_content"] = df["sequence"].apply(compute_at_content)
    df["dinuc_repeats"] = df["sequence"].apply(dinucleotide_repeat_score)
    df["has_tata"] = df["sequence"].apply(has_tata_box)
    df["entropy"] = df["sequence"].apply(shannon_entropy)
    df["length"] = df["sequence"].apply(len)
    df["kmers"] = generate_kmers_column(df, k)
    df["one_hot"] = df["sequence"].apply(lambda x: one_hot_encode(x, max_len=100))

    vectorizer = CountVectorizer()
    X_kmer = vectorizer.fit_transform(df["kmers"])

    scalar_features = df[["gc_content", "at_content", "dinuc_repeats", "has_tata", "entropy", "length"]]
    X_scalar = StandardScaler().fit_transform(scalar_features)

    
    X_onehot = np.vstack(df["one_hot"].to_numpy())  # convert to dense matrix
    X_combined = hstack([X_kmer, csr_matrix(X_scalar), csr_matrix(X_onehot)])

    return X_combined, df["label"], df
