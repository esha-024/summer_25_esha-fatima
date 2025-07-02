import sys
from Bio import pairwise2
from Bio.Seq import Seq
from Bio import SeqIO

#-------------Read fasta files--------------

def read_sequence_from_file(file_path):
    try:
        record = next(SeqIO.parse(file_path, "fasta"))
        return str(record.seq)
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        sys.exit(1)

#------------Align sequences-----------------

def align_seq(seq1, seq2, match=2, mismatch=-2, gap_open=-2, gap_extend=-2):
    alignment = pairwise2.align.globalms(seq1, seq2, match, mismatch, gap_open, gap_extend)
    best = alignment[0]
    print("Alignment score:", best.score)
    print("Aligned Seq 1:", best.seqA)
    print("Aligned Seq 2:", best.seqB)
    print("Start:", best.start, "End:", best.end)
    return best

#-----------------Caluculate similarity-------------

def similarity(alignment):
    matches = sum(a == b and a != '-' for a, b in zip(alignment.seqA, alignment.seqB))
    length = len(alignment.seqA)
    simmilarity= (matches / length) * 100 if length > 0 else 0
    print("Similarity:", round(simmilarity , 2), "%")
    return simmilarity

#---------------CalcuLate Gap penalties-----------------

def gap_penalties(alignment):
    gaps = sum(a == '-' or b == '-' for a, b in zip(alignment.seqA, alignment.seqB))
    length = len(alignment.seqA)
    gap_freq = (gaps / length) * 100 if length > 0 else 0
    print("Gaps:", gaps)
    print("Gap Frequency:", round(gap_freq, 2), "%")
    return gap_freq

#------------Identify conserved regions--------------------

def conserved_regions(alignment):
    seq1, seq2 = alignment.seqA, alignment.seqB
    regions, start, match = [], None, []

    for i, (a, b) in enumerate(zip(seq1, seq2)):
        if a == b != '-':
            if start is None:
                start = i
            match.append(a)
        else:
            if start is not None and len(match) >= 20:
                regions.append((start, i - 1, ''.join(match)))
            start, match = None, []

    if start is not None and len(match) >= 20:
        regions.append((start, len(seq1) - 1, ''.join(match)))

    if regions:
        print("Conserved regions (≥20 bp, no gaps):")
        for s, e, seq in regions:
            print(f"Start: {s}, End: {e}, Length: {e - s + 1}, Sequence: {seq}")
    else:
        print("No conserved regions of 20 or more bp.")

        #--------------------MAIN----------------------------

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py seq1.fasta seq2.fasta")
        sys.exit(1)

    file1, file2 = sys.argv[1], sys.argv[2]
    seq1 = read_sequence_from_file(file1)
    seq2 = read_sequence_from_file(file2)

    alignment = align_seq(seq1, seq2)
    similarity(alignment)
    gap_penalties(alignment)
    conserved_regions(alignment)
