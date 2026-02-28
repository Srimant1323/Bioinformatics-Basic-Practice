"""
SECTION 1: DNA / RNA / PROTEIN BASICS
Author: Srimant Bhardwaj
Program: M.Tech Bioinformatics
"""

import random
from collections import Counter


# 1. Count A, T, G, C
def count_bases(seq):
    seq = seq.upper()
    return {
        "A": seq.count("A"),
        "T": seq.count("T"),
        "G": seq.count("G"),
        "C": seq.count("C")
    }


# 2. Calculate GC content
def gc_content(seq):
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    return (gc / len(seq)) * 100


# 3. Reverse DNA sequence
def reverse_sequence(seq):
    return seq[::-1]


# 4. Reverse complement
def reverse_complement(seq):
    complement = {"A":"T", "T":"A", "G":"C", "C":"G"}
    return "".join(complement[base] for base in seq[::-1].upper())


# 5. Validate DNA sequence
def validate_dna(seq):
    return set(seq.upper()).issubset({"A", "T", "G", "C"})


# 6. Transcribe DNA → RNA
def transcribe_dna_to_rna(seq):
    return seq.upper().replace("T", "U")


# 7. Translate RNA → Protein (partial codon table)
codon_table = {
    "AUG":"M", "UUU":"F", "UUC":"F",
    "UAA":"*", "UAG":"*", "UGA":"*"
}

def translate_rna(rna):
    protein = ""
    for i in range(0, len(rna), 3):
        codon = rna[i:i+3]
        if len(codon) == 3:
            protein += codon_table.get(codon, "X")
    return protein


# 8. Find start and stop codons
def find_start_stop(seq):
    starts = []
    stops = []
    for i in range(len(seq)-2):
        codon = seq[i:i+3]
        if codon == "ATG":
            starts.append(i)
        if codon in ["TAA", "TAG", "TGA"]:
            stops.append(i)
    return starts, stops


# 9. Split into codons
def split_codons(seq):
    return [seq[i:i+3] for i in range(0, len(seq), 3)]


# 10. Find longest ORF
def longest_orf(seq):
    longest = ""
    for i in range(len(seq)):
        if seq[i:i+3] == "ATG":
            for j in range(i, len(seq), 3):
                codon = seq[j:j+3]
                if codon in ["TAA", "TAG", "TGA"]:
                    orf = seq[i:j+3]
                    if len(orf) > len(longest):
                        longest = orf
                    break
    return longest


# 11. Molecular weight of protein (simplified)
aa_weights = {
    "A": 89, "M": 149, "F": 165
}

def protein_mw(protein):
    return sum(aa_weights.get(aa, 0) for aa in protein)


# 12. Amino acid frequency
def aa_frequency(protein):
    return Counter(protein)


# 13. Detect palindromic sequence
def is_palindrome(seq):
    return seq == reverse_complement(seq)


# 14. Find motif occurrences
def find_motif(seq, motif):
    return [i for i in range(len(seq)-len(motif)+1)
            if seq[i:i+len(motif)] == motif]


# 15. Overlapping motifs
def overlapping_motifs(seq, motif):
    positions = []
    for i in range(len(seq)):
        if seq.startswith(motif, i):
            positions.append(i)
    return positions


# 16. Compare two sequences (mismatch count)
def mismatch_count(seq1, seq2):
    return sum(a != b for a, b in zip(seq1, seq2))


# 17. Percentage identity
def percent_identity(seq1, seq2):
    matches = sum(a == b for a, b in zip(seq1, seq2))
    return (matches / len(seq1)) * 100


# 18. Detect low-complexity region (simple version)
def low_complexity(seq):
    counts = Counter(seq)
    most_common = counts.most_common(1)[0][1]
    return most_common / len(seq) > 0.7


# 19. Random DNA generator
def random_dna(length):
    return "".join(random.choice("ATGC") for _ in range(length))


# 20. Hamming distance
def hamming_distance(seq1, seq2):
    return sum(a != b for a, b in zip(seq1, seq2))


# Example execution
if __name__ == "__main__":
    dna = "ATGCGTAA"
    print("Base Count:", count_bases(dna))
    print("GC Content:", gc_content(dna))
    print("Reverse Complement:", reverse_complement(dna))
    print("Random DNA:", random_dna(10))
