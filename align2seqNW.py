from Bio import SeqIO
import numpy as np
import pandas as pd

BLOSUM_62 = pd.read_csv("blosum62.csv")
GAP_PENALTY = -8

def read_fasta(fasta_file):
    seqs = {}
    i = 1
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_name = f"seq_{i}"
        seqs[seq_name] = record.seq
        i += 1
    return seqs

def get_score(aa1, aa2):
    return BLOSUM_62.loc[aa1, aa2]

def two_seq_align(seq1, seq2):
    score_matrix = np.zeros((len(seq1)+1, len(seq2)+1))
    traceback_matrix = np.empty((len(seq1)+1, len(seq2)+1), dtype='U10')

    # Filling the first line
    for i in range(1, len(seq1)+1):
        score_matrix[i, 0] += score_matrix[i-1, 0] + GAP_PENALTY
        traceback_matrix [i, 0] = "UP"
    #Filling the first colomun
    for j in range(1, len(seq2)+1):
        score_matrix[0, j] += score_matrix[0, j-1] + GAP_PENALTY
        traceback_matrix [0, j] = "LEFT"
    
    #Filling the score matrix
    for i in range(len(seq1)):
        for j in range(len(seq2)):
            match = score_matrix[i, j] + get_score(seq1[i],seq2[j])
            delete = score_matrix[i, j+1] + GAP_PENALTY
            insert = score_matrix[i+1, j] + GAP_PENALTY

            #print(get_score(seq2[j], seq1[i]))
            score_matrix[i+1, j+1] = max(match, delete, insert)

            # Fiiling the traceback_matrix
            if max(match, delete, insert) == match:
                traceback_matrix[i+1][j+1] = "DIAG" 
            elif max(match, delete, insert) == delete:
                traceback_matrix[i+1][j+1] = "UP"
            else:
                traceback_matrix[i+1][j+1] = "LEFT"

    # Trace alignment from last cell in matrix
    aligned_seq1 = ""
    aligned_seq2 = ""
    
    i = len(seq1)
    j = len(seq2)
    
    while i > 0 or j > 0:
        if traceback_matrix[i, j] == "DIAG":
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif traceback_matrix[i, j] == "UP":
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2  # UP = Gap dans seq2 (en colonne)
            i -= 1
        else:
            aligned_seq1 = "-" + aligned_seq1  # LEFT = Gap dans seq1 (en ligne)
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1
    
    # optimal alignment with score
    return aligned_seq1, aligned_seq2, score_matrix[len(seq1), len(seq2)]

def show_alignement(aligned_seq1, aligned_seq2, score):
    print("Sequence 1 : ", aligned_seq1)
    print("Sequence 2 : ", aligned_seq2)
    print("Alignement score : ", score)
if __name__== "__main__":
    test = read_fasta("./test2.fasta")
    a_seq1, a_seq2, score = two_seq_align(test['seq_1'], test['seq_2'])
    show_alignement(a_seq1, a_seq2, score)