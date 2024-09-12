"""
Multiple Sequence Alignment

This script performs multiple sequence alignment of protein sequences using the Needleman-Wunsch algorithm combined with UPGMA and other methods.

Usage:
======
    python multiple_alignement_NW.py file.fasta
"""

__author__ = ("Oumayma Meskini")
__contact__ = ("oumayma.meskini@etu.u-paris.fr")
__date__ = "2024-09-12"

from Bio import SeqIO
import numpy as np
import pandas as pd
import sys
import os
import timeit

BLOSUM_62 = pd.read_csv("blosum62.csv")
GAP_PENALTY = -8


def read_fasta(fasta_file):
    """
    Reads sequences from a FASTA file and returns them as a dictionary.

    Args:
        fasta_file (str): The path to the FASTA file.

    Returns:
        dict: A dictionary where keys are sequence identifiers (IDs) and values are sequences 
              (as Biopython Seq objects). Each key-value pair represents one sequence in the FASTA file
    """
    seqs = {}
    i = 1
    for record in SeqIO.parse(fasta_file, "fasta"):
        seqs[record.id] = record.seq
        i += 1
    return seqs

def get_score(aa1, aa2):
    """
    Retrieves the score for a pair of amino acids from the BLOSUM62 substitution matrix.

    Args:
        aa1 (str): The first amino acid in the pair, represented by its one-letter code.
        aa2 (str): The second amino acid in the pair, represented by its one-letter code.

    Returns:
        int: The score for the pair of amino acids as defined in the BLOSUM62 matrix.

    Raises:
        KeyError: If either amino acid is not found in the BLOSUM62 matrix.
    """
    return BLOSUM_62.loc[aa1, aa2]

def get_aa(cluster, pos) :
    """
    Extract amino acids (or gap characters) at a specific position from a cluster of aligned sequences.
    
    Args:
        cluster (list of str): A list of aligned sequences (strings) where each sequence is expected 
                               to have the same length.
        pos (int): The position from which to extract the amino acids in each sequence.

    Returns:
        list of str: A list containing the amino acids (or gap characters) at the specified position 
                     for each sequence in the cluster.
    """
    aa = []
    for seq in cluster:
        aa.append(seq[pos])
    return aa

def calc_multi_score(aa_cluster1, aa_cluster2):
    """
    Calculate the average alignment score between two clusters of amino acids for a certain position.
    The score is calculated as follows:
    - If either amino acid in a pair is a gap ('-'), a gap penalty (GAP_PENALTY) is applied.
    - Otherwise, the pairwise score is obtained from the substitution matrix BLOSUM62.

    Args:
        aa_cluster1 (list of str): A list of amino acids (or gaps) from the first cluster.
        aa_cluster2 (list of str): A list of amino acids (or gaps) from the second cluster.

    Returns:
        float: The average alignment score between the two clusters.
    """
    score = 0
    counter = 0
    for aa1 in aa_cluster1 : 
        for aa2 in aa_cluster2 :
            if aa1 == "-" or aa2 == "-":
                score += GAP_PENALTY
            else :
                score += get_score(aa1, aa2)
            counter += 1
    return score/counter

def calc_distance(distance1, distance2, counter):
    """
    Calculate distance between two clusters or one cluster and a sequence.

    Args:
        distance1 (float): The distance between the first cluster and an external reference (another cluster or sequence).
        distance2 (float): The distance between the second cluster and the same external reference.
        counter (int): The number of elements in the first cluster (used to weight the distance).

    Returns:
        float: The calculated average distance between the two clusters or between a cluster and a sequence.
    """
    distance = ((counter * distance1 + distance2) / (counter + 1))
    return distance

def get_min(distance_matrix):
    """
    Find the minimum value in a distance pandas matrix and its corresponding row and column.

    Args:
        distance_matrix (pd.DataFrame): A pandas DataFrame representing the distance matrix, 
                                        where the rows and columns correspond to sequence or cluster identifiers.

    Returns:
        tuple: A tuple containing:
            - min (float): The minimum value in the distance matrix.
            - row_min (str): The label of the row containing the minimum value.
            - col_min (str): The label of the column containing the minimum value.
            - min_row_index (int): The index of the row containing the minimum value.
            - min_col_index (int): The index of the column containing the minimum value.
    """
    min = distance_matrix.min().min()
    col_min = distance_matrix.min().idxmin()
    row_min = distance_matrix[col_min].idxmin()
    min_row_index = distance_matrix.index.get_loc(row_min)
    min_col_index = distance_matrix.columns.get_loc(col_min)
    return min, row_min, col_min, min_row_index, min_col_index

def two_seq_align(seq1, seq2):
    """
    Perform global sequence alignment between two sequences using dynamic programming - Needleman-wunsch.

    Args:
        seq1 (str): The first sequence to be aligned.
        seq2 (str): The second sequence to be aligned.

    Returns:
        tuple: A tuple containing:
            - aligned_seq1 (str): The aligned version of the first sequence with gaps (`-`) inserted where needed.
            - aligned_seq2 (str): The aligned version of the second sequence with gaps (`-`) inserted where needed.
            - score (float): The alignement score.
    """
    # Initialize a score matrix and a traceback matrix
    score_matrix = np.zeros((len(seq1)+1, len(seq2)+1))
    traceback_matrix = np.empty((len(seq1)+1, len(seq2)+1), dtype='U10')

    # Filling the first line with gaps
    for i in range(1, len(seq1)+1):
        score_matrix[i, 0] += score_matrix[i-1, 0] + GAP_PENALTY
        traceback_matrix [i, 0] = "UP"
    # Filling the first colomun with gaps
    for j in range(1, len(seq2)+1):
        score_matrix[0, j] += score_matrix[0, j-1] + GAP_PENALTY
        traceback_matrix [0, j] = "LEFT"
    
    # Filling the score matrix
    for i in range(len(seq1)):
        for j in range(len(seq2)):
            match = score_matrix[i, j] + get_score(seq1[i],seq2[j])
            delete = score_matrix[i, j+1] + GAP_PENALTY
            insert = score_matrix[i+1, j] + GAP_PENALTY

            score_matrix[i+1, j+1] = max(match, delete, insert)

            # Fiiling the traceback_matrix
            if max(match, delete, insert) == match:
                traceback_matrix[i+1][j+1] = "DIAG" 
            elif max(match, delete, insert) == delete:
                traceback_matrix[i+1][j+1] = "UP"
            else:
                traceback_matrix[i+1][j+1] = "LEFT"

    # Trace alignment from last cell in matrix
    # Initialise aligned sequences
    aligned_seq1 = ""
    aligned_seq2 = ""
    
    # index of the score alignement in the matrix
    i = len(seq1)
    j = len(seq2)
    
    while i > 0 or j > 0:
        # Matching amino acids
        if traceback_matrix[i, j] == "DIAG":
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        # UP = Gap in the second sequence (in column)
        elif traceback_matrix[i, j] == "UP":
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2  
            i -= 1
        # LEFT = Gap in the first sequence (in line)
        else:
            aligned_seq1 = "-" + aligned_seq1  
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1

    return aligned_seq1, aligned_seq2, score_matrix[len(seq1), len(seq2)]

def pairwise_alignement_matrix(dict_sequences):
    """
    Compute a pairwise alignment score matrix for a set of sequences.

    Args:
        dict_sequences (dict): A dictionary where the keys are sequence identifiers and the values are the 
                               sequences (str) to be aligned.

    Returns:
        pd.DataFrame: A pandas DataFrame containing the pairwise alignment scores. 
                      The matrix is lower triangular since scores are symmetric.
    """
    # initialize alignement score dictionnary and the names of sequences 
    alignement_score_matrix = {}
    sequences = list(dict_sequences.keys())
    # Aligning pair of sequences
    for index, (key, value), in enumerate(dict_sequences.items()):
        # adding Nan values because the matrix will be symetric
        alignement_score_matrix[key] = [np.nan]*(len(dict_sequences) - index)
        if index > 0:
            # Reverse loop to get scores in the right order
            for i in range(index-1, -1, -1):
                a_seq1, a_seq2, score = two_seq_align(dict_sequences[key], dict_sequences[sequences[i]])
                alignement_score_matrix[key] = [score] + alignement_score_matrix[key]
    # Tranform the alignement score dictionnary to a dataframe
    alignement_score_matrix = pd.DataFrame(alignement_score_matrix, index = sequences)
    return alignement_score_matrix

def score_to_distance(score_matrix):
    """
    Convert a pairwise alignment score matrix into a normalized distance matrix.
    The scores are first normalized (between 0 and 1) based on the minimum and maximum values in the matrix '(score - minimum score)/(maximum score - minimum score),
    and then converted to distances by applying the formula `distance = 1 - normalized_score`. 
    This method ensures that higher alignment scores correspond to shorter distances.

    Args:
        score_matrix (pd.DataFrame): A pandas DataFrame representing the pairwise alignment score matrix, 
                                     where rows and columns correspond to sequence identifiers.

    Returns:
        pd.DataFrame: A pandas DataFrame representing the normalized distance matrix. 
        The distance values range between 0 and 1, with 0 indicating identical sequences and 1 indicating maximal distance.
    """
    normalized_matrix = score_matrix.copy() 
    global_min = normalized_matrix.min().min()
    global_max = normalized_matrix.max().max() 
    # Normalisation of scores
    if global_max != global_min:
        normalized_matrix = (normalized_matrix - global_min) / (global_max - global_min)
    # Transform to distance
    distance_matrix = 1 - normalized_matrix
    return distance_matrix

def extract_clusters(cluster):
    """
    Extract all leaf clusters from a hierarchical cluster structure by recursion.

    Args:
        cluster (tuple): A hierarchical cluster structure, which can be a nested tuple or a string (representing a sequence).

    Returns:
        list: A flat list containing all leaf clusters. Each leaf cluster is either a sequence (string) or a tuple of two sequences.
    """
    result = []

    def process(subcluster):
        if isinstance(subcluster, tuple):
            if len(subcluster) == 2 and all(isinstance(item, str) for item in subcluster):
                result.append(subcluster)
            else:
                for item in subcluster:
                    process(item)
        else:
            result.append(subcluster)
    process(cluster)
    return result

def upgma(distance_matrix):
    """
    Perform hierarchical clustering using the UPGMA (Unweighted Pair Group Method with Arithmetic Mean) algorithm.

    Args:
        distance_matrix (pd.DataFrame): A square pandas DataFrame where rows and columns represent sequences, 
                                        and the values represent pairwise distances between sequences.

    Returns:
        list: A list of tuples or individual sequences in the order they were clustered. Each tuple represents a pair of sequences or clusters that were merged at a given step.
    """
    # set counter for distances weight
    counter = 1
    while len(distance_matrix) > 1:
        # Get information on the minimum distance in the matrix
        min, row, col, row_idx, col_idx = get_min(distance_matrix)
        # Form the cluster
        cluster = (row, col)
        # Prepare the new line in the matrix for the new cluster
        new_row = {}
        # Starts with Nan because it's the distance of the cluster with it self.
        new_row[cluster] = [np.nan]
        # Get the rest of sequences with which the distances shoukd be re-calculated
        seqs = [seq for seq in distance_matrix.columns if seq not in cluster]
        for seq in seqs :
            # Ensure that the distance is retrieved from the filled part of the matrix.
            if row_idx < distance_matrix.columns.get_loc(seq) and col_idx < distance_matrix.columns.get_loc(seq) : 
                distance = calc_distance(distance_matrix.iloc[row_idx, distance_matrix.columns.get_loc(seq)], distance_matrix.iloc[col_idx, distance_matrix.columns.get_loc(seq)], counter)
            elif row_idx < distance_matrix.columns.get_loc(seq) and col_idx > distance_matrix.columns.get_loc(seq) :
                distance = calc_distance(distance_matrix.iloc[row_idx, distance_matrix.columns.get_loc(seq)], distance_matrix.iloc[distance_matrix.columns.get_loc(seq), col_idx], counter)
            elif row_idx > distance_matrix.columns.get_loc(seq) and col_idx < distance_matrix.columns.get_loc(seq) :
                distance = calc_distance(distance_matrix.iloc[distance_matrix.columns.get_loc(seq), row_idx], distance_matrix.iloc[col_idx, distance_matrix.columns.get_loc(seq)], counter)
            else :
                distance = calc_distance(distance_matrix.iloc[distance_matrix.columns.get_loc(seq), row_idx], distance_matrix.iloc[distance_matrix.columns.get_loc(seq), col_idx], counter)

            # add the distance to the new line
            new_row[cluster].append(distance)
            
        # Transforming the matrix
        # Drop the two grouped sequences from lines and columns
        distance_matrix.drop([row, col], axis = 0, inplace = True)
        distance_matrix.drop([row, col], axis = 1, inplace = True)
        # adding the new line and column
        # empty part of the matrix
        distance_matrix.insert(0, cluster, [np.nan]*len(distance_matrix))
        # full part of the matrix
        new_row_df = pd.DataFrame([new_row[cluster]], columns=distance_matrix.columns, index=[cluster])
        distance_matrix = pd.concat([new_row_df, distance_matrix])

        counter +=1
    # The final cluster is a tuple containing all the clusters formed during the previous steps
    # Extracting leafs of the clustering so we can iterate on them during multiple alignment
    order_seqs = extract_clusters(cluster)
    return order_seqs

def multi_seq_align(cluster1, cluster2):
    """
    Perform multiple sequence alignment between two clusters of sequences.

    Args:
        cluster1 (list of str): A list of sequences (strings) representing the first cluster of sequences to align.
        cluster2 (list of str): A list of sequence(s) (strings) representing the second cluster of sequences to align.

    Returns:
        tuple: A tuple containing two lists:
            - aligned_seqs1 (list of str): The aligned sequences from `cluster1`, where each sequence has been aligned with the sequences in `cluster2`.
            - aligned_seqs2 (list of str): The aligned sequences from `cluster2`, where each sequence has been aligned with the sequences in `cluster1`.
    """
    n = len(cluster1[0])
    m = len(cluster2[0])
    # Initaite the score matrix and the traceback_matrix
    score_matrix = np.zeros((n+1, m+1))
    traceback_matrix = np.empty((n+1, m+1), dtype='U10')

    # Filling the first line
    for i in range(1, n+1):
        score_matrix[i, 0] += score_matrix[i-1, 0] + GAP_PENALTY
        traceback_matrix [i, 0] = "UP"
    #Filling the first colomun
    for j in range(1, m+1):
        score_matrix[0, j] += score_matrix[0, j-1] + GAP_PENALTY
        traceback_matrix [0, j] = "LEFT"

    #Filling the score matrix
    for i in range(n):
        for j in range(m):
            aa_cluster1 = get_aa(cluster1, i)
            aa_cluster2 = get_aa(cluster2, j)
            match = score_matrix[i, j] + calc_multi_score(aa_cluster1, aa_cluster2)
            delete = score_matrix[i, j+1] + GAP_PENALTY
            insert = score_matrix[i+1, j] + GAP_PENALTY

            score_matrix[i+1, j+1] = max(match, delete, insert)

            # Fiiling the traceback_matrix
            if max(match, delete, insert) == match:
                traceback_matrix[i+1][j+1] = "DIAG" 
            elif max(match, delete, insert) == delete:
                traceback_matrix[i+1][j+1] = "UP"
            else:
                traceback_matrix[i+1][j+1] = "LEFT"

    # Initiate aligned clusters, containing the futur aligned sequences for each cluster
    aligned_seqs1 = []
    aligned_seqs2 = []
    for i in range(len(cluster1)):
        aligned_seqs1.append("")
    for i in range(len(cluster2)):
        aligned_seqs2.append("")

    # Trace alignment from last cell in matrix
    i = n
    j = m

    while i>0 or j>0:
        # Matching amino acids
        if traceback_matrix[i, j] == "DIAG":
            for k in range(len(aligned_seqs1)):
                aligned_seqs1[k] = cluster1[k][i-1] + aligned_seqs1[k]
            for k in range(len(aligned_seqs2)):
                aligned_seqs2[k] = cluster2[k][j-1] + aligned_seqs2[k]
            i -= 1
            j -= 1
        # UP = Gap in the second sequence (in column)
        elif traceback_matrix[i, j] == "UP":
            for k in range(len(aligned_seqs1)):
                aligned_seqs1[k] = cluster1[k][i-1] + aligned_seqs1[k]
            for k in range(len(aligned_seqs2)):
                aligned_seqs2[k] = "-" + aligned_seqs2[k]
            i -= 1
        # LEFT = Gap in the first sequence (in line)
        else:
            for k in range(len(aligned_seqs1)):
                aligned_seqs1[k] = "-" + aligned_seqs1[k]
            for k in range(len(aligned_seqs2)):
                aligned_seqs2[k] = cluster2[k][j-1] + aligned_seqs2[k]
            j -= 1
    # Returning aligned clusters
    return aligned_seqs1, aligned_seqs2

def count_non_dash_chars(sequence):
    """
    Count the number of characters in a sequence that are not dashes '-' (gap).

    Args:
        sequence (str): The sequence string.

    Returns:
        int: The count of characters in the sequence that are not dashes.
    """
    return sum(1 for char in sequence if char != '-')

def print_aligned_sequences(sequence_names, aligned_sequences, line_length=70):
    """
    Prints aligned sequences in a formated manner

    Args:
        sequence_names (list): a list of sequence names in the alignement order
        aligned_sequences (list): a list containing aligned sequences in the alignement order
        line_length (int): Number of caracters printed per line

        Retruns:
            None: Prints the aligned sequences
    """
    # Verifying that all sequences have a name and vice versa
    if len(sequence_names) != len(aligned_sequences):
        raise ValueError("Les listes de noms de séquences et de séquences alignées doivent avoir la même longueur")
    print("\n")
    num_chars_list = [0] * len(sequence_names)
    for i in range(0, len(aligned_sequences[0]), line_length):
            k = 0
            for name, seq in zip(sequence_names, aligned_sequences):
                # Sequence name
                print(f"{name:<30}", end="")
                # Sequence
                print(seq[i:i+line_length], end=" ")
                # The number of amino acids in the line
                num_chars = count_non_dash_chars(seq[i:i+line_length]) + num_chars_list[k]
                num_chars_list[k] = num_chars
                k += 1
                print(f"{num_chars}") 
            print("\n")

def clust_align(fasta_file):
    """
    Perform alignement of sequences provided.
    This function takes a FASTA file as input,if more than 2 sequences given, performs pairwise alignments to build a score matrix, 
    converts the score matrix to a distance matrix, and then uses UPGMA to generate a hierarchical clustering of sequences. 
    It then performs multiple sequence alignments following the clustering results and prints the aligned sequences in a formatted manner.
    If only to sequences given, a two sequence alignement is performed directly.

    Args:
        fasta_file (str): The path to the FASTA file containing the sequences to be aligned.

    Returns:
        None: Prints final aligned sequences to the console.
    """
    # Get sequences to align
    sequences = read_fasta(fasta_file)
    if len(sequences) == 1:
        print("Only one sequence provided")
        return
    elif len(sequences) == 2 :
        aligned_seq1, aligned_seq2, score = two_seq_align(sequences[list(sequences.keys())[0]], sequences[list(sequences.keys())[1]])
        print_aligned_sequences(list(sequences.keys()), [aligned_seq1, aligned_seq2])
        return
    else :
        # Generate the score matrix
        score_matrix = pairwise_alignement_matrix(sequences)
        # Generate the distance matrix
        distance_matrix = score_to_distance(score_matrix)
        # Perform UPGMA clustering
        ordered_sequences = upgma(distance_matrix)
        # Perform multiple alignement
        aligned_sequences_name = [ordered_sequences[0][0], ordered_sequences[0][1]]
        #First alignement
        aligned_seq1, aligned_seq2, score = two_seq_align(sequences[ordered_sequences[0][0]], sequences[ordered_sequences[0][1]])
        base_cluster = [aligned_seq1, aligned_seq2]

        for cluster in ordered_sequences[1:] :
            # If the cluster contains two sequences, align them together before aligning to the previous cluster
            if len(cluster) == 2:
                aligned_seq1, aligned_seq2, score = two_seq_align(sequences[cluster[0]], sequences[cluster[1]])
                aligned_cluster1, aligned_cluster2 = multi_seq_align(base_cluster, [aligned_seq2, aligned_seq2])
                base_cluster = aligned_cluster1 + aligned_cluster2
                aligned_sequences_name.append(cluster[0])
                aligned_sequences_name.append(cluster[1])
            # If only one sequence represent the cluster, align it directly to the cluster
            else:
                aligned_cluster1, aligned_cluster2 = multi_seq_align(base_cluster, [sequences[cluster]])
                base_cluster = aligned_cluster1 + aligned_cluster2
                aligned_sequences_name.append(cluster)
        
        # Show the alignement
        print_aligned_sequences(aligned_sequences_name, base_cluster)


if __name__== "__main__":
    if len(sys.argv) != 2:
        sys.exit("ERROR: Exactly one argument is required.")
    
    file_path = sys.argv[1]

    if not os.path.exists(file_path):
        sys.exit(f"ERROR: The file '{file_path}' does not exist.")

    execution_time = timeit.timeit("clust_align(file_path)", globals=globals(), number=1)
    print(f"Execution time: {execution_time:.2f} seconds")