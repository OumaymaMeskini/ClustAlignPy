import pandas as pd
import numpy as np

DISTANCE_MATRIX = {'Hu': [np.nan, np.nan, np.nan, np.nan, np.nan],
        'Ch': [15, np.nan, np.nan, np.nan, np.nan],
        'Go': [45, 30, np.nan, np.nan, np.nan],
        'Or': [143, 126, 92, np.nan, np.nan],
        'Gi': [198, 179, 179, 179, np.nan]}

def calc_distance(distance1, distance2, counter):
    """Calculate distance between two clusters.

    Parameters
    ----------
    distance1 : float
    distance2 : float
    counter : int

    Returns
    -------
    float
        a distance between the two clusters
    """
    distance = ((counter * distance1 + distance2) / (counter + 1))
    return distance

def get_min(distance_matrix):
    """Get the minimum distance in a distance matrix

    Parameters
    ----------
    distance_matrix : Dataframe

    Retruns
    -------
    float
        a distance between two clusters
    """
    col_min = distance_matrix.min().idxmin()
    row_min = distance_matrix[col_min].idxmin()
    # Trouver les indices de la ligne et de la colonne du minimum
    min_row_index = distance_matrix.index.get_loc(row_min)
    min_col_index = distance_matrix.columns.get_loc(col_min)
    return row_min, col_min, min_row_index, min_col_index

def upgma(distance_matrix):
    counter = 1
    while len(distance_matrix) > 1:
        row, col, row_idx, col_idx = get_min(distance_matrix)
        cluster = (row, col)
        new_row = {}
        new_row[cluster] = [np.nan]
        seqs = [seq for seq in distance_matrix.columns if seq not in cluster]
        # Vérifier que la distance est récupérer dans la partie pleine de la matrice
        for seq in seqs :
            if row_idx < distance_matrix.columns.get_loc(seq) and col_idx < distance_matrix.columns.get_loc(seq) : 
                distance = calc_distance(distance_matrix.iloc[row_idx, distance_matrix.columns.get_loc(seq)], distance_matrix.iloc[col_idx, distance_matrix.columns.get_loc(seq)], counter)
            elif row_idx < distance_matrix.columns.get_loc(seq) and col_idx > distance_matrix.columns.get_loc(seq) :
                distance = calc_distance(distance_matrix.iloc[row_idx, distance_matrix.columns.get_loc(seq)], distance_matrix.iloc[distance_matrix.columns.get_loc(seq), col_idx], counter)
            elif row_idx > distance_matrix.columns.get_loc(seq) and col_idx < distance_matrix.columns.get_loc(seq) :
                distance = calc_distance(distance_matrix.iloc[distance_matrix.columns.get_loc(seq), row_idx], distance_matrix.iloc[col_idx, distance_matrix.columns.get_loc(seq)], counter)
            else :
                distance = calc_distance(distance_matrix.iloc[distance_matrix.columns.get_loc(seq), row_idx], distance_matrix.iloc[distance_matrix.columns.get_loc(seq), col_idx], counter)
            new_row[cluster].append(distance)
            
        print(new_row)
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

    #return distance_matrix

    

if __name__ == "__main__":
    test = pd.DataFrame(DISTANCE_MATRIX, index=['Hu', 'Ch', 'Go', 'Or', 'Gi'])
    test2 = pd.read_csv("./matrice_dist_UPGMA.csv")
    test2 = pd.DataFrame(test2, index=['Bsu', 'Bst', 'Lvi', 'Amo', 'Mlu'])
    upgma(test2)