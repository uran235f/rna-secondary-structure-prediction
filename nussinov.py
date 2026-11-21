import numpy as np
import pandas as pd


def init_matrix(sequence):
    size = len(sequence)
    matrix = np.full((size, size), 0)
    return matrix


def has_connection(pair) -> bool:
    connections_list = [("A", "U"),
                        ("U", "A"),
                        ("G", "C"),
                        ("C", "G"),
                        ("G", "U"),
                        ("U", "G")]
    if pair in connections_list:
        return True
    return False


def fill_matrix(matrix, rna):
    for k in range(1, len(rna)):
        for i in range(len(rna) - k):
            j = i + k
            down = matrix[i+1][j]
            left = matrix[i][j-1]
            diag = matrix[i+1][j-1]

            if has_connection((rna[i], rna[j])):
                matrix[i][j] = diag + 1
            else:
                matrix[i][j] = max(down, left, diag)
    return matrix

rna = "GGGGUGACCC"
result = init_matrix(rna)
result = fill_matrix(result, rna)

nucliodites = [nucliotide for nucliotide in rna]

result = pd.DataFrame(result)
result.columns = nucliodites
result.index = nucliodites
print(result)