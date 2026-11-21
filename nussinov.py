import numpy as np
import pandas as pd

MIN_LOOP_LENGTH = 3


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



def fill_matrix(matrix, rna, min_loop_length=MIN_LOOP_LENGTH):
    n = len(rna)

    for k in range(1, n):
        for i in range(n - k):
            j = i + k

            # Case 1: nucleotide i remains unpaired
            down = matrix[i + 1][j]

            # Case 2: nucleotide j remains unpaired
            left = matrix[i][j - 1]

            # Case 3: nucleotides i and j form a pair
            diag = matrix[i + 1][j - 1]

            best = max(down, left, diag)

            if (has_connection((rna[i], rna[j])) and j - i > min_loop_length):
                best = max(best, diag + 1)

            # Case 4: split the sequence into two subsequences
            for split in range(i, j):
                score = (
                    matrix[i][split]
                    + matrix[split + 1][j]
                )
                best = max(best, score)

            matrix[i][j] = best

    return matrix


def traceback(matrix, rna, i, j, pairs, min_loop_length=MIN_LOOP_LENGTH):
    if i >= j:
        return

    # Case 1: nucleotide i remains unpaired
    if matrix[i][j] == matrix[i + 1][j]:
        traceback(matrix, rna, i + 1, j, pairs, min_loop_length)

    # Case 2: nucleotide j remains unpaired
    elif matrix[i][j] == matrix[i][j - 1]:
        traceback(matrix, rna, i, j - 1, pairs, min_loop_length)

    # Case 3: nucleotides i and j form a pair
    elif (has_connection((rna[i], rna[j])) and j - i > min_loop_length and matrix[i][j] == matrix[i + 1][j - 1] + 1):
        pairs.append((i, j))
        traceback(matrix, rna, i + 1, j - 1, pairs, min_loop_length)

    # Case 4: split into two subsequences
    else:
        for k in range(i, j):
            if matrix[i][j] == (
                matrix[i][k] + matrix[k + 1][j]
            ):
                traceback(matrix, rna, i, k, pairs, min_loop_length)
                traceback(matrix, rna, k + 1, j, pairs, min_loop_length)
                break

    return pairs


rna = "GGGGUGACCC"
result = init_matrix(rna)
result = fill_matrix(result, rna)

nucliodites = [nucliotide for nucliotide in rna]

pairs = traceback(result, rna, 0, len(rna) - 1, [])

print("\nOptimal base pairs:")
for i, j in sorted(pairs):
    print(f"{rna[i]}{i + 1} - {rna[j]}{j + 1}")

result = pd.DataFrame(result)
result.columns = nucliodites
result.index = nucliodites
print(result)