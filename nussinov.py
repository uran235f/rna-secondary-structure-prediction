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


def to_dot_bracket(rna, pairs):
    structure = ["."] * len(rna)

    for i, j in pairs:
        if not (0 <= i < j < len(rna)):
            raise ValueError(f"Invalid base pair: ({i}, {j})")

        if structure[i] != "." or structure[j] != ".":
            raise ValueError("A nucleotide cannot belong to multiple pairs.")

        structure[i] = "("
        structure[j] = ")"

    return "".join(structure)


def main():
    # Input RNA sequence
    rna = "GGGGUGACCC"

    # Calculate the dynamic programming matrix
    matrix = init_matrix(rna)
    matrix = fill_matrix(matrix, rna)

    # Reconstruct the optimal base pairs
    pairs = traceback(matrix, rna, 0, len(rna) - 1, [])

    # Generate dot-bracket representation
    dot_bracket = to_dot_bracket(rna, pairs)

    # Calculate prediction summary
    n = len(rna)
    max_pairs = int(matrix[0][n - 1])

    # Display results
    print("=" * 50)
    print("       RNA SECONDARY STRUCTURE PREDICTION")
    print("=" * 50)

    print("\n1. RNA SEQUENCE")
    print(f"   {rna}")
    print(f"   Length: {n} nucleotides")

    print("\n2. PREDICTED SECONDARY STRUCTURE")
    print(f"   {dot_bracket}")

    print("\n3. PREDICTION SUMMARY")
    print(f"   Maximum number of base pairs: {max_pairs}")
    print(f"   Unpaired nucleotides: {n - 2 * len(pairs)}")

    print("\n4. OPTIMAL BASE PAIRS")
    if pairs:
        print(f"   {'Position 1':<15}{'Position 2':<15}{'Pair':<10}")
        print("   " + "-" * 38)

        for i, j in sorted(pairs):
            pair = f"{rna[i]}-{rna[j]}"
            print(f"   {i + 1:<15}{j + 1:<15}{pair:<10}")
    else:
        print("   No base pairs found.")

    print("\n5. DYNAMIC PROGRAMMING MATRIX")
    labels = [
        f"{i + 1}:{nucleotide}"
        for i, nucleotide in enumerate(rna)
    ]

    matrix_df = pd.DataFrame(matrix, index=labels, columns=labels)
    print(matrix_df.to_string())
    
    print("\n" + "=" * 50)


if __name__ == "__main__":
    main()
