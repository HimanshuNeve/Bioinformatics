def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-1):
   
    n = len(seq2)
    m = len(seq1)

    # --- Initialization ---
    score_matrix = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    pointer_matrix = [["" for _ in range(m + 1)] for _ in range(n + 1)]

    # first row/col zeros
    for i in range(1, n + 1):
        score_matrix[i][0] = 0
        pointer_matrix[i][0] = ''
    for j in range(1, m + 1):
        score_matrix[0][j] = 0
        pointer_matrix[0][j] = ''

    pointer_matrix[0][0] = 'S'

    # --- Fill the matrix ---
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            mm = match if seq2[i - 1] == seq1[j - 1] else mismatch
            diag_score = score_matrix[i - 1][j - 1] + mm
            up_score = score_matrix[i - 1][j] + gap
            left_score = score_matrix[i][j - 1] + gap

            if diag_score >= up_score and diag_score >= left_score:
                score_matrix[i][j] = diag_score
                pointer_matrix[i][j] = 'D'
            elif up_score >= left_score:
                score_matrix[i][j] = up_score
                pointer_matrix[i][j] = 'U'
            else:
                score_matrix[i][j] = left_score
                pointer_matrix[i][j] = 'L'

    # --- Choose traceback start to allow trailing end-gaps free ---
    last_row_max = max(score_matrix[n])
    last_col_max = max(score_matrix[i][m] for i in range(n + 1))
    final_score = last_row_max if last_row_max >= last_col_max else last_col_max

    if score_matrix[n][m] == final_score:
        i, j = n, m
    else:
        found = False
        for j_candidate in range(m, -1, -1):
            if score_matrix[n][j_candidate] == final_score:
                i, j = n, j_candidate
                found = True
                break
        if not found:
            for i_candidate in range(n, -1, -1):
                if score_matrix[i_candidate][m] == final_score:
                    i, j = i_candidate, m
                    found = True
                    break
        if not found:
            i, j = n, m

    # Record start position — needed to append trailing characters later
    start_i, start_j = i, j

    # --- Traceback (use pointer matrix but also be robust) ---
    align1 = ""
    align2 = ""
    traceback = []

    while i > 0 or j > 0:
        move = pointer_matrix[i][j] if 0 <= i <= n and 0 <= j <= m else ""

        if move == 'D':
            align1 = seq1[j - 1] + align1
            align2 = seq2[i - 1] + align2
            if seq1[j - 1] == seq2[i - 1]:
                traceback.append(f"Match: {seq1[j - 1]} with {seq2[i - 1]}")
            else:
                traceback.append(f"Mismatch: {seq1[j - 1]} with {seq2[i - 1]}")
            i -= 1
            j -= 1
            continue

        if move == 'U' and i > 0:
            align1 = "-" + align1
            align2 = seq2[i - 1] + align2
            traceback.append(f"Gap in Seq1: {seq2[i - 1]} with -")
            i -= 1
            continue

        if move == 'L' and j > 0:
            align1 = seq1[j - 1] + align1
            align2 = "-" + align2
            traceback.append(f"Gap in Seq2: {seq1[j - 1]} with -")
            j -= 1
            continue

        # Boundary: insert remaining end-gaps until indices consumed
        if j > 0 and i == 0:
            align1 = seq1[j - 1] + align1
            align2 = "-" + align2
            traceback.append(f"End-gap (Seq2): {seq1[j - 1]} with -")
            j -= 1
            continue
        if i > 0 and j == 0:
            align1 = "-" + align1
            align2 = seq2[i - 1] + align2
            traceback.append(f"End-gap (Seq1): - with {seq2[i - 1]}")
            i -= 1
            continue

        # Fallback if pointers missing: consume diagonal
        if i > 0 and j > 0:
            align1 = seq1[j - 1] + align1
            align2 = seq2[i - 1] + align2
            if seq1[j - 1] == seq2[i - 1]:
                traceback.append(f"Match(FO): {seq1[j - 1]} with {seq2[i - 1]}")
            else:
                traceback.append(f"Mismatch(FO): {seq1[j - 1]} with {seq2[i - 1]}")
            i -= 1
            j -= 1
            continue

        break

    # Reverse traceback for human-readable order
    traceback.reverse()

    # --- APPEND TRAILING CHARACTERS AS END-GAPS IF TRACEBACK STARTED EARLIER ---
    # If we started at (start_i, start_j) not equal to (n,m), append any characters
    # to the right of the alignment that were skipped — these are trailing chars.
    # If start_j < m -> seq1 has trailing chars that must be added, paired with gaps in seq2.
    if start_j < m:
        # chars in seq1 from index start_j (0-based: start_j ... m-1)
        tail = seq1[start_j:m]
        align1 = align1 + tail
        align2 = align2 + ("-" * len(tail))
        # add messages for clarity (optional)
        for ch in tail:
            traceback.append(f"Appended tail (Seq1): {ch} with -")

    # If start_i < n -> seq2 has trailing chars that must be added, paired with gaps in seq1.
    if start_i < n:
        tail = seq2[start_i:n]
        align1 = align1 + ("-" * len(tail))
        align2 = align2 + tail
        for ch in tail:
            traceback.append(f"Appended tail (Seq2): - with {ch}")

    # final_score already computed above (final_score variable)
    return align1, align2, final_score, score_matrix, traceback


# print the matrix 
def print_matrix(matrix, seq1, seq2):
    max_val_width = 1 
    for row in matrix:
        for val in row:
            val_width = len(str(val))
            if val_width > max_val_width:
                max_val_width = val_width
    
    col_width = max_val_width + 2

   
    print(" " * 5, end="") 
    print(f"{'':>{col_width}}", end="")
    for char in seq1:
        print(f"{char:>{col_width}}", end="")
    print()

    side_chars = " " + seq2
    for i in range(len(matrix)):
        print(f"  {side_chars[i]}  ", end="") 
        for val in matrix[i]:
            print(f"{val:>{col_width}}", end="")
        print() 


if __name__ == "__main__":
    
    seq_1 = input("Enter first sequence: ").upper()
    seq_2 = input("Enter second sequence: ").upper()
    
    match_score = 1
    mismatch_score = -2
    gap_score = -4
    
    print(f"\nSeq 1: {seq_1}")
    print(f"Seq 2: {seq_2}")
    
    
    alignment_1, alignment_2, final_score, matrix_to_print, traceback = needleman_wunsch(
        seq_1, 
        seq_2, 
        match_score, 
        mismatch_score, 
        gap_score
    )
    
 
    print("\n--- Final Score Matrix ---")
    print_matrix(matrix_to_print, seq_1, seq_2)

    print("\n--- Traceback Steps ---")
    print("\n- Seq1 with Seq2 -")
    for step in traceback:
        print(step)
    
    print("\n--- Results ---")
    print(f"Alignment 1: {alignment_1}")
    print(f"Alignment 2: {alignment_2}")
    print(f"Final Score: {final_score}")