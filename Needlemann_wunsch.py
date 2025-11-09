def needleman_wunsch(seq1, seq2, match=1, mismatch=-1, gap=-1):
    """
    Strict global Needleman–Wunsch.
    - First row and first column are initialized with cumulative gap penalties (global).
    - Tie-breaking prefers diagonal when equal (standard).
    - Traceback always starts at bottom-right (n,m).
    Returns: align1, align2, final_score, score_matrix, traceback_list
    """
    n = len(seq2)   # rows
    m = len(seq1)   # cols

    # --- Initialization ---
    score_matrix = [[0 for _ in range(m + 1)] for _ in range(n + 1)]
    pointer_matrix = [["" for _ in range(m + 1)] for _ in range(n + 1)]

    # Global initialization: cumulative gap penalties on first row & column
    for i in range(1, n + 1):
        score_matrix[i][0] = score_matrix[i - 1][0] + gap
        pointer_matrix[i][0] = 'U'   # came from up (gap in seq1)
    for j in range(1, m + 1):
        score_matrix[0][j] = score_matrix[0][j - 1] + gap
        pointer_matrix[0][j] = 'L'   # came from left (gap in seq2)

    pointer_matrix[0][0] = 'S'  # start/origin

    # --- Fill the matrix (standard recurrence) ---
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            mm = match if seq2[i - 1] == seq1[j - 1] else mismatch
            diag_score = score_matrix[i - 1][j - 1] + mm
            up_score = score_matrix[i - 1][j] + gap
            left_score = score_matrix[i][j - 1] + gap

            # Standard tie-handling: prefer diagonal if equal (diag >= up, left)
            if diag_score >= up_score and diag_score >= left_score:
                score_matrix[i][j] = diag_score
                pointer_matrix[i][j] = 'D'
            elif up_score >= left_score:
                score_matrix[i][j] = up_score
                pointer_matrix[i][j] = 'U'
            else:
                score_matrix[i][j] = left_score
                pointer_matrix[i][j] = 'L'

    # --- Traceback: start at bottom-right for strict global alignment ---
    i, j = n, m
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

        # Defensive fallback: this should not happen in strict global, but handle gracefully
        if i > 0 and j > 0:
            align1 = seq1[j - 1] + align1
            align2 = seq2[i - 1] + align2
            traceback.append(f"Mismatch(FO): {seq1[j - 1]} with {seq2[i - 1]}")
            i -= 1
            j -= 1
            continue
        elif j > 0:
            align1 = seq1[j - 1] + align1
            align2 = "-" + align2
            traceback.append(f"End-gap (Seq2): {seq1[j - 1]} with -")
            j -= 1
            continue
        elif i > 0:
            align1 = "-" + align1
            align2 = seq2[i - 1] + align2
            traceback.append(f"End-gap (Seq1): - with {seq2[i - 1]}")
            i -= 1
            continue
        else:
            break

    traceback.reverse()
    final_score = score_matrix[n][m]
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