
import pandas as pd
import numpy as np


def global_msa_sub(align1, align2, subs_mat, gap_penalty):
    '''Global alignment method aka Needleman-Wunsch alignment generalised for MSA'''

    # alignments come with list of lists
    # create two matrices (one for storing the values and the other for traceback)
    n_align1 = len(align1[0])
    n_align2 = len(align2[0])

    nrow = len(align1) + 1
    ncol = len(align2) + 1

    # generate row and column names
    def generate_names(align):
        names = ['']
        names_sub = [''.join(el) for el in align]
        names.extend(names_sub)
        return names, names_sub

    # we need string representation of lists that are inside list
    rownames, align1_str = generate_names(align1)
    colnames, align2_str = generate_names(align2)

    # initialize pandas df
    df1 = pd.DataFrame(np.zeros((nrow, ncol)), index=rownames, columns=colnames)

    # calculate substitution score (each letter against each letter in an alignment position)
    def calculate_subs_score(align1_pos, align2_pos, subs_mat):
        subs_score = 0
        for el1 in align1_pos:
            for el2 in align2_pos:
                if el1 == '_' or el2 == '_':
                    subs_score += gap_penalty
                else:
                    subs_score += subs_mat.loc[el1, el2]
        return subs_score / (len(align1_pos) * len(align2_pos))

    # fill the first row and the first column
    df1.iloc[0, 0] = 0
    for i in range(1, nrow):
        df1.iloc[i, 0] = df1.iloc[i - 1, 0] + gap_penalty
    for j in range(1, ncol):
        df1.iloc[0, j] = df1.iloc[0, j - 1] + gap_penalty

    # fill the rest of the table
    for j in range(1, ncol):
        for i in range(1, nrow):
            subs_score = calculate_subs_score(align1[i - 1], align2[j - 1], subs_mat)

            df1.iloc[i, j] = max(df1.iloc[i - 1, j] + gap_penalty,
                                 df1.iloc[i, j - 1] + gap_penalty,
                                 df1.iloc[i - 1, j - 1] + subs_score)

    final_score = df1.iloc[nrow - 1, ncol - 1]

    # reconstruct individual sequences
    def traceback(df, align1, align2):

        reconst_align1 = [align1_str[-1]]
        reconst_align2 = [align2_str[-1]]

        # i for row number and j for column number
        j = df.shape[1] - 1
        i = df.shape[0] - 1

        while i > 1 and j > 1:

            # if score in the diagonal is the highest
            if max(df.iloc[i - 1, j - 1], df.iloc[i - 1, j], df.iloc[i, j - 1]) == df.iloc[i - 1, j - 1]:
                reconst_align1.append(align1_str[i - 2])
                reconst_align2.append(align2_str[j - 2])
                i -= 1
                j -= 1

                # if score in the left is the highest
            elif max(df.iloc[i - 1, j - 1], df.iloc[i - 1, j - 1], df.iloc[i, j - 1]) == df.iloc[i, j - 1]:
                reconst_align1.append('_' * n_align1)
                reconst_align2.append(align2_str[j - 2])
                j -= 1

            # if score above is the highest
            else:
                reconst_align1.append(align1_str[i - 2])
                reconst_align2.append('_' * n_align2)
                i -= 1

        # print('\n'.join(get_individual_seqs(reconst_align1, reconst_align2)))
        return reconst_align1[::-1], reconst_align2[::-1]

    def combine_reconst_align(reconst1, reconst2):
        combined = []
        for i in range(len(reconst1)):
            combined_sub = []
            for el in reconst1[i]:
                combined_sub.append(el)
            for el in reconst2[i]:
                combined_sub.append(el)
            combined.append(combined_sub)
        return combined

    reconst1, reconst2 = traceback(df1, align1, align2)
    combined = combine_reconst_align(reconst1, reconst2)

    return combined


def progressive_msa(sequences, aligner, guide_tree, subs_mat, seq_names, gap_penalty):
    node1 = guide_tree.get_left()
    node2 = guide_tree.get_right()
    if node1.is_leaf():
        node1_align = sequences[seq_names[node1.get_id()]]
    else:
        node1_align = progressive_msa(sequences, aligner, node1, subs_mat, seq_names, gap_penalty)
    if node2.is_leaf():
        node2_align = sequences[seq_names[node2.get_id()]]
    else:
        node2_align = progressive_msa(sequences, aligner, node2, subs_mat, seq_names, gap_penalty)
    alignment = aligner(node1_align, node2_align, subs_mat, gap_penalty)
    return alignment