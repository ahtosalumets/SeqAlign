import pandas as pd
import numpy as np

def global_alignment(seqs, subs_mat, gap_penalty, sequence_type):
    '''Global alignment method aka Needleman-Wunsch alignment'''

    # unpack dict
    names_l, seqs_l = zip(*seqs.items())
    seq1 = seqs_l[0]
    seq2 = seqs_l[1]

    # create two matrices (one for storing the values and the other for traceback)
    nrow = len(seq1) + 1
    ncol = len(seq2) + 1

    rownames = ['']
    rownames.extend(list(seq1))

    colnames = ['']
    colnames.extend(list(seq2))

    df1 = pd.DataFrame(np.zeros((nrow, ncol)), index=rownames, columns=colnames)
    df2 = pd.DataFrame(np.zeros((nrow, ncol)), index=rownames, columns=colnames)

    # fill the first row and the first column
    df1.iloc[0, 0] = 0
    for i in range(1, nrow):
        df1.iloc[i, 0] = df1.iloc[i - 1, 0] + gap_penalty
    for j in range(1, ncol):
        df1.iloc[0, j] = df1.iloc[0, j - 1] + gap_penalty

    # make trace df
    df2.iloc[0, 0] = '•'
    for i in range(1, nrow):
        df2.iloc[i, 0] = '↑'
    for j in range(1, ncol):
        df2.iloc[0, j] = '←'

    # fill the rest of the table
    for j in range(1, ncol):
        for i in range(1, nrow):
            if sequence_type == 'protein':
                subs_score = subs_mat.loc[seq1[i - 1], seq2[j - 1]]
            else:
                if seq1[i - 1] == seq2[j - 1]:
                    subs_score = 1
                else: subs_score = -1

            df1.iloc[i, j] = max(df1.iloc[i - 1, j] + gap_penalty,
                                 df1.iloc[i, j - 1] + gap_penalty,
                                 df1.iloc[i - 1, j - 1] + subs_score)

            if df1.iloc[i, j] == (df1.iloc[i - 1, j - 1] + subs_score):
                df2.iloc[i, j] = '↖'
            elif df1.iloc[i, j] == (df1.iloc[i - 1, j] + gap_penalty):
                df2.iloc[i, j] = '←'
            elif df1.iloc[i, j] == (df1.iloc[i, j - 1] + gap_penalty):
                df2.iloc[i, j] = '↑'

    final_score = df1.iloc[nrow - 1, ncol - 1]

    # print('{0} \n\n {1}'.format(df1, df2))

    # traceback function (function inside a function to avoid causing problems with a scope)
    def traceback(df, seq1, seq2, names_l):

        reconst_seq1 = [seq1[-1]]
        reconst_seq2 = [seq2[-1]]
        # i for row number and j for column number
        j = df.shape[1] - 1
        i = df.shape[0] - 1

        while i > 1 and j > 1:
            # if score in the diagonal is the highest
            if max(df.iloc[i - 1, j - 1], df.iloc[i - 1, j], df.iloc[i, j - 1]) == df.iloc[i - 1, j - 1]:
                reconst_seq1.append(seq1[i - 2])
                reconst_seq2.append(seq2[j - 2])
                i -= 1
                j -= 1

                # if score in the left is the highest
            elif max(df.iloc[i - 1, j - 1], df.iloc[i - 1, j - 1], df.iloc[i, j - 1]) == df.iloc[i, j - 1]:
                reconst_seq1.append('_')
                reconst_seq2.append(seq2[j - 2])
                j -= 1

            # if score above is the highest
            else:
                reconst_seq1.append(seq1[i - 2])
                reconst_seq2.append('_')
                i -= 1

        return '\n'.join([names_l[0] + ': ' + ''.join(reconst_seq1[::-1]),
                          names_l[1] + ': ' + ''.join(reconst_seq2[::-1])])

    # print('Sequence alignment score:', final_score, '\n')
    traceback_seqs = traceback(df1, seq1, seq2, names_l)
    return traceback_seqs, '{0} \n\n {1}'.format(df1, df2), final_score
