import sys
import numpy as np
import pandas as pd
from Bio.SubsMat import MatrixInfo as matrixFile

available_mat = [
 'benner6',
 'benner22',
 'benner74',
 'blosum30',
 'blosum35',
 'blosum40',
 'blosum45',
 'blosum50',
 'blosum55',
 'blosum60',
 'blosum62',
 'blosum65',
 'blosum70',
 'blosum75',
 'blosum80',
 'blosum85',
 'blosum90',
 'blosum95',
 'blosum100',
 'feng',
 'fitch',
 'genetic',
 'gonnet',
 'grant',
 'ident',
 'johnson',
 'levin',
 'mclach',
 'miyata',
 'nwsgappep',
 'pam30',
 'pam60',
 'pam90',
 'pam120',
 'pam180',
 'pam250',
 'pam300',
 'rao',
 'risler',
 'structure'
]


def select_substitution_mat(matrix):
    try:
        subs_matrix = eval('matrixFile.' + matrix)
    except AttributeError:
        print('Select different subsitution matrix:' + '\n'.join(available_mat))
    except:
        print("Unexpected error:", sys.exc_info()[0])
        raise

    def make_df(subs_matrix):
        # get individual aminoacids
        aminoacids = sorted(list(set([pair[0] for pair in subs_matrix.keys()])))

        # initialize matrix with zeros
        df = pd.DataFrame(np.zeros((len(aminoacids), len(aminoacids))),
                          index=aminoacids, columns=aminoacids)

        # fill the table with substitution penalties
        for pair in subs_matrix.keys():
            penalty = subs_matrix[pair]
            df.loc[pair[0], pair[1]] = penalty
            df.loc[pair[1], pair[0]] = penalty
        return df

    subs_mat = make_df(subs_matrix)
    return subs_mat


