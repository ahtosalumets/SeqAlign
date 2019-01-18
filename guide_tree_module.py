import os
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import average, dendrogram, to_tree

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from time import gmtime, strftime


# guide tree generation given tree will guide
# the alignment of sequences and alignments
def generate_guide_tree(seqs, k=3):
    # calculate kmer composition
    def kmer(seq, k):
        return set([seq[i:i + k] for i in range(len(seq) - (k - 1))])

    # calculate distance between two sequences
    def kmer_dist(seq1, seq2, k=3):
        kmer1 = kmer(seq1, k)
        kmer2 = kmer(seq2, k)
        all_kmers = kmer1 | kmer2
        shared_kmers = kmer1 & kmer2
        uniq_kmers = len(all_kmers) - len(shared_kmers)
        return uniq_kmers / len(all_kmers)

    # generate distance matrix
    
    def dist_matrix(seqs, k=3):
        names = list(seqs.keys())
        n = len(names)
        dm = pd.DataFrame(np.zeros((n, n)), index=names, columns=names)
        for j in range(n):
            for i in range(n):
                dm.iloc[i, j] = kmer_dist(seqs[names[i]], seqs[names[j]], k)
        return dm

    # the first step is to generate distance matrix
    # the distance is defined as the fraction of unique
    # kmers between two sequences
    dm = dist_matrix(seqs, k)

    

    # make condensed distance matrix
    def condensed_dm(dm):
        cdm = []
        n = dm.shape[0]
        for j in range(n - 1):
            for i in range(j + 1, n):
                cdm.append(dm.iloc[i, j])
        return np.array(cdm)

    # the second step is to form a condensed distance matrix
    # i.e. 1 dimensional distance matrix. E.g. let's take 3 sequences
    # then the first element is a distance between seq1 and seq2,
    # second is between seq1 and seq3 and third one is distance
    # between seq2 and seq3
    cdm = condensed_dm(dm)
    

    def make_guide_tree(cdm, dm):
        guide = average(cdm)
        guide_tree = to_tree(guide)
        guide_d = dendrogram(guide, labels=list(dm.columns.values), orientation='right',
                             link_color_func=lambda x: 'black')
        current_time = strftime('%d%m%Y_%H%M%S', gmtime())
        if not os.path.isdir('./guide_tree_image'):
            os.makedirs('./guide_tree_image')
        plt.savefig('./guide_tree_image/guidetree_' + current_time + '.png')
        return guide_tree

    # third step is to make a tree based on
    # condensed distance matrix
    tree = make_guide_tree(cdm, dm)

    # traverse given tree in postorder and produce
    # a list of nodes
    def postorder(tree):
        ordered = []

        def postorder_sub(node):
            if node:
                postorder_sub(node.get_right())
                postorder_sub(node.get_left())
                if node.is_leaf():
                    ordered.append((node.get_id(), 'Leaf'))
                else:
                    ordered.append((node.get_id(), 'Node'))

        postorder_sub(tree)

        return ordered

    # generate postorder sequence of sequences
    ordered_seqs = []
    node_ind = postorder(tree)
    seq_names = list(dm.columns.values)
    for el in node_ind:
        if el[1] == 'Leaf':
            ordered_seqs.append((seq_names[el[0]], 'Leaf'))
        else:
            ordered_seqs.append(el)

    # in addition to the tree I also return column values
    # because in the tree the id corresponds to the position
    # in the dm.column.values which inturn helps to get
    # the name of the sequence. For debugging purposes
    # I also output the postorder structure (ordered_seqs)
    # of the given guide tree
    return tree, seq_names, ordered_seqs
