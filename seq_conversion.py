
import numpy as np


# this script converts protein sequence to dna based on human codon usage
# codon usage (probabilities) are taken from https://www.genscript.com/tools/codon-frequency-table

def convert_to_prot(seq):

    '''converts DNA intp protein sequence'''
    # it may chop off some dna letters if n%3 !=0

    conversion_table_DNA_to_prot = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '-', 'TAG': '-',
        'TGC': 'C', 'TGT': 'C', 'TGA': '-', 'TGG': 'W',
    }

    # '_' denotes stop codon

    protein = []
    for i in range(0, len(seq), 3):
        codon = seq[i:i + 3]
        aa = conversion_table_DNA_to_prot[codon]

        # if stop codon encountered then break the cycle
        if aa == '-':
            break

        protein.append(aa)

    return ''.join(protein)




def convert_to_dna(seq):

    ''' Converts DNA to protein sequence (based on codon useage probabilities)
    probabiliteis are taken from https://www.genscript.com/tools/codon-frequency-table '''

    conversion_table_prot_to_DNA = {
        'A': (['GCT', 'GCC', 'GCA', 'GCG'], [0.26, 0.40, 0.23, 0.11]),
        'R': (['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], [0.09, 0.19, 0.11, 0.21, 0.20, 0.20]),
        'N': (['AAT', 'AAC'], [0.46, 0.54]),
        'D': (['GAT', 'GAC'], [0.46, 0.54]),
        'C': (['TGT', 'TGC'], [0.45, 0.55]),
        'Q': (['CAA', 'CAG'], [0.25, 0.75]),
        'E': (['GAA', 'GAG'], [0.42, 0.58]),
        'G': (['GGT', 'GGC', 'GGA', 'GGG'], [0.16, 0.34, 0.25, 0.25]),
        'H': (['CAT', 'CAC'], [0.41, 0.59]),
        'I': (['ATT', 'ATC', 'ATA'], [0.36, 0.48, 0.16]),
        'L': (['TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'], [0.07, 0.13, 0.13, 0.20, 0.07, 0.40]),
        'K': (['AAA', 'AAG'], [0.42, 0.58]),
        'M': (['ATG'], [1]),
        'F': (['TTC', 'TTT'], [0.55, 0.45]),
        'P': (['CCT', 'CCC', 'CCA', 'CCG'], [0.28, 0.33, 0.28, 0.11]),
        'S': (['TCT', 'TCC', 'TCA', 'TCG', 'AGC', 'AGT'], [0.18, 0.22, 0.15, 0.06, 0.15, 0.24]),
        'T': (['ACT', 'ACC', 'ACA', 'ACG'], [0.24, 0.36, 0.28, 0.12]),
        'W': (['TGG'], [1]),
        'Y': (['TAC', 'TAT'], [0.57, 0.43]),
        'V': (['GTT', 'GTC', 'GTA', 'GTG'], [0.18, 0.24, 0.11, 0.47]),
    }

    converted_seq = []
    for aa in seq:
        triplets, prob = conversion_table_prot_to_DNA[aa]
        dna_triplet = np.random.choice(triplets, p=prob)
        converted_seq.append(dna_triplet)

    return ''.join(converted_seq)

