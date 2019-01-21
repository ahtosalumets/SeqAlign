import click
from pyfiglet import Figlet
import time
from functools import wraps
import substitution_matrix_module
import pairwise_global_module
import pairwise_local_module
import msa_global_module
import msa_local_module
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import guide_tree_module
import antibody_fr_extraction
import seq_conversion
import numpy as np


def timer(function):
    '''
    Running time measurement
    :param function
    :return running time and task result
    '''

    @wraps(function)
    def func_time(*args, **kwargs):
        t0 = time.clock()
        result = function(*args, **kwargs)
        running_time = time.clock() - t0
        # print("Total running time: %s s" % (str(running_time)))
        return [running_time, result]

    return func_time


def plotting(cost, title='Performance Analysis'):
    '''
    plotting function
    :param time cost
    :return output image
    '''

    plt.style.use('seaborn-whitegrid')
    fig = plt.figure(figsize=(10,8),dpi=100)
    ax = plt.axes()
    #plt.ylim(0, max(cost)*1.2)
    ax.plot(cost, color='red')
    plt.xticks(range(len(cost)))
    plt.xlabel("loop times", fontsize=18) 
    plt.ylabel("running time", fontsize=18) 
    plt.title(title, fontsize=18)
    plt.savefig("analysis.jpg")


def read_file(file, top):
    '''
    fasta file read function
    :param filename
    :return file
    '''

    seqs = dict()
    with open(file, 'r') as f:
        content = f.readlines()
    n_seqs = len(content) // 2
    if top < n_seqs:
        n_seqs = top
    for i in range(0, 2*n_seqs-1, 2):
        seq_name = content[i].strip().replace('>', '')
        seqs[seq_name] = content[i+1].strip()
    return seqs

def print_seqs(description, input_seqs):
    '''print dictionary of sequences in a nice way'''

    max_len_name = max(len(name) for name in input_seqs.keys())
    max_len_seq = max(len(seq) for seq in input_seqs.values())
    print(description)
    for k, v in input_seqs.items():
        print('{:>{max_name}} {:>{max_seq}}'.format(k, v, max_name=max_len_name, max_seq=max_len_seq))
    print()
    return

def replace_BZ(seq):
    '''replace not unique amino acid letters'''

    while 'B' in seq:
        seq = seq.replace('B', np.random.choice(['D', 'N'], p=[0.56, 0.44]), 1)
    while 'Z' in seq:
        seq = seq.replace('Z', np.random.choice(['Q', 'E'], p=[0.43, 0.57]), 1)
    return seq


def pretty_print_msa(alignment, sequences, alignment_type, ordered):
    '''
    result output function
    :param alignment result, original sequence list, alignment type, ordered sequence names
    :return
    '''
    sequences_inv = {v : k for k, v in sequences.items()}
    seq_strings = []
    max_len_name = 0
    max_len_seq = 0
    for j in range(len(alignment[0])):
        seq_l = []
        for i in range(len(alignment)):
            seq_l.append(alignment[i][j])
            seq_name = ordered[j]
            
        # update max len of name and seq for better formatting
        if len(seq_name) > max_len_name:
            max_len_name = len(seq_name)
        if len(seq_l) > max_len_seq:
            max_len_seq = len(seq_name)

        # width based formatting
        seq_strings.append('{:>{max_name}}: {:>{max_seq}}'.format(seq_name, ''.join(seq_l), max_name=max_len_name,
                                                                  max_seq=max_len_seq))

    return '\n'.join(seq_strings)


@timer
def align(alignment_type, input_seqs, subs_matrix, gap_penalty, sequence_type):
    '''
    alignment work
    :param alignment type, original sequence list, substitution matrix, penalty for a gap
    :return alignment result
    '''

    if sequence_type == 'protein':
        k_value = 2
    elif sequence_type == 'dna' or sequence_type == 'rna':
        k_value = 3
    n_seqs = len(input_seqs)
    # if only 1 seq, then raise error, if two then pairwise, if more then multiple
    # sequence alignment (msa)
    if n_seqs < 2:
        print('There must be at least 2 sequences in an input file!')
        raise IOError
    elif n_seqs == 2:
        if alignment_type == 'global':
            alignment = pairwise_global_module.global_alignment(input_seqs, subs_matrix, gap_penalty, sequence_type)
        else:
            alignment = pairwise_local_module.local_alignment(input_seqs, subs_matrix, gap_penalty, sequence_type)
        return alignment
    else:
        tree, seq_names, postorder_seq = guide_tree_module.generate_guide_tree(input_seqs, k=k_value)
        ordered = []
        for leaf in postorder_seq:
            if leaf[1] == 'Leaf':
                ordered.append(leaf[0])
        ordered.reverse()
        if alignment_type == 'global':
            alignment = msa_global_module.progressive_msa(
                input_seqs, msa_global_module.global_msa_sub, tree, subs_matrix, seq_names, gap_penalty, sequence_type)
        else:
            alignment = msa_local_module.progressive_msa(
                input_seqs, msa_local_module.local_msa_sub, tree, subs_matrix, seq_names, gap_penalty, sequence_type)
        return pretty_print_msa(alignment, input_seqs, alignment_type, ordered)



@click.command()
@click.argument('file')
@click.option('--alignment_type', default = 'global',
              help='Optional argument for alignment type: global or local',
              show_default=True)
@click.option('--substitution_matrix', default = 'blosum50',
              help='Optional argument for substitution matrix,'
                   'available matrixes are:\n {0}'.format('\n'.join(substitution_matrix_module.available_mat)),
              show_default=True)
@click.option('--gap_penalty', default = -8,
              help='Optional argument for gap penalty',
              show_default=True)
@click.option('--top', default = 'all',
              help='Use only top sequences from file',
              show_default=True)
@click.option('--sequence_type', default = 'protein',
              help='Optional argument for sequence type (alternative - dna)',
              show_default=True)
@click.option('--extract_fw', is_flag=True,
              help='Extract antibody\'s framework region (works only for antibodies)')
@click.option('--subregion_size', default = 'all' ,
              help='Optional argument for subregion size of sequences',
              show_default=True)
@click.option('--convert_to_dna', is_flag=True,
              help='Convert protein sequences to DNA based on codon usage fractions',
              show_default=True)
@click.option('--convert_to_protein', is_flag=True,
              help='Convert DNA sequences to protein sequences')
def main(file, alignment_type, substitution_matrix, gap_penalty, top, sequence_type,  subregion_size,
         extract_fw, convert_to_dna, convert_to_protein):
    """ FILE: Input file that contains sequences in fasta format """
    
    # change gap penalty for dna tyep
    if sequence_type == 'dna':
        gap_penalty = -1
    
    f = Figlet(font='slant')
    
    click.echo(f.renderText('SeqAlign'))
    click.echo('Using file {0} as an input'.format(file))
    click.echo('Selected alignment type is {0}'.format(alignment_type))
    click.echo('Selected substitution matrix is {0}\n'.format(substitution_matrix))
    click.echo('Selected gap penalty is {0}\n'.format(str(gap_penalty)))
    click.echo('Sub-region size is {0}'.format(subregion_size))
    click.echo('Antibody framework extraction is set to {0}'.format(extract_fw))
    click.echo('Convert to DNA: {0}'.format(convert_to_dna))
    click.echo('Convert to protein {0}'.format(convert_to_protein))
    click.echo('Using {0} sequences'.format(top))
    print()

    # Meaningful gap penalty is only negative since
    # the higher score the more similar the sequences are
    if gap_penalty > 0:
        raise ValueError('Gap penalty needs to be <= 0!')

    # retrive subsitution matrix as pandas dataframe
    subs_matrix = substitution_matrix_module.select_substitution_mat(substitution_matrix)

    # read input file
    if top != 'all':
        try:
            top = int(top)
        except:
            raise ValueError('top needs to be positive integer (at least 2)!')
        if top < 2:
            raise ValueError('top needs to be positive integer!')

    input_seqs = read_file(file, top)

    # print input sequences in a nice way
    print_seqs('Input sequences:', input_seqs)

    # in the protein sequences there might me B or Z
    # standing for asparagine/aspartic acid - asx - B - N/D
    # or glutamine/glutamic acid - glx - Z - Q/E
    # since ANARCI and substitution matrices won't work with those then
    # they must be changed. I used their frequences to calculate
    # suitable subsitution probabilities

    if sequence_type == 'protein':
        input_seqs = {k: replace_BZ(v) for k, v in input_seqs.items()}

    # extract antibody framework
    if extract_fw:
        try:
            input_seqs_fw = dict()
            for k, v in input_seqs.items():
                fw = antibody_fr_extraction.run_anarci_tool(v)
                if fw != '':
                    input_seqs_fw[k] = fw
                else:
                    print('Removed sequence {0} due to ANARCI incompatibility!'.format(k))
            # if didn't work then don't do next part!
            test_seq_len = [1 if el != '' else 0 for el in input_seqs_fw.values()]
            if sum(test_seq_len) > 2:
                input_seqs = input_seqs_fw
                del(input_seqs_fw)
            else:
                print('Too many sequences were incompatible, using original sequences instead!')
        except Exception as e:
            print('Couldn\'t extract antibody\'s framework!: {0}'.format(e))

    # sub-region size
    if subregion_size != 'all':
        try:
            subregion_size = int(subregion_size)
            if subregion_size < 1:
                raise ValueError('Subreggion size needs to be a positive integer!')
        except:
            raise ValueError('Subreggion size needs to be an integer!')
        try:
            input_seqs = {k: v[:subregion_size + 1] for k, v in input_seqs.items()}
        except Exception as e:
            print('Couldn\'t subset sequences: {0}'.format(e))

    # check whether it should convert input to DNA or not
    if sequence_type == 'protein' and convert_to_dna:
        try:
            input_seqs = {k: seq_conversion.convert_to_dna(v) for k, v in input_seqs.items()}
            sequence_type = 'dna'
        except Exception as e:
            print('Couldn\'t convert to DNA: {0}'.format(e))

    # check if needs to be converted into protein
    if sequence_type != 'protein' and convert_to_protein:
        try:
            input_seqs = {k: seq_conversion.convert_to_prot(v) for k, v in input_seqs.items()}
            sequence_type = 'protein'
        except Exception as e:
            print('Couldn\'t convert to protein: {0}'.format(e))

    # print sequences that go into alignment
    print()
    print_seqs('Input sequences (processed) for alignment:', input_seqs)

    # records in experiments for running time and gap number
    times_experiments = []
    gaps_experiments = []

    loop_times = 1
    times = []
    while loop_times:
        loop_times -= 1
        t, result = align(alignment_type, input_seqs, subs_matrix, gap_penalty, sequence_type)
        times.append(t)
    if len(result) == 3:
        print('Alignment:')
        print(result[1], end='\n\n')
        print('Sequence alignment score:', result[2], end='\n\n')
        print(result[0], end='\n\n')
        print('Average running time: %f' %(np.mean(times)))
        print('Average gap number: %d' %(result[0].count('_')/len(input_seqs)), end='\n\n')
    else:
        print('Alignment:')
        print(result, end='\n\n')
        gaps_experiments.append(int(result.count('_')/len(input_seqs)))
        print('Average running time: %f' %(np.mean(times)))
        print('Average gap number: %d' %(result.count('_')/len(input_seqs)), end='\n\n')
    times_experiments.append(np.mean(times))

    # plot output analysis
    plotting(times_experiments, '%s type alignment performance analysis' %(alignment_type))

if __name__ == '__main__':
    main()
