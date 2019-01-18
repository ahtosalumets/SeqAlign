import click
import time
import numpy as np
from functools import wraps
from pyfiglet import Figlet
import substitution_matrix_module
import pairwise_global_module
import pairwise_local_module
import msa_global_module
import msa_local_module
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import guide_tree_module


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


def read_file(file, max_len=100):
    '''
    fasta file read function
    :param filename
    :return file
    '''

    seqs = dict()
    with open(file, 'r') as f:
        content = f.read().split('>')[1:]
    count = 0
    for sequence in content:
        count += 1
        if count > max_len:
            break
        sequence = sequence.split('\n')
        seq_tmp = ''
        for seq in sequence[1:]:
            seq_tmp += seq
        seqs[sequence[0].split(' ')[0]] = seq_tmp

    return seqs


def pretty_print_msa(alignment, sequences, alignment_type, ordered):
    '''
    result output function
    :param alignment result, original sequence list, alignment type, ordered sequence names
    :return
    '''
    sequences_inv = {v : k for k, v in sequences.items()}
    seq_strings = []
    for j in range(len(alignment[0])):
        seq_l = []
        for i in range(len(alignment)):
            seq_l.append(alignment[i][j])
        try:
            if alignment_type == 'global':
                seq_name = sequences_inv[''.join(seq_l).replace('_', '')]
            else:
                seq_name = ordered[j]
        # during testing this current version of msa did chop some parts of few seqences
        # and I couldn't find the bug causing it
        except:
            seq_name = 'wrong'
        
        seq_strings.append(seq_name + ': ' + ''.join(seq_l))
    return '\n'.join(seq_strings)


@timer
def align(alignment_type, input_seqs, subs_matrix, gap_penalty, sequence_type):
    '''
    alignment work
    :param alignment type, original sequence list, substitution matrix, penalty for a gap
    :return alignment result
    '''

    if sequence_type == 'protein':
        k_value = 3
    elif sequence_type == 'dna' or sequence_type == 'rna':
        k_value = 5
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
## Where there are two required parameters
#@click.argument('file', 'alignment_type')
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
@click.option('--length', nargs=3, type=int, default = [5,5,1],
              help='Optional argument for maximum file read length',
              show_default=True)
@click.option('--sequence_type', default = 'protein',
              help='Optional argument for sequence type',
              show_default=True)
def main(file, alignment_type, substitution_matrix, gap_penalty, length, sequence_type):
    """ FILE: Input file that contains sequences in fasta format """
    
    # change gap penalty for dna tyep
    if sequence_type == 'dna' or sequence_type == 'rna':
        gap_penalty = -1
    
    f = Figlet(font='slant')
    
    click.echo(f.renderText('SeqAlign'))
    click.echo('Using file {0} as an input'.format(file))
    click.echo('Selected alignment type is {0}'.format(alignment_type))
    click.echo('Selected subsitution matrix is {0}'.format(substitution_matrix))
    click.echo('Selected gap peanlty is {0}'.format(str(gap_penalty)))
    click.echo('Selected sequence type is {0}\n'.format(str(sequence_type)))

    # Meaningful gap penalty is only negative since
    # the higher score the more similar the sequences are
    if gap_penalty > 0:
        raise ValueError('Gap penalty needs to be <= 0!')

    # retrive subsitution matrix as pandas dataframe
    subs_matrix = substitution_matrix_module.select_substitution_mat(substitution_matrix)
    
    # records in experiments for running time and gap number
    times_experiments = []
    gaps_experiments = []
    
    
    for i in range(length[0], length[1]+1, length[2]):
        # read input file
        input_seqs = read_file(file, i)
        print('File read finished, length:%d\n'%(len(input_seqs)))
        # loop time for getting average time
        loop_times = 1
        times = []
        while loop_times:
            loop_times -= 1
            t, result = align(alignment_type, input_seqs, subs_matrix, gap_penalty, sequence_type)
            times.append( t )
        if len(result) == 3:
            print(result[1], end='\n\n')
            print('Sequence alignment score:' ,result[2], end='\n\n')
            print(result[0], end='\n\n')
            print('Average running time: %f' %(np.mean(times)))
            print('Average gap number: %d' %(result[0].count('_')/len(input_seqs)), end='\n\n')
        else: 
            print(result, end='\n\n')
            gaps_experiments.append(int(result.count('_')/i))
            print('Average running time: %f' %(np.mean(times)))
            print('Average gap number: %d' %(result.count('_')/len(input_seqs)), end='\n\n')
        times_experiments.append(np.mean(times))

    # plot output analysis
    plotting(times_experiments, '%s type alignment performance analysis' %(alignment_type))

if __name__ == '__main__':
    main()
