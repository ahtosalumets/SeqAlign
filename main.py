
import click
from pyfiglet import Figlet
import substitution_matrix_module
import pairwise_global_module
import pairwise_local_module
import msa_global_module
import guide_tree_module



def read_file(file):
    seqs = dict()
    with open(file, 'r') as f:
        content = f.readlines()
    for i in range(0, len(content)-1, 2):
        seq_name = content[i].strip().replace('>', '')
        seqs[seq_name] = content[i+1].strip()
    return seqs

def pretty_print_msa(alignment, sequences):
    sequences_inv = {v : k for k, v in sequences.items()}
    seq_strings = []
    for j in range(len(alignment[0])):
        seq_l = []
        for i in range(len(alignment)):
            seq_l.append(alignment[i][j])
        try:
            seq_name = sequences_inv[''.join(seq_l).replace('_', '')]
        # during testing this current version of msa did chop some parts of few seqences
        # and I couldn't find the bug causing it
        except:
            seq_name = 'wrong'
        seq_strings.append(seq_name + ': ' + ''.join(seq_l))
    return '\n'.join(seq_strings)


@click.command()

@click.argument('file', 'alignment_type')

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

def main(file, alignment_type, substitution_matrix, gap_penalty):
    """ FILE: Input file that contains sequences in fasta format """
    f = Figlet(font='slant')
    click.echo(f.renderText('SeqAlign'))
    click.echo('Using file {0} as an input'.format(file))
    click.echo('Selected alignment type is {0}'.format(alignment_type))
    click.echo('Selected subsitution matrix is {0}\n'.format(substitution_matrix))
    click.echo('Selected gap peanlty is {0}\n'.format(str(gap_penalty)))

    # Meaningful gap penalty is only negative since
    # the higher score the more similar the sequences are
    if gap_penalty > 0:
        raise ValueError('Gap penalty needs to be <= 0!')

    # retrive subsitution matrix as pandas dataframe
    subs_matrix = substitution_matrix_module.select_substitution_mat(substitution_matrix)

    # read input file
    input_seqs = read_file(file)

    # if only 1 seq, then raise error, if two then pairwise, if more then multiple
    # sequence alignment (msa)
    n_seqs = len(input_seqs)
    if n_seqs < 2:
        print('There must be at least 2 sequences in an input file!')
        raise IOError
    elif n_seqs == 2:
        if alignment_type == 'global':
            alignment = pairwise_global_module.global_alignment(input_seqs, subs_matrix, gap_penalty)
        else:
            alignment = pairwise_local_module.local_alignment(input_seqs, subs_matrix, gap_penalty)
        print(alignment)
        print('\n')
    else:
        tree, seq_names, postorder_seq = guide_tree_module.generate_guide_tree(input_seqs, k=3)
        alignment = msa_global_module.progressive_msa(input_seqs, msa_global_module.global_msa_sub, tree, subs_matrix, seq_names, gap_penalty)
        print(pretty_print_msa(alignment, input_seqs))

if __name__=='__main__':
    main()