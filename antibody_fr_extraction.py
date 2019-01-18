
import subprocess

# given script uses anarci tool (can be downloaded from http://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/ANARCI.php)
# to extract antibody framework sequences from entire sequence by using IMGT nomenclature

# from IMGT (http://www.imgt.org/IMGTScientificChart/Nomenclature/IMGT-FRCDRdefinition.html)
# we can see that framework numbers are following:
# FR1 - 1...26
# FR2 - 39...55
# FR3 - 66...104

#ANARCI -i EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSA --scheme chothia


def process_output(output):

    framework_region = []

    # framework positions
    FR_range = []
    FR_range.extend(list(range(1, 27)))
    FR_range.extend(list(range(39, 56)))
    FR_range.extend(list(range(66, 105)))

    # if position in framework region then include this aa to final sequence
    for line in output.split('\n'):
        if line.startswith('H '):
            l = line.split()
            pos = int(l[1])
            aa = l[2]
            if pos in FR_range and aa != '-':
                framework_region.append(aa)

    return ''.join(framework_region)


def run_anarci_tool(input_seq, scheme = 'imgt'):
    try:
        completed = subprocess.run(
            ['ANARCI', '-i', input_seq, '--scheme', scheme],
            check=True, stdout=subprocess.PIPE
            )
    except subprocess.CalledProcessError as error:
        raise Exception('ERROR: {0}'.format(error))
    output = completed.stdout.decode('utf-8')
    framework_region = process_output(output)
    return framework_region

