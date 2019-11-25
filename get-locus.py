#!/usr/bin/env python3
import sys, os, argparse


def GetOptParser():

    optionParser = argparse.ArgumentParser()

    optionParser.add_argument('--i', '--in',
        action='store',
        help="Input .fasta file name (or stdin if not set)")

    optionParser.add_argument('--o', '--out',
        action='store',
        help="Output .fasta file name (or stdout if not set)")

    optionParser.add_argument('--l', '--locus',
        action='store',
        help="Locus description in format Chr1,start,end or just Chr1.\nIn 1-based coordinate system")

    optionParser.add_argument('--s', '--space',
        action='store_true',
        help="Store only prefix of fasta entry name before first space character")

    return optionParser


def ReadFasta(input_lines, split_space = False):

    fasta_data = []

    for line in input_lines:
        line = line.rstrip()
        if len(line) < 1:
            continue
        if line[0] == '>':
            if split_space == False:
                fasta_data.append([line[1:],''])
            else:
                fasta_data.append([line[1:].split(' ')[0],''])
        else:
            fasta_data[-1][1] += line

    for i, entry in enumerate(fasta_data):
        fasta_data[i].append(len(entry[1]))

    return fasta_data

if __name__ == '__main__':

    runArgs = GetOptParser().parse_args(sys.argv[1:])
    inputSource = sys.stdin
    if runArgs.i is not None:
        if os.path.isfile(runArgs.i):
            inputSource = open(runArgs.i, 'r')
        else:
            sys.stderr.write('No such file {}.\nExiting with an error!'.format(unArgs.i))
            exit()

    outputDest = sys.stdout
    if runArgs.o is not None:
        try:
            outputDest = open(runArgs.o, 'w')
        except IOError:
            sys.stderr.writelines('Error: unable to open output file {},\
            writing to stdout instead\n'.format(runArgs.o))
            outputDest = sys.stdout

    entry_name = None
    cut_from   = 0
    cut_to     = 0

    if runArgs.l is not None:
        locus_details = runArgs.l.split(',')
        if len(locus_details) == 1:
            entry_name = locus_details[0]
        elif len(locus_details) == 3:
            entry_name = locus_details[0]
            cut_from   = int(locus_details[1])
            cut_to     = int(locus_details[2])
        else:
            sys.stderr.writelines('Please specify locus in proper format: \'chrname,start,end\' or \'chrname\'\nExiting with an error!')
            exit()

    else:
        sys.stderr.write('No locus specified!\nExiting with an error!')
        exit()

    split_space = False
    if runArgs.s is not None:
        split_space = True
    fasta_entries = ReadFasta(inputSource, split_space)

    selected_locus_seq = None
    for entry in fasta_entries:
        if entry[0] == entry_name:
            selected_locus_seq = entry[1]

    if selected_locus_seq is None:
        sys.stderr.write('No scaffold with name {}!'.format(entry_name))
        exit()

    if cut_from == 0 and cut_to == 0:
        outputDest.writelines('>{}\n{}\n'.format(entry_name, selected_locus_seq))
    elif cut_to >= cut_from and cut_to < len(selected_locus_seq):
        outputDest.writelines('>{}\n{}\n'.format(entry_name, selected_locus_seq[cut_from-1: cut_to]))
    elif cut_to >= len(selected_locus_seq):
        sys.stderr.write('End coordinate {} is out of scaffold\'s length {}!'.format(cut_to, len(selected_locus_seq)))
        exit()
    else:
        sys.stderr.write('Can get such locus.\nPlease check --locus key again.')
        exit()
