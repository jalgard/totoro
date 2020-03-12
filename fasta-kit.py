#!/usr/bin/env python3
'''
    written by J@1gard, 2019 (c)

    jalgard@gmail.com

    Operate with .fasta files:
    - pick sequences
    - remove sequences
    - rename entries
    - upper\lowercase convert

    ...

'''


import sys, random, os, argparse

from collections import defaultdict

'''
    Summon our options parser
'''

def GetOptParser():

    optionParser = argparse.ArgumentParser()

    optionParser.add_argument('--i', '--in',
        action='store',
        help="Input .fasta file name (or stdin if not set)")

    optionParser.add_argument('--o', '--out',
        action='store',
        help="Output .fasta file name (or stdout if not set)")

    optionParser.add_argument('--l', '--list',
        action='store',
        help="File with entry names to process [if applicable]")

    optionParser.add_argument('--a', '--action',
        action='append',
        help="Action to perform on entries from .fasta\
        (for some actions user can provide a list file with entry names)\n \
        Actions implemented:\t \
        * 'upper' make sequences uppercase, * \n\
        * 'lower' make sequences lowercase, * \n\
        * 'remove' remove listed sequences [require --list option], * \n\
        * 'keep' keep listed sequences [require --list option], *\n\
        * 'pick' pick sequences specified be index range, * \n\
        * 'lrename' rename using list, * \n\
        * 'srename' modify sequence title with splitter function, * \n\
        * 'prename' rename using prefix and increment, *\n\
        * 'crotate' circular rotate sequences so to start with leftmost 'pattern' occurance, *\n\
        * 'falength' print length of entries in stdout in format 'enrty name   entry length', * \n\
        * 'startswith' keep entries that start with 'pattern', * \n\
        * 'haspat' keep entries that contain 'pattern', * \n\
        * 'revcom' reverse-and-complement sequences. * \n\
        * 'pickname' pick sequence with 'pattern' in name. * \n\
        * 'cutbefore' cuts sequence before 'pattern' occurance. * \n\
        * 'croptail' keeps only 'range_begin' nucleotides from the start of each sequence. *\n\
        * 'ccrotate' circular rotate sequences so to start with 'range_begin' coordinate *\n\
        * 'listrename' rename fasta entries using list of new names [require --list option] *\n\
        * 'listsort' sorts fasta entries in order provided in list [require --list option] *\n\
        ")

    optionParser.add_argument('--range_begin',
        type=int, default=0,
        help="For 'pick' action, set first index, for ccrotate set the begin coordinate")
    optionParser.add_argument('--range_end',
        type=int, default=0,
        help="For 'pick' action, set last index")

    optionParser.add_argument('--splitter',
        type=str, default=' ',
        help="For 'srename' action, set splitter character")
    optionParser.add_argument('--token',
        type=int, default=0,
        help="For 'srename' action, set token index")

    optionParser.add_argument('--prefix',
        type=str, default='seq',
        help="For 'prename' action, set prefix")

    optionParser.add_argument('--pattern',
        type=str, default='',
        help="For 'crotate' | 'startswith' | 'haspat' | 'pickname' | 'cutbefore' action, set pattern")

    return optionParser

'''
    Utility loaders
'''

def ReadFasta(input_lines):

    fasta_data = []

    for line in input_lines:
        line = line.rstrip()
        if len(line) < 1:
            continue
        if line[0] == '>':
            fasta_data.append([line[1:],''])
        else:
            fasta_data[-1][1] += line

    for i, entry in enumerate(fasta_data):
        fasta_data[i].append(len(entry[1]))

    return fasta_data

def LoadList(input_list_file):

    selected_entries = defaultdict(lambda: '')

    with open(input_list_file, 'r') as list:
        for line in list:
            line = line.rstrip().split('\t')
            if len(line) > 1:
                selected_entries[line[0]] = line[1]
            else:
                selected_entries[line[0]] = ''

    return selected_entries

def WriteFasta(fasta_data, fasta_file):

    for entry in fasta_data:
        fasta_file.writelines('>{}\n{}\n'.format(entry[0], entry[1]))

revNuc = {
'A': 'T',
'T': 'A',
'G': 'C',
'C': 'G',
'a': 't',
't': 'a',
'g': 'c',
'c': 'g'
}
revNucDict = defaultdict(lambda: 'N')
revNucDict.update(revNuc)

def RevcomDNA(dna):
    revcom = ''
    global revNucDict

    for n in dna[::-1]:
        revcom += revNucDict[n]

    return revcom

'''
    Supported actions
'''

def Remove(input_fasta_entries, input_defdict, **options):

    filtered_fasta_entries = []

    for entry in input_fasta_entries:
        if input_defdict[entry[0]] == '':
             filtered_fasta_entries.append(entry)

    input_fasta_entries[:] = list(filtered_fasta_entries)

def Keep(input_fasta_entries, input_defdict, **options):

    filtered_fasta_entries = []

    for entry in input_fasta_entries:
        if input_defdict[entry[0]] != '':
             filtered_fasta_entries.append(entry)

    input_fasta_entries[:] = list(filtered_fasta_entries)


def RenameDict(input_fasta_entries, input_defdict, **options):

    for i, entry in enumerate(input_fasta_entries):
        if input_defdict[entry[0]] != '':
            input_fasta_entries[i][0] = input_defdict[entry[0]]

def RenamePrefix(input_fasta_entries, input_defdict, **options):

    entry_counter = 0
    prefix = options.get('prefix')
    for i, entry in enumerate(input_fasta_entries):
        entry_counter += 1
        input_fasta_entries[i][0] = prefix + str(entry_counter)

def RenameSplit(input_fasta_entries, input_defdict, **options):

    splitter  = options.get('splitter')
    token_num = int(options.get('token'))

    for i, entry in enumerate(input_fasta_entries):
        new_name = entry[0].split(splitter)[token_num]
        input_fasta_entries[i][0] = new_name

def UppercaseEntry(input_fasta_entries, input_defdict, **options):

    for i, entry in enumerate(input_fasta_entries):
        input_fasta_entries[i][1] = input_fasta_entries[i][1].upper()

def LowercaseEntry(input_fasta_entries, input_defdict, **options):

    for i, entry in enumerate(input_fasta_entries):
        input_fasta_entries[i][1] = input_fasta_entries[i][1].lower()

def PickRange(input_fasta_entries, input_defdict, **options):

    filtered_fasta_entries = []
    range_from = int(options.get('range_begin'))
    range_to   = int(options.get('range_end'))

    for i, entry in enumerate(input_fasta_entries):
        if i+1 >= range_from and i+1 <= range_to:
            filtered_fasta_entries.append(entry)

    input_fasta_entries[:] = list(filtered_fasta_entries)

def CircularRotate(input_fasta_entries, input_defdict, **options):

    pattern = options.get('pattern')
    for i, entry in enumerate(input_fasta_entries):
        pos = entry[1].lower().find(pattern.lower())
        if pos > 0:
            rotated = entry[1][pos:] +  entry[1][:pos]
            input_fasta_entries[i][1] = rotated

def CircularRotateByCoordinate(input_fasta_entries, input_defdict, **options):

    coord = int(options.get('range_begin'))
    for i, entry in enumerate(input_fasta_entries):
        if len(entry[1]) > coord and coord > 0:
            rotated = entry[1][coord:] +  entry[1][:coord]
            input_fasta_entries[i][1] = rotated

def PrintLength(input_fasta_entries, input_defdict, **options):

    for entry in input_fasta_entries:
        print('{}\t{}'.format(entry[0], entry[2]))
    exit()

def StartsWith(input_fasta_entries, input_defdict, **options):

    filtered_fasta_entries = []
    pattern = options.get('pattern')

    for i, entry in enumerate(input_fasta_entries):
        if entry[1][0:len(pattern)].lower() == pattern.lower():
            filtered_fasta_entries.append(entry)

    input_fasta_entries[:] = list(filtered_fasta_entries)

def HasPattern(input_fasta_entries, input_defdict, **options):

    filtered_fasta_entries = []
    pattern = options.get('pattern')

    for i, entry in enumerate(input_fasta_entries):
        if pattern.lower() in entry[1].lower():
            filtered_fasta_entries.append(entry)

    input_fasta_entries[:] = list(filtered_fasta_entries)

def RevcomEntries(input_fasta_entries, input_defdict, **options):

    for i, entry in enumerate(input_fasta_entries):
        input_fasta_entries[i][1] = RevcomDNA(entry[1])

def Pickname(input_fasta_entries, input_defdict, **options):

    filtered_fasta_entries = []
    pattern = options.get('pattern')

    for entry in input_fasta_entries:
        if pattern in entry[0]:
             filtered_fasta_entries.append(entry)

    input_fasta_entries[:] = list(filtered_fasta_entries)

def Cutbefore(input_fasta_entries, input_defdict, **options):

    filtered_fasta_entries = []
    pattern = options.get('pattern')


    for entry in input_fasta_entries:
        pcount = int(options.get('range_begin'))
        if entry[1].count(pattern) >= pcount:
            start = entry[1].find(pattern)
            while start >= 0 and pcount > 1:
                start = entry[1].find(pattern, start+len(pattern))
                pcount -= 1
            left_n = entry[0] + "_Left"
            right_n = entry[0] + "_Right"
            left = entry[1][:start]
            right = entry[1][start:]
            filtered_fasta_entries.append([left_n, left])
            filtered_fasta_entries.append([right_n, right])
        else:
            filtered_fasta_entries.append(entry)

    input_fasta_entries[:] = list(filtered_fasta_entries)

def Croptail(input_fasta_entries, input_defdict, **options):

    plen = int(options.get('range_begin'))
    for i, entry in enumerate(input_fasta_entries):
        input_fasta_entries[i][1] = input_fasta_entries[i][1][:plen]

def RenameList(input_fasta_entries, input_defdict, **options):

    for i, entry in enumerate(input_fasta_entries):
        if i < len(input_defdict):
            input_fasta_entries[i][0] = input_defdict[i]

def SortList(input_fasta_entries, input_defdict, **options):
    order = []
    with open(options.get('list'), 'r') as listfile:
        for line in listfile:
            order.append(line.rstrip())
    sorted_fasta_entries = []
    for name in order:
        for entry in input_fasta_entries:
            if entry[0] == '>' + name:
                sorted_fasta_entries.append(entry)
    input_fasta_entries[:] = list(sorted_fasta_entries)



ACTIONS = {
'upper' :    UppercaseEntry,
'lower' :    LowercaseEntry,
'remove' :   Remove,
'keep' :     Keep,
'lrename' :  RenameDict,
'srename' :  RenameSplit,
'prename' :  RenamePrefix,
'pick' :     PickRange,
'crotate' :  CircularRotate,
'falength':  PrintLength,
'haspat' :   HasPattern,
'startswith':StartsWith,
'revcom' :   RevcomEntries,
'pickname' : Pickname,
'cutbefore' : Cutbefore,
'croptail' : Croptail,
'ccrotate' : CircularRotateByCoordinate,
'listrename' : RenameList,
'listsort' : SortList
}

if __name__ == '__main__':

    runArgs = GetOptParser().parse_args(sys.argv[1:])

    inputSource = sys.stdin
    if runArgs.i is not None:
        if os.path.isfile(runArgs.i):
            inputSource = open(runArgs.i, 'r')

    outputDest = sys.stdout
    if runArgs.o is not None:
        try:
            outputDest = open(runArgs.o, 'w')
        except IOError:
            sys.stderr.writelines('Error: unable to open output file {},\
            writing to stdout instead\n'.format(runArgs.o))
            outputDest = sys.stdout

    # load fasta entries
    fasta_data = ReadFasta(inputSource)
    options = {
        'splitter' : runArgs.splitter,
        'token' : runArgs.token,
        'prefix' : runArgs.prefix,
        'range_begin': runArgs.range_begin,
        'range_end': runArgs.range_end,
        'pattern' : runArgs.pattern,
        'list' : runArgs.l
    }

    # here is a good place to sort entries
    # if necessary


    selectedEntries = defaultdict(lambda: 0)
    [selectedEntries.update({name[0]:1}) for name in fasta_data]

    if runArgs.l is not None:
        try:
            selected_list = LoadList(runArgs.l)
            selectedEntries = selected_list
        except:
            sys.stderr.writelines('Error: unable read list file {},\
            acting on all entries in .fasta\n'.format(runArgs.l))

    actionFunc = ACTIONS[runArgs.a[0]]

    actionFunc(fasta_data, selectedEntries, **options)

    WriteFasta(fasta_data, outputDest)

    #outputDest.close()
