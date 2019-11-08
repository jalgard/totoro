#!/usr/bin/env python3

'''
    Written by J@lgard, 2019

    Plot histogram of incoming list of values

    Current version uses Seaborn package
'''

import sys, os, argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

def GetOptParser():

    optionParser = argparse.ArgumentParser()

    optionParser.add_argument('--i', '--in',
        action='store',
        help="Input .fasta file name (or stdin if not set)")

    optionParser.add_argument('--o', '--out',
        action='store',
        help="Output .fasta file name (or stdout if not set)")

    return optionParser


if __name__ == '__main__':

    runArgs = GetOptParser().parse_args(sys.argv[1:])
    inputSource = sys.stdin
    if runArgs.i is not None:
        if os.path.isfile(runArgs.i):
            inputSource = open(runArgs.i, 'r')
        else:
            sys.stderr.write('No such file {}.\nExiting with an error!'.format(unArgs.i))
            exit()

    outputDest = None
    if runArgs.o is not None:
        outputDest = runArgs.o
    else:
        sys.stderr.write('No output file specified.\nExiting with an error!')
        exit()


    data_values = []
    for line in inputSource:
        data_values.append(float(line.rstrip()))

    q25 = np.quantile(data_values, 0.25)
    q75 = np.quantile(data_values, 0.75)
    iqr = q75 - q25

    filtered_data = [x for x in data_values if (x > q25 - 1.5*iqr) and (x < q75 + 1.5*iqr)]

    f, (ax_box, ax_hist) = plt.subplots(2, sharex=True,
                                    gridspec_kw={"height_ratios": (.05, .95)})

    sns.boxplot(filtered_data, ax=ax_box)
    sns.distplot(filtered_data, ax=ax_hist)
    ax_box.set(yticks=[])
    sns.despine(ax=ax_hist)
    sns.despine(ax=ax_box, left=True)
    plt.savefig(outputDest, dpi=200)
