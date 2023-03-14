'''
    This script is used to plot PCA for kmer frequecny vectors
    It uses .fasta file as input and plots a dot per fasta entry
    You can supply multiple fasta files, script will use separate
    dot color per fasta file.

    Author: jalgard
    Date:   March 2021
'''

# Libraries
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from math import pi
from mpl_toolkits.mplot3d import Axes3D
import sys, os, argparse, random
from collections import defaultdict
import os.path
import seaborn as sns

# for PCA
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

def GetOptParser():

    optionParser = argparse.ArgumentParser()

    optionParser.add_argument('--f', '--fasta',
        action='append',
        help="Input fasta file with sequences")

    optionParser.add_argument('--o', '--out',
        action='store',
        help="Output in png format")

    optionParser.add_argument('--k', '--kmer',
        action='store',
        help="Kmer size in range 2-12 [default: 3]")

    return optionParser

def GenerateKmerOrder(kmer):
    kmer_order_dict = {'':''}
    sorted_kmer_order = []
    for n in range(kmer):
        for k in list(kmer_order_dict.keys()):
            for l in 'AGTC':
                kmer_order_dict[k+l] = k+ l
    for i in sorted(kmer_order_dict.keys()):
        if len(i) == kmer:
            sorted_kmer_order.append(i)
    return sorted_kmer_order


def ProcessKmersPerEntry(fasta, kmerl, sorted_kmer_order):
    freqVecs = []
    freqVecLabels = []
    countItem = defaultdict(lambda : 0.0)
    totalItem = 0.0
    labelItem = ''
    with open(fasta, 'r') as data_input_file:
        for line in data_input_file:
            if line[0] == '>':
                if len(countItem) > 1:
                    item_vec = []
                    for k in sorted_kmer_order:
                        item_vec.append(countItem[k] / totalItem)
                    freqVecs.append(item_vec)
                    freqVecLabels.append(labelItem)
                countItem = defaultdict(lambda : 0.0)
                totalItem = 0.0
                labelItem = line[1:].rstrip()


            else:
                for pos in range(0, len(line.lstrip().rstrip())-kmerl):
                    countItem[line[pos:pos+kmerl]] += 1
                    totalItem += 1

    if len(countItem) > 1:
        item_vec = []
        for k in sorted_kmer_order:
            item_vec.append(countItem[k] / totalItem)
        freqVecs.append(item_vec)
        freqVecLabels.append(labelItem)

    return freqVecs, freqVecLabels

if __name__ == '__main__':

    runArgs = GetOptParser().parse_args(sys.argv[1:])

    input_fasta_files = []
    if runArgs.f is not None:
        for ffile in runArgs.f:
            input_fasta_files.append(ffile)
    else:
        sys.stderr.write('No input files... Exiting...')
        exit()


    output_file = 'kmer_pca'
    if runArgs.o is not None:
        output_file = runArgs.o


    kmer = 3
    if runArgs.k is not None:
        kmer = int(runArgs.k)

    sorted_kmer_order = GenerateKmerOrder(kmer)
    kfqTable = pd.DataFrame(columns=['entry', 'color'] + sorted_kmer_order)
    palette = sns.color_palette(None, len(input_fasta_files))

    for idx, ff in enumerate(input_fasta_files):
        freqs, labs = ProcessKmersPerEntry(ff, kmer, sorted_kmer_order)
        for entry, fvector in zip(labs, freqs):
            kfqRow = {'color' : palette[idx], 'entry' : entry}
            for k, f in zip(sorted_kmer_order, fvector):
                kfqRow[k] = f
            kfqTable = kfqTable.append(kfqRow, ignore_index=True)


    frequencyData = kfqTable.loc[:, sorted_kmer_order].values
    frequencyData = StandardScaler().fit_transform(frequencyData)
    pca = PCA(n_components=3)
    pComponents = pca.fit_transform(frequencyData)
    pDf = pd.DataFrame(data = pComponents, columns = ['PC 1', 'PC 2', 'PC 3'])
    finalDf = pd.concat([pDf, kfqTable[['color', 'entry']]], axis = 1)
    print(finalDf)
    plt.clf()
    plt.cla()
    fig = plt.figure(figsize = (12,12))
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('Principal Component 1', fontsize = 15)
    ax.set_ylabel('Principal Component 2', fontsize = 15)
    ax.set_zlabel('Principal Component 3', fontsize = 15)
    ax.set_title('3 component PCA', fontsize = 20)

    for color in palette:
        indicesToKeep = finalDf['color'] == color
        ax.scatter(finalDf.loc[indicesToKeep, 'PC 1'], finalDf.loc[indicesToKeep, 'PC 2'], finalDf.loc[indicesToKeep, 'PC 3'],  c = color, s = 90)
    #ax.legend(dataset.keys())
    ax.grid()
    plt.tight_layout()
    plt.savefig('{}_PCA.png'.format(output_file), format='png', dpi=300)
    plt.savefig('{}_PCA.pdf'.format(output_file), format='pdf', dpi=300)
