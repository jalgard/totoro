# This script plots codon usage as radar plot

# Libraries
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from math import pi
import sys, os, argparse, random
from collections import defaultdict
import os.path
import seaborn as sns

# for PCA
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA



def GetOptParser():

    optionParser = argparse.ArgumentParser()

    optionParser.add_argument('--sl', '--species_list',
        action='store',
        help="Species list, tab-separated file linking species name to CDS fasta file")

    optionParser.add_argument('--o', '--out',
        action='store',
        help="Output in png format")

    return optionParser


if __name__ == '__main__':

    runArgs = GetOptParser().parse_args(sys.argv[1:])
    speciesListFile = None
    if runArgs.sl is not None:
        if os.path.isfile(runArgs.sl):
            speciesListFile = runArgs.sl
        else:
            sys.stderr.write('Wrong species list file or species list file not set!')
            exit()

    outputPrefix = "CUplot"
    if runArgs.o is not None:
        outputPrefix = runArgs.o
    # read species file
    # tab-separated line "species \t cds fasta file"

    codonsLeu = ['CTA','CTT','CTG','CTC','TTA','TTG']
    codonsArg = ['CGA','CGT','CGG','CGC','AGA','AGG']
    codonsSer = ['TCA','TCT','TCG','TCC','AGC','AGT']

    dataset = {}
    with open(speciesListFile, 'r') as species_file:
        for line in species_file:
            toks = line.rstrip().split('\t')
            if os.path.isfile(toks[1]):
                dataset[toks[0]] = toks[1]
            else:
                print('Can\'t find file {}, skipping'.format(toks[1]))

    # initialize codons
    codons = []
    for i in 'AGTC':
        for j in 'AGTC':
            for k in 'AGTC':
                codons.append(i+j+k)

    cuTable = pd.DataFrame(columns=['species'] + codons)
    # count codon usage for each species

    palette = sns.color_palette(None, len(dataset.keys()))
    palette5 = sns.color_palette(None, 5)

    for species in dataset:
        print('Counting codon usage table for {}'.format(species))
        cuSpecies = {'species': species }
        cu = defaultdict(lambda : 0.0)
        ct = 0.0
        cdsEntries = []
        ignoredCds = 0
        with open(dataset[species], 'r') as inputCdsFile:
            for line in inputCdsFile:
                if line[0] == '>':
                    cdsEntries.append('')
                else:
                    cdsEntries[-1] += line.rstrip().upper()
        for cds in cdsEntries:
            if len(cds) % 3 == 0:
                for p in range(0, len(cds), 3):
                    ct += 1
                    cu[cds[p:p+3]] += 1
            else:
                ignoredCds += 1
        for codon in codons:
            cuSpecies[codon] = cu[codon] / ct
        cuTable = cuTable.append(cuSpecies, ignore_index=True)
        print('Ignored CDS: {}'.format(ignoredCds))


    cuTable.set_index('species', inplace=True)
    cuTable.sort_values(by=[dataset.keys()[0]],axis=1, ascending=True, inplace=True, kind='quicksort', na_position='last')
    cuTable.reset_index(inplace=True)
    print(cuTable)

    categories = list(cuTable)[1:]
    N = len(categories)

    angles = [n / float(N) * 2 * pi for n in range(N)]
    angles += angles[:1]

    my_dpi=300
    fig=plt.figure(figsize=(4000/my_dpi, 4000/my_dpi), dpi=my_dpi)
    ax = plt.subplot(111, polar=True)
    plt.rcParams["font.family"] = "Times New Roman"
    ax.set_theta_offset(pi / 2)
    ax.set_theta_direction(-1)

    plt.xticks(angles[:-1], categories)
    for i, label in enumerate(ax.get_xticklabels()):
        label.set_rotation(i*45)
        label.set_fontsize(12)
        #label.set_fontproperties('ticks_font')
        label.set_horizontalalignment("center")
        label.set_rotation_mode("anchor")


    ax.set_rlabel_position(0)
    plt.yticks([1.0/64, 2.0/64], ["1/64", "2/64"], color="black", size=10)
    plt.ylim(0,0.055)

    for i in range(len(dataset.keys())):
        values = cuTable.loc[i].drop('species').values.flatten().tolist()
        values += values[:1]
        #if i == 0:
        #    ax.plot(angles, values, linewidth=4, linestyle='solid', color='#87CEFA', label = r"$\it{C. bombi}$")
        #    ax.fill(angles, values, '#87CEFA', alpha=0.1)
        #else:
        ax.plot(angles, values, linewidth=1, linestyle='solid', color=palette[i], label=dataset.keys()[i])

    plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))
    plt.tight_layout()
    plt.savefig('{}_Radar.png'.format(outputPrefix), format='png', dpi=300)
    plt.savefig('{}_Radar.pdf'.format(outputPrefix), format='pdf', dpi=300)

    codonData = cuTable.loc[:, codons].values
    codonData = StandardScaler().fit_transform(codonData)
    pca = PCA(n_components=2)
    pComponents = pca.fit_transform(codonData)
    pDf = pd.DataFrame(data = pComponents, columns = ['PC 1', 'PC 2'])
    finalDf = pd.concat([pDf, cuTable[['species']]], axis = 1)
    print(finalDf)
    plt.clf()
    plt.cla()
    fig = plt.figure(figsize = (12,12))
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel('Principal Component 1', fontsize = 15)
    ax.set_ylabel('Principal Component 2', fontsize = 15)
    ax.set_title('2 component PCA', fontsize = 20)

    for target, color in zip(dataset.keys(), palette):
        indicesToKeep = finalDf['species'] == target
        ax.scatter(finalDf.loc[indicesToKeep, 'PC 1'], finalDf.loc[indicesToKeep, 'PC 2'], c = color, s = 90)
    ax.legend(dataset.keys())
    ax.grid()
    plt.tight_layout()
    plt.savefig('{}_PCA.png'.format(outputPrefix), format='png', dpi=300)
    plt.savefig('{}_PCA.pdf'.format(outputPrefix), format='pdf', dpi=300)

    #  plot codon family radars
    for codonsFamily, familyName in zip([codonsSer, codonsLeu, codonsArg], ['Ser', 'Leu', 'Arg']):
        plt.clf()
        plt.cla()

        categories = codonsFamily
        N = len(categories)

        angles = [n / float(N) * 2 * pi for n in range(N)]
        angles += angles[:1]

        my_dpi=300
        fig=plt.figure(figsize=(4000/my_dpi, 4000/my_dpi), dpi=my_dpi)
        ax = plt.subplot(111, polar=True)
        plt.rcParams["font.family"] = "Times New Roman"
        ax.set_theta_offset(pi / 2)
        ax.set_theta_direction(-1)

        plt.xticks(angles[:-1], categories)
        for i, label in enumerate(ax.get_xticklabels()):
            label.set_rotation(i*45)
            label.set_fontsize(22)
            label.set_horizontalalignment("center")
            label.set_rotation_mode("anchor")


        ax.set_rlabel_position(0)
        plt.yticks([1.0/64, 2.0/64], ["1/64", "2/64"], color="black", size=20)
        plt.ylim(0,0.055)

        #for i in range(len(dataset.keys())):
        for j, i in enumerate(['B.saltans','T.brucei','L.tarentolae','E.monterogeii']):
            values = cuTable.loc[cuTable['species']==i, codonsFamily].values.flatten().tolist()
            values += values[:1]
            #ax.plot(angles, values, linewidth=4, linestyle='solid', color=palette[i], alpha=0.5, label=dataset.keys()[i])
            ax.plot(angles, values, linewidth=4, linestyle='solid', color=palette5[j], alpha=0.5, label=i)

        plt.legend(loc='upper right', fontsize=32, bbox_to_anchor=(0.1, 0.1))
        plt.tight_layout()
        plt.savefig('{}_{}_Radar.png'.format(outputPrefix, familyName), format='png', dpi=300)
        plt.savefig('{}_{}_Radar.pdf'.format(outputPrefix, familyName), format='pdf', dpi=300)
