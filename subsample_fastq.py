import sys, random, os, argparse


from collections import defaultdict

def GetOptParser():

    optionParser = argparse.ArgumentParser()

    optionParser.add_argument('--i', '--fastq1',
        action='store',
        help="Input FASTQ with SE reads or R1 of PE reads")

    optionParser.add_argument('--r', '--fastq2',
        action='store',
        help="Input FASTQ with R2 only for PE reads")

    optionParser.add_argument('--o', '--out1',
        action='store',
        help="Output FASTQ for SE reads or R1 of PE reads")

    optionParser.add_argument('--u', '--out2',
        action='store',
        help="Output FASTQ for R2 of PE reads")

    optionParser.add_argument('--s', '--sample',
        action='store',
        help="Subsample this fraction of reads, % [e.g. 40]")

    return optionParser



def SubsampleReadsPairedEnd(fraq, r1, r2, o1, o2):
    rf1 = open(r1, 'r').readlines()
    rf2 = open(r2, 'r').readlines()
    of1 = open(o1, 'w')
    of2 = open(o2, 'w')

    i = 0
    while i < len(rf1):
        s = random.randint(1, 100)
        if s <= fraq:
            of1.writelines(rf1[i])
            of1.writelines(rf1[i+1])
            of1.writelines(rf1[i+2])
            of1.writelines(rf1[i+3])

            of2.writelines(rf2[i])
            of2.writelines(rf2[i+1])
            of2.writelines(rf2[i+2])
            of2.writelines(rf2[i+3])
        i += 4

    of1.close()
    of2.close()


def SubsampleReadsSingleEnd(fraq, r1, o1):
    rf1 = open(r1, 'r').readlines()
    of1 = open(o1, 'w')

    i = 0
    while i < len(rf1):
        s = random.randint(1, 100)
        if s <= fraq:
            of1.writelines(rf1[i])
            of1.writelines(rf1[i+1])
            of1.writelines(rf1[i+2])
            of1.writelines(rf1[i+3])

        i += 4

    of1.close()



if __name__ == '__main__':

    runArgs = GetOptParser().parse_args(sys.argv[1:])

    subsample_fraction = 100

    if runArgs.s is not None:
        subsample_fraction = int(runArgs.s)

    fastq_in1  = ''
    fastq_in2  = ''


    if runArgs.i is not None and os.path.isfile(runArgs.i):
        fastq_in1 = runArgs.i
    else:
        sys.stderr.writelines('Input file not set or does not exist!')
        exit()

    if runArgs.r is not None and os.path.isfile(runArgs.r):
        fastq_in2 = runArgs.r
        print('Mode is set to paired-end')

    fastq_out1 = fastq_in1 + '.subsample'
    fastq_out2 = fastq_in2 + '.subsample'

    if runArgs.o is not None:
        fastq_out1 = runArgs.o
    if runArgs.u is not None:
        fastq_out2 = runArgs.u

    if fastq_in1 != '' and fastq_in2 != '':
        SubsampleReadsPairedEnd(subsample_fraction, fastq_in1, fastq_in2, fastq_out1, fastq_out2)
    else:
        SubsampleReadsSingleEnd(subsample_fraction, fastq_in1, fastq_out1)
