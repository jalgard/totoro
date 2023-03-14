'''
    this simple script takes only
    that reads, not listed in input
    reads list file


    run example (results are printed in stdout):
    ./filter-mapped_fastq.py input.fastq input.list

'''

import sys
from collections import defaultdict
from itertools import islice

listed_read_counter = 0
mapped_reads_list = defaultdict(lambda : 0)

with open(sys.argv[2], 'r') as imapped:
    for line in imapped:
        mapped_reads_list[line.rstrip()] = 1
        listed_read_counter += 1


total_read_counter   = 0
printed_read_counter = 0
with open(sys.argv[1], 'r') as infile:
    while True:
        read_data = list(islice(infile, 4))
        if not read_data:
            break
        read_name = read_data[0].split(' ')[0][1:]
        if mapped_reads_list[read_name] == 0:
            sys.stdout.write('{}{}{}{}'.format(read_data[0], read_data[1], read_data[2], read_data[3]))
            printed_read_counter += 1
        total_read_counter += 1

sys.stderr.write('Total reads: {}\n'.format(total_read_counter))
sys.stderr.write('Selected reads: {}\n'.format(printed_read_counter))
sys.stderr.write('Size of list: {}\n'.format(listed_read_counter))
