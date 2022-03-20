import argparse
import os
import re
from random import random

import numpy as np
from Bio import SeqIO
from tqdm import tqdm


def get_step(reference_length, length_mean, depth):
    num_reads = reference_length // length_mean * depth
    step_mean = (reference_length - length_mean) / (num_reads - 1)
    step_mean = round(step_mean / 100) * 100
    step_std = step_mean * 0.1
    return step_mean, step_std


def introduce_errors(subs, indel):
    def f(read):
        i = 0
        condition = subs > 0.0 or indel > 0.0
        while condition and i < len(read):
            if random() < subs:
                r = random()
                if r < 0.25:
                    c = 'A'
                elif r < 0.5:
                    c = 'C'
                elif r < 0.75:
                    c = 'G'
                else:
                    c = 'T'
                read = read[:i] + c + read[i+1:]
                i += 1
                continue

            if random() < indel:
                r = random()
                if r < indel / 2:
                    read = read[:i] + read[i+1:]
                    continue
                if r > 1 - indel / 2:
                    read = read[:i] + 2 * read[i] + read[i+1:]  # Only homopolymer insertion, not a random one
                    i += 1
                    continue
            i += 1

        return read
    return f

class RandomSequenceGenerator:
    def __init__(self, length_mean, length_std, step_mean, step_std, stop):
        self.length_mean = length_mean
        self.length_std = length_std
        self.step_mean = step_mean
        self.step_std = step_std
        self.stop = stop

    def random_sequence(self):
        start = 0
        step = 0
        while start + step  < self.stop:
            start = start + step
            length = int(np.random.normal(self.length_mean, self.length_std))
            step = int(np.random.normal(self.step_mean, self.step_std))
            end = start + length
            if end < self.stop:
                yield (start, end)
            else:
                break

    def total(self):
        return int(self.stop/self.step_mean)


def sample_strand(reference, reads_list, strand, g, introduce_errors):
    idx = len(reads_list)
    stop = len(reference)

    for i, (start, end) in tqdm(enumerate(g.random_sequence()), total=g.total()):
        read = reference[start:end]
        read.id = str(idx + i)
        read.seq = introduce_errors(read.seq)
        if strand == '+':
            read.description = f'idx={read.id}, strand=+, start={start}, end={end}'
        else:
            read.description = f'idx={read.id}, strand=-, start={stop-start}, end={stop-end}'

        read.letter_annotations = {'phred_quality': [50] * len(read)}
        reads_list.append(read)

    # return reads_list


def main(args):
    reference_path = os.path.abspath(args.ref)
    out_path = os.path.abspath(args.out)
    assert re.findall('.*\.([A-Za-z]*)', out_path)[-1] in ('fastq', 'fq'), \
            "Output path should be in the FASTQ format"
    depth = args.depth
    length_mean = args.length_mean
    length_std = args.length_std if args.length_std is not None else length_mean * 0.075

    subs = args.subs
    indel = args.indel

    reference = next(SeqIO.parse(reference_path, 'fasta'))
    reference_rc = reference.reverse_complement()

    step_mean, step_std = get_step(len(reference), length_mean, depth)
    reads_list = []

    g_error = introduce_errors(subs, indel)
    # Sample positive and negative strand
    g_seq = RandomSequenceGenerator(length_mean, length_std, step_mean, step_std, len(reference))
    sample_strand(reference, reads_list, '+', g_seq, g_error)
    g_seq = RandomSequenceGenerator(length_mean, length_std, step_mean, step_std, len(reference))
    sample_strand(reference_rc, reads_list, '-', g_seq, g_error)
    SeqIO.write(reads_list, out_path, 'fastq')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Takes reference path and path where to store the generated reads.')
    parser.add_argument('--ref', type=str, default='debug/reference.fasta', help='reference path')
    parser.add_argument('--out', type=str, default='debug/reads.fastq', help='path where to store the reads')
    parser.add_argument('--length-mean', type=int, default=20000, help='mean length of the simulated reads')
    parser.add_argument('--length-std', type=int, help='standard deviation in length of the simulated reads')
    parser.add_argument('--depth', type=int, default=20, help='sequencing depth to be simulated')
    parser.add_argument('--subs', type=float, default=0.0, help='substitution percentage')
    parser.add_argument('--indel', type=float, default=0.0, help='indel percentage')
    args = parser.parse_args()
    main(args)

