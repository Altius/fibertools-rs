#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""read fiberseq bam file and export tables"""
import sys
from os import path, makedirs
from timeit import default_timer as timer
import subprocess

import pysam
import numpy as np
# import decode_m6A_events
from bisect import bisect_left


DEBUG = False
debugVerbose = False  # set to catch bad reads

# symbol representation
nucSymbol = 'n'       # '-'
mspSymbol = 'm'       # '+'
irrelSymbol = '.'     # '.'
emptySymbol = ' '     # ' '
unmarkedSymbol = '_'  # 'o'

input_file = None
region = None
outputFolder = './output'
fa_file = '/Volumes/photo2/fiberseq_data/hg38/hg38.fa'
# fa_file = '/home/ehaugen/refseq/hg38/hg38.fa'

# input_file = '/net/seq/data2/projects/Palladium_Dataset/stream/alld2rest.rerun.nokinetics.rg.bam.allchr.dipcall.haplotagged.HAP1.bam'
# region = 'chrX:0-45000000'

# input_file = '/net/seq/data2/projects/ehaugen/fiberseq/2023mar08/indexed_bed/m6A.D2.rest.hg38.bed.gz'
# region = 'chrX:49200000-49300000'

# input_file = '/Volumes/photo2/fiberseq_data/chr21_chr22_resdir_bed.aligned.m6a.bed.gz'
# input_file = '/Volumes/photo2/fiberseq_data/chr21_chr22_resdir_bed.aligned.m6a.bam'
# region = 'chr21:0-45000000'
# input_file = '/Volumes/photo2/fiberseq_data/chrX_10Mbp_resdir_bed.aligned.m6a.bed.gz'
# input_file = '/Volumes/photo2/fiberseq_data/chrX_10Mbp_resdir_bed.aligned.m6a.bam'
# region = 'chrX:40000000-50000000'
# input_file = '/Volumes/photo2/fiberseq_data/chrX_100Kb_resdir_bed.aligned.m6a.bed.gz'
# input_file = '/Volumes/photo2/fiberseq_data/chrX_100Kb_resdir_bed.aligned.m6a.bam'
# region = 'chrX:49200000-49300000'
# region = 'chrX:49202000-49203000'


def fibseq_bam():

    global input_file, region, outputFolder, fa_file

    if len(sys.argv) < 3:
        py_file = path.basename(sys.argv[0])
        print('Usage:\tpython {} input_file region [-o output_folder] [-f fa_file] [-s cpg,nuc,msp]'.format(py_file))
        print('      \t input file can be indexed .bed.gz or indexed .bam')
        print('      \t e.g. python {} chrX_100Kb_resdir_bed.aligned.m6a.bam chrX:49002000-49003000 -s cpg,nuc'.format(py_file))
        print('      \t   runs m6a output over the region specified and also shows cpg and nuc data')
        exit(0)

    if not input_file:
        input_file = sys.argv[1]
    if not region:
        region = sys.argv[2]
    if '-o' in sys.argv:
        outputFolder = sys.argv[sys.argv.index('-o') + 1]
    if '-f' in sys.argv:
        fa_file = sys.argv[sys.argv.index('-f') + 1]

    show_list = 'm6A,m5C,nuc,msp'
    if '-s' in sys.argv:
        show_list = sys.argv[sys.argv.index('-s') + 1]
    show_m6A = True  # always on
    show_m5C = True if 'm5C' in show_list else False
    show_nuc = True if 'nuc' in show_list else False
    show_msp = True if 'msp' in show_list else False

    start_time = timer()

    # get gene region
    tmp = pysam.FastaFile(fa_file)
    chromosome, highlight = region.split(':')
    highlightMin0, highlightMax1 = [int(x) for x in highlight.split('-')]
    highlightMin1 = highlightMin0 + 1
    region1 = '{}:{}-{}'.format(chromosome, highlightMin1, highlightMax1)
    referenceString = tmp.fetch(region=region1)
    # print(referenceString)

    expectedReferenceLength = highlightMax1 - highlightMin0
    foundReferenceLength = len(referenceString)
    print(foundReferenceLength)
    if foundReferenceLength != expectedReferenceLength:
        print('Expecting {} bp but read {}'.format(expectedReferenceLength, foundReferenceLength))
        exit(-1)

    # Mark the A's and T's that could get m6A calls
    if show_m6A:
        refBasesAT = []
        for i in range(expectedReferenceLength):
            if referenceString[i].upper() in ['A', 'T']:
                refBasesAT.append(highlightMin0 + i)
        countATsInRange = len(refBasesAT)
        print('AT count in region:{}'.format(countATsInRange))

    # Mark the C's and G's that could get m5C calls
    if show_m5C:
        refBasesCG = []
        for i in range(expectedReferenceLength):
            if referenceString[i].upper() in ['C', 'G']:
                refBasesCG.append(highlightMin0 + i)
        countCGsInRange = len(refBasesCG)
        print('CG count in region:{}'.format(countCGsInRange))

    hashes = []
    seg_len = highlightMax1 - highlightMin0

    records = []
    ext = path.splitext(input_file)[1]
    if ext not in ['.bam']:
        print('Invalid input file')
        exit(-1)
    # command_line = '../target/debug/ft extract {} --region {} --m6a stdout -r'.format(input_file, region1)
    # output = subprocess.getoutput(command_line)
    # m6a_bed_records = [x.split('\t') for x in output[:-1].split('\n')]
    # Test against bed values
    # blockStarts = m6a_bed_records[0][-1]
    # starts_bed = [int(x) for x in blockStarts.strip(',').split(',')[1:-1]]
    # temp = [x for x in m6a_starts if x not in starts_bed]

    exe_path = '/net/photo/photo1/Keith/ft'
    exe_path = '/Users/kgrochow/dev/fibertools-rs/target/debug/ft'
    command_line = '{} extract {} --region {} -a stdout -s -r'.format(exe_path, input_file, region1)
    output = subprocess.getoutput(command_line)
    output = output[output.index('#ct'):]  # strip warnings
    lines = [x.split('\t') for x in output[:-1].split('\n')]
    records = [dict(zip(lines[0], x)) for x in lines[1:]]
    for cnt, record in enumerate(records):
        chrom = record['#ct']
        chromStart = int(record['st'])
        chromEnd = int(record['en'])
        name = record['fiber']
        strand = record['strand']
        line = [chrom, chromStart, chromEnd, name, strand]

        sort_vector = []
        m6A_qual_vector = []
        if show_m6A and record['ref_m6a'] != '.':
            ref_starts = [int(x) - chromStart for x in record['ref_m6a'].strip(',').split(',')]
            quals = [int(x) for x in record['m6a_qual'].strip(',').split(',')]
            start_dicts = [{'r_pos': x, 'qual': quals[n]} for n, x in enumerate(ref_starts) if x != -1]

            if not start_dicts:
                continue
            hash_m6A = dict([(x['r_pos'], x['qual']) for x in start_dicts])

            m6A_qual_vector = [float('nan')] * seg_len
            idx0 = bisect_left(refBasesAT, chromStart + 2)
            idx1 = bisect_left(refBasesAT, chromEnd - 1)
            for refBase0 in refBasesAT[idx0:idx1]:
                if refBase0 > chromStart + 1 and refBase0 < chromEnd - 1:
                    offset = refBase0 - chromStart
                    if offset in hash_m6A:  # methylated
                        sort_vector.append(1.0)
                        m6A_qual_vector[refBase0 - highlightMin0] = float(hash_m6A[offset])
                    else:
                        sort_vector.append(0.0)
                        m6A_qual_vector[refBase0 - highlightMin0] = -1.0
                else:
                    # Will include with NA's for missing extent out of range for this molecule, incomplete overlap
                    print('should not get here')
                    sort_vector.append(float('nan'))
            # print(sort_vector)

        # Skip this if no A/T bases in range
        if not len(sort_vector):
            continue

        m5C_qual_vector = []
        if show_m5C and record['ref_5mC'] != '.':
            ref_starts = [int(x) - chromStart for x in record['ref_5mC'].strip(',').split(',')]
            quals = [int(x) for x in record['5mC_qual'].strip(',').split(',')]
            start_dicts = [{'r_pos': x, 'qual': quals[n]} for n, x in enumerate(ref_starts) if x != -1]

            if not start_dicts:
                continue
            hash_m5C = dict([(x['r_pos'], x['qual']) for x in start_dicts])

            m5C_qual_vector = [float('nan')] * seg_len
            idx0 = bisect_left(refBasesCG, chromStart + 2)
            idx1 = bisect_left(refBasesCG, chromEnd - 1)
            for refBase0 in refBasesCG[idx0:idx1]:
                if refBase0 > chromStart + 1 and refBase0 < chromEnd - 1:
                    offset = refBase0 - chromStart
                    if offset in hash_m5C:  # methylated
                        m5C_qual_vector[refBase0 - highlightMin0] = float(hash_m5C[offset])
                    else:
                        m5C_qual_vector[refBase0 - highlightMin0] = -1.0
                else:
                    # Will include with NA's for missing extent out of range for this molecule, incomplete overlap
                    print('should not get here')

        nuc_starts = []
        nuc_lengths = []
        if show_nuc and record['ref_nuc_starts'] != '.':
            nuc_starts = [int(x) - chromStart for x in record['ref_nuc_starts'].strip(',').split(',')]
            nuc_lengths = [int(x) for x in record['ref_nuc_lengths'].strip(',').split(',')]

        msp_starts = []
        msp_lengths = []
        if show_msp and record['ref_msp_starts'] != '.':
            msp_starts = [int(x) - chromStart for x in record['ref_msp_starts'].strip(',').split(',')]
            msp_lengths = [int(x) for x in record['ref_msp_lengths'].strip(',').split(',')]

        mean = np.nanmean(sort_vector)
        hash = {
            'chrom': chrom,
            'name': name,
            'line': line,
            'vector': sort_vector,
            'qual_m6A': m6A_qual_vector,
            'qual_m5C': m5C_qual_vector,
            'nuc_starts': nuc_starts,
            'nuc_lens': nuc_lengths,
            'msp_starts': msp_starts,
            'msp_lens': msp_lengths,
            'start': chromStart,
            'end': chromEnd,
            'strand': strand,
            'inputorder': cnt + 1,
            'meanmeth': mean
        }

        hashes.append(hash)

    print('Records processed: {}'.format(len(records)))

    makedirs(outputFolder, exist_ok=True)

    # Matrix of the m6A statuses within excerpt region
    file_matrix = path.join(outputFolder, 'matrix_{}_{}_{}.tsv'.format(chromosome, highlightMin0, highlightMax1))
    out_matrix = open(file_matrix, 'w')
    matrix_header = '\t'.join(['chrom', 'start', 'end', 'ID', 'strand', 'type', referenceString]) + '\n'
    out_matrix.write(matrix_header)

    # Original m6A BED12 lines, just sorted (should I redo them with the header?)
    file_sorted = path.join(outputFolder, 'included_m6A_bed_tracks_sorted_by_methylation.txt')
    out_sorted = open(file_sorted, 'w')
    included_header = '\t'.join(['chrom', 'start', 'end', 'ID', 'strand', 'mean_meth']) + '\n'
    out_sorted.write(included_header)

    # Save out the reports
    methsorted = sorted(hashes, key=lambda x: x.get('meanmeth'), reverse=True)
    for hashref in methsorted:
        out_sorted.write('{}\n'.format('\t'.join([str(x) for x in hashref['line'] + [hashref['meanmeth']]])))

        out_matrix.write(
            '{}\t{}\t{}\t{}\t{}\n'.format(hashref['name'], hashref['chrom'], hashref['start'], hashref['end'], hashref['strand']))

        chromStart = hashref['start']
        chromEnd = hashref['end']
        disp_string = [emptySymbol] * seg_len
        start = max(0, chromStart - highlightMin0)
        end = min(seg_len, chromEnd - highlightMin0)
        disp_string[start: end] = [irrelSymbol] * (end - start)

        def update_qual_string(qual_vector):
            my_disp = disp_string[:]
            for n, v in enumerate(qual_vector):
                if not np.isnan(v):
                    if v == -1.0:
                        my_disp[n] = unmarkedSymbol
                    else:
                        my_disp[n] = str(int(v) // 26)
            return ''.join(my_disp)

        def update_nuc_msp_string(ref_starts, ref_lens, symbol):
            my_disp = disp_string[:]
            for n, x in enumerate(ref_starts):
                start = chromStart - highlightMin0 + x
                end = min(seg_len, start + ref_lens[n])
                start = max(0, start)
                if start < seg_len and end > 0:
                    my_disp[start: end] = [symbol] * (end - start)
            return ''.join(my_disp)

        if show_m6A:
            m6A_disp = update_qual_string(hashref['qual_m6A'])
            out_matrix.write('\t\t\t\t\tm6A\t{}\n'.format(m6A_disp))
        if show_m5C:
            m5C_disp = update_qual_string(hashref['qual_m5C'])
            out_matrix.write('\t\t\t\t\tm5C\t{}\n'.format(m5C_disp))
        if show_nuc:
            nuc_disp = update_nuc_msp_string(hashref['nuc_starts'], hashref['nuc_lens'], nucSymbol)
            out_matrix.write('\t\t\t\t\tnuc\t{}\n'.format(nuc_disp))
        if show_msp:
            msp_disp = update_nuc_msp_string(hashref['msp_starts'], hashref['msp_lens'], mspSymbol)
            out_matrix.write('\t\t\t\t\tmsp\t{}\n'.format(msp_disp))

    print('Record count (min, max): {} ({}, {})'.format(len(methsorted),
                                methsorted[-1]['meanmeth'] if len(methsorted) else 0.0,
                                methsorted[0]['meanmeth'] if len(methsorted) else 0.0))

    print('Completed : {:.1f} min'.format((timer() - start_time) / 60))

    if methsorted:
        m6A_vals = np.array([x for x in hashref['qual_m6A'] for hashref in methsorted])
        m6A_vals = m6A_vals[~np.isnan(m6A_vals)]  # remove placeholder nan's
        m6A_vals = m6A_vals[~(m6A_vals == -1)]    # remove unmethylated
        m6a_bins = np.histogram(m6A_vals, 10, (0, 255))
        print('m6A bins {}'.format(m6a_bins[0]))

        m5C_vals = np.array([x for x in hashref['qual_m5C'] for hashref in methsorted])
        m5C_vals = m5C_vals[~np.isnan(m5C_vals)]
        m5C_vals = m5C_vals[~(m5C_vals == -1)]
        m5C_bins = np.histogram(m5C_vals, 10, (0, 255))
        print('m5C bins {}'.format(m5C_bins[0]))

    print('Completed : {:.1f} sec'.format((timer() - start_time)))

if __name__ == '__main__':
    fibseq_bam()
