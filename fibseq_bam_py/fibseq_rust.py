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
input_file = '/Volumes/photo2/fiberseq_data/chrX_100Kb_resdir_bed.aligned.m6a.bam'
# region = 'chrX:49200000-49300000'
region = 'chrX:49202000-49203000'


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

    show_m6a = True  # always on for now
    show_list = ''
    if '-s' in sys.argv:
        show_list = sys.argv[sys.argv.index('-s') + 1]
    show_cpg = True if 'cpg' in show_list else False
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
    refBasesAT = []
    refBases_string = []
    for i in range(expectedReferenceLength):
        if referenceString[i].upper() in ['A', 'T']:
            refBasesAT.append(highlightMin0 + i)
            refBases_string.append(referenceString[i].upper())
    countATsInRange = len(refBasesAT)
    print(countATsInRange)

    hashes = []

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

    command_line = '../target/debug/ft extract {} --region {} -a stdout -s -r'.format(input_file, region1)
    output = subprocess.getoutput(command_line)
    lines = [x.split('\t') for x in output[:-1].split('\n')]
    records = [dict(zip(lines[0], x)) for x in lines[1:]]
    for cnt, record in enumerate(records):
        chrom = record['#ct']
        chromStart = int(record['st'])
        chromEnd = int(record['en'])
        name = record['fiber']
        line = [chrom, chromStart, chromEnd, name]

        seg_len = highlightMax1 - highlightMin0
        disp_string = [' '] * seg_len
        start = max(0, chromStart - highlightMin0)
        end = min(seg_len, chromEnd - highlightMin0)
        disp_string[start: end] = ['.'] * (end - start)

        if show_m6a:
            ref_starts = [int(x) - chromStart for x in record['ref_m6a'].strip(',').split(',')]
            quals = [int(x) for x in record['m6a_qual'].strip(',').split(',')]
            start_dicts = [{'r_pos': x, 'qual': quals[n]} for n, x in enumerate(ref_starts) if x != -1]

            if not start_dicts:
                continue
            hash_m6A = dict([(x['r_pos'], x['qual']) for x in start_dicts])

            sort_vector = []
            qual_vector = []
            disp_m6a = disp_string[:]
            idx0 = bisect_left(refBasesAT, chromStart + 2)
            idx1 = bisect_left(refBasesAT, chromEnd - 1)
            for refBase0 in refBasesAT[idx0:idx1]:
                if refBase0 > chromStart + 1 and refBase0 < chromEnd - 1:
                    offset = refBase0 - chromStart
                    if offset in hash_m6A:  # methylated
                        sort_vector.append(1.0)
                        qual_vector.append(hash_m6A[offset])
                        disp_m6a[refBase0 - highlightMin0] = str(hash_m6A[offset] // 26)
                    else:
                        sort_vector.append(0.0)
                        qual_vector.append(-1)
                        disp_m6a[refBase0 - highlightMin0] = 'o'
                else:
                    # Will include with NA's for missing extent out of range for this molecule, incomplete overlap
                    print('should not get here')
                    sort_vector.append(float('nan'))
            # print(sort_vector)

        # Skip this if no A/T bases in range
        if not len(sort_vector):
            continue

        if show_cpg:
            ref_starts = [int(x) - chromStart for x in record['ref_m6a'].strip(',').split(',')]
            quals = [int(x) for x in record['m6a_qual'].strip(',').split(',')]
            start_dicts = [{'r_pos': x, 'qual': quals[n]} for n, x in enumerate(ref_starts) if x != -1]

            if not start_dicts:
                continue
            hash_cpg = dict([(x['r_pos'], x['qual']) for x in start_dicts])

            disp_cpg = disp_string[:]
            idx0 = bisect_left(refBasesAT, chromStart + 2)
            idx1 = bisect_left(refBasesAT, chromEnd - 1)
            for refBase0 in refBasesAT[idx0:idx1]:
                if refBase0 > chromStart + 1 and refBase0 < chromEnd - 1:
                    offset = refBase0 - chromStart
                    if offset in hash_cpg:  # methylated
                        disp_cpg[refBase0 - highlightMin0] = str(hash_cpg[offset] // 26)
                    else:
                        disp_cpg[refBase0 - highlightMin0] = 'o'
                else:
                    # Will include with NA's for missing extent out of range for this molecule, incomplete overlap
                    print('should not get here')

        if show_nuc:
            ref_starts = [int(x) - chromStart for x in record['ref_nuc_starts'].strip(',').split(',')]
            ref_lens = [int(x) for x in record['ref_nuc_lengths'].strip(',').split(',')]
            disp_nuc = disp_string[:]
            for n, x in enumerate(ref_starts):
                start = chromStart - highlightMin0 + x
                end = min(seg_len, start + ref_lens[n])
                start = max(0, start)
                if start < seg_len and end > 0:
                    disp_nuc[start: end] = ['-'] * (end - start)

        if show_msp:
            ref_starts = [int(x) - chromStart for x in record['ref_msp_starts'].strip(',').split(',')]
            ref_lens = [int(x) for x in record['ref_msp_lengths'].strip(',').split(',')]
            disp_msp = disp_string[:]
            for n, x in enumerate(ref_starts):
                start = chromStart - highlightMin0 + x
                end = min(seg_len, start + ref_lens[n])
                start = max(0, start)
                if start < seg_len and end > 0:
                    disp_msp[start: end] = ['-'] * (end - start)

        mean = np.nanmean(sort_vector)
        hash = {
            'chrom': chrom,
            'name': name,
            'line': line,
            'vector': sort_vector,
            'disp_m6a': ''.join(disp_m6a),
            'disp_cpg': ''.join(disp_cpg),
            'disp_nuc': ''.join(disp_nuc),
            'disp_msp': ''.join(disp_msp),
            'start': chromStart,
            'end': chromEnd,
            'inputorder': cnt + 1,
            'meanmeth': mean
        }
        # print(hash)

        hashes.append(hash)

    print('Records processed: {}'.format(len(records)))

    makedirs(outputFolder, exist_ok=True)

    # Matrix of the m6A statuses within excerpt region
    fileMatrix = path.join(outputFolder, 'matrix_{}_{}_{}.tsv'.format(chromosome, highlightMin0, highlightMax1))
    out_matrix = open(fileMatrix, 'w')
    matrixHeader = '\t'.join(['ID', 'chrom', 'start', 'end', 'type', referenceString]) + '\n'
    out_matrix.write(matrixHeader)

    # Original m6A BED12 lines, just sorted (should I redo them with the header?)
    fileSorted = path.join(outputFolder, 'included_m6A_bed_tracks_sorted_by_methylation.txt')
    out_sorted = open(fileSorted, 'w')

    # Same order for outputs
    methsorted = sorted(hashes, key=lambda x: x.get('meanmeth'), reverse=True)
    for hashref in methsorted:
        # print(hashref['meanmeth'])
        out_sorted.write('{}\n'.format('\t'.join([str(x) for x in hashref['line']])))

        out_matrix.write(
            '{}\t{}\t{}\t{}\tm6a\t{}\n'.format(hashref['name'], hashref['chrom'], hashref['start'], hashref['end'],
                                          hashref['disp_m6a']))
        if show_cpg:
            out_matrix.write('\t\t\t\tcpg\t{}\n'.format(hashref['disp_cpg']))
        if show_nuc:
            out_matrix.write('\t\t\t\tnuc\t{}\n'.format(hashref['disp_nuc']))
        if show_msp:
            out_matrix.write('\t\t\t\tmsp\t{}\n'.format(hashref['disp_msp']))

    print('Record count (min, max): {} ({}, {})'.format(len(methsorted),
                                methsorted[-1]['meanmeth'] if len(methsorted) else 0.0,
                                methsorted[0]['meanmeth'] if len(methsorted) else 0.0))

    print('Completed : {:.1f} min'.format((timer() - start_time) / 60))

if __name__ == '__main__':
    fibseq_bam()
