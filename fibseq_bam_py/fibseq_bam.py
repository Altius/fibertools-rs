#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""read fiberseq bam file and export tables"""
import sys
from os import path, makedirs
from timeit import default_timer as timer

import pysam
import numpy as np
# import decode_m6A_events
from bisect import bisect_left
import re


input_file = None
region = None
outputFolder = './output'
# fa_file = '/Volumes/photo2/fiberseq_data/hg38/hg38.fa'
fa_file = '/home/ehaugen/refseq/hg38/hg38.fa'

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


def cigar2end(left, cigar):
    """Return right-most position of aligned read."""
    # store info about each CIGAR category
    cigar_pat = re.compile(r"\d+[MIDNSHP=X]{1}")
    counts = {"M": 0,  # M 0 alignment match (can be a sequence match or mismatch)
              'I': 0,  # I 1 insertion to the reference
              'D': 0,  # D 2 deletion from the reference
              'N': 0,  # N 3 skipped region from the reference
              'S': 0,  # S 4 soft clipping (clipped sequences present in SEQ)
              'H': 0,  # H 5 hard clipping (clipped sequences NOT present in SEQ)
              'P': 0,  # P 6 padding (silent deletion from padded reference)
              '=': 0,  # = 7 sequence match
              'X': 0,  # X 8 sequence mismatch
              }
    # split cigar entries
    for centry in cigar_pat.findall(cigar):
        ccount = int(centry[:-1])
        csymbol = centry[-1]
        counts[csymbol] = ccount
    # get number of aligned 'reference' bases
    aligned = counts["M"] + counts["D"] + counts["N"]  # + counts["="] + counts["X"]
    right = left + aligned
    return right


def realign_pos(aligned_pairs, pos_dicts):
    # find the shared positions in the reference
    ret = []
    cur_pos = 0
    for q_pos, r_pos in aligned_pairs:
        val_to_match = q_pos
        if val_to_match is None:
            continue
        # iterate over positions until we find the exact position or move past it
        while cur_pos < len(pos_dicts) and pos_dicts[cur_pos]['q_pos'] <= val_to_match:
            if pos_dicts[cur_pos]['q_pos'] == val_to_match and r_pos:
                ret_dict = {'q_pos': q_pos, 'r_pos': r_pos, 'qual': pos_dicts[cur_pos]['qual']}
                ret.append(ret_dict)
            cur_pos += 1
        if cur_pos == len(pos_dicts):
            break
    return ret


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
    print('Region bases: {}'.format(foundReferenceLength))
    if foundReferenceLength != expectedReferenceLength:
        print('Expecting {} bp but read {}'.format(expectedReferenceLength, foundReferenceLength))
        exit(-1)


    # Mark the A's and T's that could get m6A calls
    if show_m6a:
        refBasesAT = []
        for i in range(expectedReferenceLength):
            if referenceString[i].upper() in ['A', 'T']:
                refBasesAT.append(highlightMin0 + i)
        countATsInRange = len(refBasesAT)
        print('AT count in region:{}'.format(countATsInRange))

    # Mark the C's and G's that could get m6A calls
    if show_cpg:
        refBasesCG = []
        for i in range(expectedReferenceLength):
            if referenceString[i].upper() in ['C', 'G']:
                refBasesCG.append(highlightMin0 + i)
        countCGsInRange = len(refBasesCG)
        print('CG count in region:{}'.format(countCGsInRange))

    hashes = []
    disp_m6a = []
    disp_cpg = []
    disp_nuc = []
    disp_msp = []

    ext = path.splitext(input_file)[1]
    if ext not in ['.bam']:
        print('Invalid input file')
        exit(-1)

    samfile = pysam.AlignmentFile(input_file, "rb")
    records = samfile.fetch(region=region1)
    tot = 0
    for cnt, record in enumerate(records):
        tot += 1
        chrom = record.reference_name
        chromStart = record.reference_start
        chromEnd = record.reference_end
        name = record.query_name

        line = [chrom, chromStart, chromEnd, name]

        seg_len = highlightMax1 - highlightMin0
        disp_string = [' '] * seg_len
        start = max(0, chromStart - highlightMin0)
        end = min(seg_len, chromEnd - highlightMin0)
        disp_string[start: end] = ['.'] * (end - start)

        if show_m6a:
            keys = [x for x in record.modified_bases.keys() if 'a' in x]
            base_starts = []
            for key in keys:
                base_starts += [{'q_pos': x, 'qual': v} for x, v in record.modified_bases[key]]
            sorted_starts = sorted(base_starts, key=lambda x: x['q_pos'])
            starts = realign_pos(record.aligned_pairs, sorted_starts)
            hash_m6A = dict([(x['r_pos'] - chromStart, x['qual']) for x in starts])

            sort_vector = []
            qual_vector = []
            disp_m6a = disp_string[:]
            idx0 = bisect_left(refBasesAT, chromStart + 2)
            idx1 = bisect_left(refBasesAT, chromEnd - 1)
            for refBase0 in refBasesAT[idx0:idx1]:
            # for refBase0 in refBasesAT:
                # double check that ignoring the BED12 end positions that aren't considered in the track
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

        # Skip this if no A/T bases in range
        if not len(sort_vector):
            continue

        if show_cpg:
            keys = [x for x in record.modified_bases.keys() if 'm' in x]  # cpg change to 'm'
            base_starts = []
            for key in keys:
                base_starts += [{'q_pos': x, 'qual': v} for x, v in record.modified_bases[key]]
            sorted_starts = sorted(base_starts, key=lambda x: x['q_pos'])
            starts = realign_pos(record.aligned_pairs, sorted_starts)
            hash_cpg = dict([(x['r_pos'] - chromStart, x['qual']) for x in starts])

            disp_cpg = disp_string[:]
            idx0 = bisect_left(refBasesCG, chromStart + 2)
            idx1 = bisect_left(refBasesCG, chromEnd - 1)
            for refBase0 in refBasesCG[idx0:idx1]:
            # for refBase0 in refBasesCG:
                # double check that ignoring the BED12 end positions that aren't considered in the track
                if refBase0 > chromStart + 1 and refBase0 < chromEnd - 1:
                    offset = refBase0 - chromStart
                    if offset in hash_cpg:  # methylated
                        disp_cpg[refBase0 - highlightMin0] = str(hash_cpg[offset] // 26)
                    else:
                        disp_cpg[refBase0 - highlightMin0] = 'o'
                else:
                    # Will include with NA's for missing extent out of range for this molecule, incomplete overlap
                    print('should not get here')

        def get_tag(code):
            try:
                return record.tags[[x[0] for x in record.tags].index(code)][1]
            except:
                return None

        if show_nuc:
            nuc_starts = list(get_tag('ns'))
            nuc_lengths = list(get_tag('nl'))
            if record.is_reverse:
                nuc_ends = [x + nuc_lengths[n] for n, x in enumerate(nuc_starts)]
                nuc_starts = [seg_len - x - 1 for x in reversed(nuc_ends)]
                nuc_lengths.reverse()

            disp_nuc = disp_string[:]
            for n, x in enumerate(nuc_starts):
                start = chromStart - highlightMin0 + x
                end = min(seg_len, start + nuc_lengths[n])
                start = max(0, start)
                if start < seg_len and end > 0:
                    disp_nuc[start: end] = ['-'] * (end - start)

        if show_msp:
            msp_starts = list(get_tag('as'))
            msp_lengths = list(get_tag('al'))
            if record.is_reverse:
                msp_ends = [x + msp_lengths[n] for n, x in enumerate(msp_starts)]
                msp_starts = [seg_len - x - 1 for x in reversed(msp_ends)]
                msp_lengths.reverse()

            disp_msp = disp_string[:]
            for n, x in enumerate(msp_starts):
                start = chromStart - highlightMin0 + x
                end = min(seg_len, start + msp_lengths[n])
                start = max(0, start)
                if start < seg_len and end > 0:
                    disp_msp[start: end] = ['+'] * (end - start)

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

    print('Records processed: {}'.format(tot))

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

    out_sorted.close()
    out_matrix.close()

    print('Completed : {:.1f} sec'.format((timer() - start_time)))


if __name__ == '__main__':
    fibseq_bam()
