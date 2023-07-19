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
        print('Usage:\t{} input_file region [-o output_folder] [-f fa_file]'.format(path.basename(sys.argv[0])))
        print('      \t input file can be indexed .bed.gz or indexed .bam')
        exit(0)

    if not input_file:
        input_file = sys.argv[1]
    if not region:
        region = sys.argv[2]
    if '-o' in sys.argv:
        outputFolder = sys.argv[sys.argv.index('-o') + 1]

    start = timer()

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

    ext = path.splitext(input_file)[1]
    if ext not in ['.bam', '.gz']:
        print('Invalid input file')
        exit(-1)
    if ext == '.bam':
        # command_line = '../target/debug/ft extract {} --region {} --m6a stdout -r'.format(input_file, region1)
        # output = subprocess.getoutput(command_line)
        # records = [x.split('\t') for x in output[:-1].split('\n')]
        samfile = pysam.AlignmentFile(input_file, "rb")
        records = samfile.fetch(region=region1)
    else:
        if not path.exists(input_file + '.tbi'):
            base_bed_file = path.splitext(input_file)[0]
            pysam.tabix_index(base_bed_file, preset='bed')
        tbx = pysam.TabixFile(input_file)
        records = tbx.fetch(region1, parser=pysam.asTuple())
    tot = 0
    for cnt, record in enumerate(records):
        tot = cnt + 1
        # _chrom, chromStart, chromEnd, name, coverage, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts = record
        # chrom, chromStart, chromEnd, name, _, _, _, _, _, _, _, blockStarts = record
        # Not relative locations, just literally keep track of every base
        # and remember to ignore the end two positions of the fiber record
        # starts = blockStarts.split(',')[1:-1]
        # line = record
        chrom = record.reference_name
        chromStart = record.reference_start
        chromEnd = record.reference_end
        name = record.query_name
        keys = [x for x in record.modified_bases.keys() if 'a' in x]
        base_starts = []
        for key in keys:
            base_starts += [{'q_pos': x, 'qual': v} for x, v in record.modified_bases[key]]
        sorted_starts = sorted(base_starts, key=lambda x: x['q_pos'])
        starts = realign_pos(record.aligned_pairs, sorted_starts)
        line = [record.reference_name, chromStart, chromEnd, name]

        if not starts:
            continue

        hash_m6A = dict([(x['r_pos'] - chromStart, x['qual']) for x in starts])

        chromStart, chromEnd = [int(x) for x in [chromStart, chromEnd]]
        # mappedLength = chromEnd - chromStart
        # print(chrom, chromStart, mappedLength)

        sort_vector = []
        qual_vector = []
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
                else:
                    sort_vector.append(0.0)
                    qual_vector.append(-1)
            else:
                # Will include with NA's for missing extent out of range for this molecule, incomplete overlap
                print('should not get here')
                sort_vector.append(float('nan'))
        # print(sort_vector)

        # Skip this if no A/T bases in range
        if not len(sort_vector):
            continue

        mean = np.nanmean(sort_vector)
        hash = {
            'chrom': chrom,
            'name': name,
            'line': line,
            'vector': qual_vector,
            'start': chromStart,
            'end': chromEnd,
            'idx0': idx0,
            'idx1': idx1,
            'inputorder': cnt + 1,
            'meanmeth': mean
        }
        # print(hash)

        hashes.append(hash)

    print('Records processed: {}'.format(tot))

    makedirs(outputFolder, exist_ok=True)
    compact_output = '-c' in sys.argv

    # Matrix of the m6A statuses within excerpt region
    fileMatrix = path.join(outputFolder, 'matrix_{}_{}_{}.tsv'.format(chromosome, highlightMin0, highlightMax1))
    out_matrix = open(fileMatrix, 'w')
    if compact_output:
        matrixHeader = '\t'.join(['chrom', 'start', 'end', 'ID', ''.join(refBases_string)]) + '\n'
    else:
        matrixHeader = '\t'.join(['ID'] + [str(x) for x in refBasesAT]) + '\n'
    out_matrix.write(matrixHeader)

    # Original m6A BED12 lines, just sorted (should I redo them with the header?)
    fileSorted = path.join(outputFolder, 'included_m6A_bed_tracks_sorted_by_methylation.txt')
    out_sorted = open(fileSorted, 'w')

    # Same order for outputs
    methsorted = sorted(hashes, key=lambda x: x.get('meanmeth'), reverse=True)
    for hashref in methsorted:
        # print(hashref['meanmeth'])
        out_sorted.write('{}\n'.format('\t'.join([str(x) for x in hashref['line']])))
        if compact_output:
            full_line = ['_'] * len(refBasesAT)
            # full_line[hashref['idx0']:hashref['idx1']] = ['.' if not x else '1' for x in hashref['vector']]
            full_line[hashref['idx0']:hashref['idx1']] = ['.' if x < 0 else str(x // 26) for x in hashref['vector']]
            out_matrix.write(
                '{}\t{}\t{}\t{}\t{}\n'.format(hashref['chrom'], hashref['start'], hashref['end'], hashref['name'],
                                              ''.join(full_line)))
        else:
            full_line = [float('nan')] * len(refBasesAT)
            full_line[idx0:idx1] = hashref['vector']
            out_matrix.write('{}\t{}\n'.format(hashref['name'], '\t'.join([str(x) for x in full_line])))

    print('Record count (min, max): {} ({}, {})'.format(len(methsorted), methsorted[-1]['meanmeth'],
                                                        methsorted[0]['meanmeth']))

    out_sorted.close()
    out_matrix.close()

    print('Completed : {:.1f} min'.format((timer() - start) / 60))


if __name__ == '__main__':
    fibseq_bam()
