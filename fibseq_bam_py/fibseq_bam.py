#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""read fiberseq bam file and export tables"""

from os import path, makedirs
from timeit import default_timer as timer

import pysam
import numpy as np
# import decode_m6A_events
from bisect import bisect_left


'''tabix $BED chrX:49200000-49300000 | test.pl $FA chrX 49200000 49300000  m6A.D2.rest.FOXP3_neighborhood/ '''

fa_file = '/home/ehaugen/refseq/hg38/hg38.fa'
fa_file = '../data/hg38/hg38.fa'
bed_file = '/net/seq/data2/projects/ehaugen/fiberseq/2023mar08/indexed_bed/m6A.D2.rest.hg38.bed.gz'
bed_file = '../data/chr21_chr22_resdir_bed.aligned.m6a.bed.gz'
# bed_file = '../data/chrX_10Mbp_resdir_bed.aligned.m6a.bed.gz'
bam_file = '../data/m64323e_220528_041330_resdir.fiberseq.chrX_10Mb.chrX_40000000_50000000.bam'
region = 'chr21:0-45000000'
# region = 'chrX:40000000-50000000'
# region = 'chrX:49200000-49300000'
outputFolder = './m6A.D2.rest.FOXP3_neighborhood_new'


def fibseq_bam():
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
    for i in range(expectedReferenceLength):
        if referenceString[i].upper() in ['A', 'T']:
            refBasesAT.append(highlightMin0 + i)
    countATsInRange = len(refBasesAT)
    print(countATsInRange)

    hashes = []
    # vectorOffsets = []

    tbx = pysam.TabixFile(bed_file)
    records = tbx.fetch(region1, parser=pysam.asTuple())
    # samfile = pysam.AlignmentFile(bam_file, "rb")
    # records = samfile.fetch(region=region1)
    # print(sum(1 for _ in records))
    for cnt, record in enumerate(records):
        # chrom, chromStart, chromEnd, name, coverage, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts = record
        chrom, chromStart, chromEnd, name, _, _, _, _, _, _, _, blockStarts = record
        chromStart, chromEnd = [int(x) for x in [chromStart, chromEnd]]
        # mappedLength = chromEnd - chromStart
        # print(chrom, chromStart, mappedLength)

        # Not relative locations, just literally keep track of every base
        # and remember to ignore the end two positions of the fiber record
        starts = blockStarts.split(',')[1:-1]
        if not starts:
            continue
        hash_m6A = dict([(int(k), 1) for k in starts])
        # my %hash_m6A = map { $_ => 1 } @starts;

        sortvector = []
        # for refBase0 in refBasesAT:
        idx0 = bisect_left(refBasesAT, chromStart + 1)
        idx1 = bisect_left(refBasesAT, chromEnd - 1)
        for refBase0 in refBasesAT[idx0:idx1]:
            # Ignoring the BED12 end positions that aren't considered in the track
            if refBase0 > chromStart + 1 and refBase0 < chromEnd - 1:
                offset = refBase0 - chromStart
                methylated = 1.0 if offset in hash_m6A else 0.0
                sortvector.append(methylated)
            else:
                # Will include with NA's for missing extent
                # out of range for this molecule, incomplete overlap
                sortvector.append(float('nan'))
        # print(sortvector)

        # Skip this if no A/T bases in range
        if not len(sortvector):
            continue

        mean = np.nanmean(sortvector)
        hash = {
            'name': name,
            'line': record,
            'vector': sortvector,
            'inputorder': cnt + 1,
            'meanmeth': mean
        }
        # print(hash)
        # unless (defined($countATsInRange)) {
        #      $countATsInRange = scalar(@sortvector);
        # }
        hashes.append(hash)
        # if cnt == 0:
        #     break

    makedirs(outputFolder, exist_ok=True)

    # Matrix of the m6A statuses within excerpt region
    fileMatrix = path.join(outputFolder, 'matrix_{}_{}_{}.txt'.format(chromosome, highlightMin0, highlightMax1))
    out_matrix = open(fileMatrix, 'w')
    matrixHeader = '\t'.join(['ID'] + [str(x) for x in refBasesAT]) + '\n'
    out_matrix.write(matrixHeader)

    # Original m6A BED12 lines, just sorted (should I redo them with the header?)
    fileSorted = path.join(outputFolder, 'included_m6A_bed_tracks_sorted_by_methylation.txt')
    out_sorted = open(fileSorted, 'w')

    # Same order for outputs
    methsorted = sorted(hashes, key=lambda x: x.get('meanmeth'), reverse=True)
    for hashref in methsorted:
        # print(hashref['meanmeth'])
        out_sorted.write('{}\n'.format('\t'.join(hashref['line'])))
        out_matrix.write('{}\t{}\n'.format(hashref['name'], '\t'.join([str(x) for x in hashref['vector']])))

    print('Record count (min, max): {} ({}, {})'.format(len(methsorted), methsorted[-1]['meanmeth'], methsorted[0]['meanmeth']))

    out_sorted.close()
    out_matrix.close()

    print('Completed : {:.1f} min'.format((timer() - start) / 60))

if __name__ == '__main__':
    fibseq_bam()
