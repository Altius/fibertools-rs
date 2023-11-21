#!/usr/bin/env python
# -*- coding: utf-8 -*-
# HISTORY:  8/10/23  Sasha looks at Keith's script
# 				added +/- to chr name to denote strand;
# 				extra line break before m6A line so that it looks collinear in monospaced font
#               output log file with counts before and after conversions: immediately clear there is a problem with msp coordinates on (-) strand reads
#
#             8/11/23  output m6A values as number symbols
#               print strand as a separate field
#               changed "cpg" to "m5C" throughout (including the switch in the command line) to reflect what is actually being done
#               declared symbols at the start
#               filtering out weird reads that throw errors (see "Yikes")
#               TODO: this works fine (if slow) on *filtered* BAM files with full range of values for m6A and m5C [0..255]
#                       the BAM files are filtered for entries that contain "?" in MM:Z:C field:   samtools view -h in.bam | awk '!/MM:Z:C\+m\?/' | samtools view -bS - > filtered.bam
#                       it would be nice to have this filter here; otherwise, filtering takes a while

#
#		      8/12/23	fixed the issue with msp on (-) strand
#               summary stats are calculated; TODO: output to a file
#                         a question re summary stats: how to output the values
#             8/14/23   this branch is for output of m6A and m5C counts in detail
#                       file format is set needed to roll back the generate_vectors; it is slow
#             8/15/23
#             8/19/23   added file naming options
#             10/9/23   get the top 5% rather than top decile for m6A


"""read fiberseq bam file and export tables"""
# import cProfile
import sys
from os import path, makedirs
from timeit import default_timer as timer

import pysam
import numpy as np
import pandas as pd
# import decode_m6A_events
from bisect import bisect_left
import re

 # AAG - bunch of declarations
DEBUG        = False
debugVerbose = False  # set to catch bad reads
unMarked     = ','      #  'o'
nucSymbol    = '7'     # 'n'
mspSymbol    = '7'     # 'm'
irrelSymbol  = '.'         #'.'
emptySymbol  = ' '
#~ fullScore = False
prefix       = ''
prefix2      = ''


input_file = None
region = None
outputFolder = './output'
fa_file = '/home/ehaugen/refseq/hg38/hg38.fa'


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


def realign_posOLD(aligned_pairs, pos_dicts):
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

def realign_pos(aligned_pairs, pos_dicts):  # very slightly opt for speed
    ret = []
    cur_pos = 0
    pos_dicts_len = len(pos_dicts)

    for q_pos, r_pos in aligned_pairs:
        if q_pos is None or r_pos is None:
            continue

        while cur_pos < pos_dicts_len and pos_dicts[cur_pos]['q_pos'] <= q_pos:
            if pos_dicts[cur_pos]['q_pos'] == q_pos:
                ret.append({
                    'q_pos': q_pos,
                    'r_pos': r_pos,
                    'qual': pos_dicts[cur_pos]['qual']
                })
                cur_pos += 1
                break
            cur_pos += 1

    return ret


def fibseq_bam():
    global input_file, region, outputFolder, fa_file, prefix, prefix2

    if len(sys.argv) < 3:
        py_file = path.basename(sys.argv[0])
        print('Usage:\tpython {} input_file region [-o output_folder] [-f fa_file] [-s m5C,nuc,msp] [--fullScore]'.format(py_file))  # AAG
        print('      \t input file can be indexed .bed.gz or indexed .bam')
        print('      \t e.g. python {} chrX_100Kb_resdir_bed.aligned.m6A.bam chrX:49002000-49003000 -s m5C,nuc'.format(py_file))
        print('      \t   runs m6A output over the region specified and also shows m5C and nuc data')
        exit(0)

    if not input_file:
        input_file = sys.argv[1]
    if not region:
        region = sys.argv[2]
    if '-o' in sys.argv:
        outputFolder = sys.argv[sys.argv.index('-o') + 1]
    if '-f' in sys.argv:
        fa_file = sys.argv[sys.argv.index('-f') + 1]
    if '--prefix2' in sys.argv:
        prefix2 = sys.argv[sys.argv.index('--prefix2') + 1]
    if '--prefix' in sys.argv:
        prefix = sys.argv[sys.argv.index('--prefix') + 1]
    if '--fullScore' in sys.argv:   # AAG  # false for [244..255]   true for [0..255]
        fullScore = True
    else:
        fullScore = False

    show_m6A      = True  # always on for now
    show_list     = ''
    if '-s' in sys.argv:
        show_list = sys.argv[sys.argv.index('-s') + 1]
    show_m5C      = True if 'm5C' in show_list else False
    show_nuc      = True if 'nuc' in show_list else False
    show_msp      = True if 'msp' in show_list else False

    start_time = timer()



    # AAG: PART 1  here we are dealing with the reference sequence only
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

    # AAG create the methylation output file
    fileMeth = path.join(outputFolder, '{}_{}_{}_{}_{}.tsv'.format(prefix, prefix2, chromosome, highlightMin0, highlightMax1))
    out_meth = open(fileMeth, 'w')


    # Mark the A's and T's that could get m6A calls
    if show_m6A:
        refBasesAT = []
        for i in range(expectedReferenceLength):
            if referenceString[i].upper() in ['A', 'T']:
                refBasesAT.append(highlightMin0 + i)
        countATsInRange = len(refBasesAT)
        print('A/T count in region:{}'.format(countATsInRange))

    # Mark the C's and G's that could get modification calls
    if show_m5C:
        refBasesCG = []
        for i in range(expectedReferenceLength):
            if referenceString[i].upper() in ['C', 'G']:
                refBasesCG.append(highlightMin0 + i)
        countCGsInRange = len(refBasesCG)
        print('C/G count in region:{}'.format(countCGsInRange))   # AAG NOTE: this is (C or G) count, not CpG count


    # AAG Log file for debug
    if DEBUG:
        fileLog = path.join(outputFolder, 'matrix_{}_{}_{}.log'.format(chromosome, highlightMin0, highlightMax1))
        out_log = open(fileLog, 'w')
        #~ logHeader = '\t'.join(['ID', 'chrom', 'start', 'end', 'type']) + '\n'
        #~ out_log.write(logHeader)

    hashes   = []
    disp_m6A = []
    disp_m5C = []
    disp_nuc = []
    disp_msp = []

    ext = path.splitext(input_file)[1]
    if ext not in ['.bam']:
        print('Invalid input file')
        exit(-1)

    # AAG: PART 2 now we are starting to go through the BAM file   #####
    samfile = pysam.AlignmentFile(input_file, "rb")

    records = samfile.fetch(region=region1)
    tot = 0
    for cnt, record in enumerate(records):
        tot += 1
        chrom = record.reference_name
        chromStart = record.reference_start
        chromEnd = record.reference_end
        name = record.query_name
        if debugVerbose:
            print(cnt, name)

        strand = "+"   # AAG: added strand check
        if record.is_reverse:
            strand = "-"

        line = [record.reference_name, chromStart, chromEnd, name]

        seg_len = highlightMax1 - highlightMin0
        disp_string = [emptySymbol] * seg_len
        start = max(0, chromStart - highlightMin0)          # start >0 if read starts inside the window
        end = min(seg_len, chromEnd - highlightMin0)     # end =seg_len if read ends outside the window
        #~ disp_string[start: end] = ['.'] * (end - start)   # original
        disp_string[start: end] = [irrelSymbol] * (end - start)


        if show_m6A:  ## supposedly optimized
            keys = [x for x in record.modified_bases.keys() if 'a' in x]
            base_starts = [item for key in keys for item in [{'q_pos': x, 'qual': v} for x, v in record.modified_bases[key]]]
            starts = realign_pos(record.aligned_pairs, sorted(base_starts, key=lambda x: x['q_pos']))
            hash_m6A = {x['r_pos'] - chromStart: x['qual'] for x in starts}

            sort_vector = []
            qual_vector = []
            disp_m6A = disp_string[:]
            idx0 = bisect_left(refBasesAT, chromStart + 2)
            idx1 = bisect_left(refBasesAT, chromEnd - 1)

            for refBase0 in refBasesAT[idx0:idx1]:
                if chromStart + 1 < refBase0 < chromEnd - 1:
                    offset = refBase0 - chromStart
                    if offset in hash_m6A:  # methylated
                        sort_vector.append(1.0)
                        qual_vector.append(hash_m6A[offset])

                        # AAG
                        if fullScore:
                            disp_m6A[refBase0 - highlightMin0] = str(hash_m6A[offset] // 26)
                            # disp_m6A[refBase0 - highlightMin0] = chr(hash_m6A[offset])
                        else:
                            disp_m6A[refBase0 - highlightMin0] = chr(max(0, (hash_m6A[offset]-240)//2+48))
                    else:
                        sort_vector.append(0.0)
                        qual_vector.append(-1)
                        disp_m6A[refBase0 - highlightMin0] = unMarked
                else:
                    print('should not get here')
                    sort_vector.append(float('nan'))


        # Skip this if no A/T bases in range
        if not len(sort_vector):
            continue

        if show_m5C:   # supposedly optimized by ChatGPT
            # List comprehension is usually faster than a for-loop in Python.
            keys = [x for x in record.modified_bases.keys() if 'm' in x]

            # Using a list comprehension directly without additional appends
            base_starts = [{'q_pos': x, 'qual': v} for key in keys for x, v in record.modified_bases[key]]

            # The lambda is fairly straightforward so it's probably not a bottleneck, but if you find it is, consider using an itemgetter.
            sorted_starts = sorted(base_starts, key=lambda x: x['q_pos'])
            starts = realign_pos(record.aligned_pairs, sorted_starts)

            # Direct dictionary comprehension
            hash_m5C = {x['r_pos'] - chromStart: x['qual'] for x in starts}

            disp_m5C = disp_string[:]

            idx0 = bisect_left(refBasesCG, chromStart + 2)
            idx1 = bisect_left(refBasesCG, chromEnd - 1)

            # Limit the loop by using the indices
            for refBase0 in refBasesCG[idx0:idx1]:
                offset = refBase0 - chromStart

                # Checking the keys of the dictionary is generally fast
                if offset in hash_m5C:  # methylated
                        # AAG
                        if fullScore:
                            # disp_m5C[refBase0 - highlightMin0] = str(hash_m5C[offset] // 26)
                            disp_m5C[refBase0 - highlightMin0] = chr(hash_m5C[offset])
                            # if (hash_m5C[offset] < 100):
                                # print(" {}\t->\t{} ".format(hash_m5C[offset],disp_m5C[refBase0 - highlightMin0]))
                        else:
                            disp_m5C[refBase0 - highlightMin0] = chr(max(0, (hash_m5C[offset]-240)//2+48))
                        # disp_m5C[refBase0 - highlightMin0] = str(hash_m5C[offset] // 26)
                else:
                    disp_m5C[refBase0 - highlightMin0] = unMarked


        def get_tag(code):
            try:
                return record.tags[[x[0] for x in record.tags].index(code)][1]
            except:
                return None

        def get_tag_opt(code):  # optimized? no improvement on 10K chunk
            for tag, value in record.tags:
                if tag == code:
                    return value
            return None


        if show_nuc:    # NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
            seq_len = record.qlen
            nus=get_tag('ns')
            if nus == None:  # AAG: filter out reads w/o nucleosomes
                print("Yikes: a read w/o nucleosome call was skipped:", name)
                continue
            nuc_starts = list(nus)
            #~ print('nuc_starts:', nuc_starts)
            nuc_lengths = list(get_tag('nl'))
            if record.is_reverse:
                nuc_ends = [x + nuc_lengths[n] for n, x in enumerate(nuc_starts)]
                nuc_starts = [seq_len - x - 1 for x in reversed(nuc_ends)]
                nuc_lengths.reverse()

            if DEBUG:   # AAG  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                #logStr = ('{}\t{}\n'.format(hashref['disp_m5C']))
                logInfo = ','.join(map(str, nuc_starts)) + '\n'
                logStr =  ('{}\t{}\t'.format(name, strand)) + 'nuc_starts:\t'+ logInfo
                out_log.write(logStr)
                #~ logInfo = ','.join(map(str, nuc_lengths)) + '\n'
                #~ logStr =  ('{}\t{}\t'.format(name, strand)) + 'nuc_lengths:\t'+ logInfo
                #~ out_log.write(logStr)

                ns = np.array(nuc_starts)
                nl = np.array(nuc_lengths)
                ne =  ns+nl
                nuc_ends = list(ne)
                nuc_endd = np.array(nuc_starts)+np.array(nuc_lengths)

                logInfo = ','.join(map(str, nuc_ends)) + '\n'
                logStr =  ('{}\t{}\t'.format(name, strand)) + '*nuc_ends: \t'+ logInfo
                out_log.write(logStr)
                logInfo = ','.join(map(str, nuc_endd)) + '\n'
                logStr =  ('{}\t{}\t'.format(name, strand)) + '*numpy   : \t'+ logInfo
                out_log.write(logStr)



            disp_nuc = disp_string[:]
            for n, x in enumerate(nuc_starts):
                start = chromStart - highlightMin0 + x
                end = min(seg_len, start + nuc_lengths[n])
                start = max(0, start)
                if start < seg_len and end > 0:
                    disp_nuc[start: end] = [nucSymbol] * (end - start)

        if show_msp:   # MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
            mss=get_tag('as')
            if mss == None:  # AAG: filter out reads w/o nucleosomes
                print("Yikes: a read w/o msp calls was skipped:", name)
                continue
            msp_starts = list(mss)
            msp_lengths = list(get_tag('al'))   # AAG: hopefully no issue if we get here
            if DEBUG:   # AAG  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                #logStr = ('{}\t{}\n'.format(hashref['disp_m5C']))
                logInfo = ','.join(map(str, msp_starts)) + '\n'
                logStr =  ('{}\t{}\t'.format(name, strand)) + 'msp_starts:\t'+ logInfo
                out_log.write(logStr)
                logInfo = ','.join(map(str, msp_lengths)) + '\n\n'
                logStr =  ('{}\t{}\t'.format(name, strand)) + 'msp_lengths:\t'+ logInfo
                out_log.write(logStr)

            # AAG: ok, reverse strand here
            if record.is_reverse:
                msp_ends = [x + msp_lengths[n] for n, x in enumerate(msp_starts)]
                msp_starts = [seq_len - x - 1 for x in reversed(msp_ends)]   # AAG seg_len is typo?? - YES; changing it fixed the problem
                msp_lengths.reverse()
                if DEBUG:   # AAG  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
                #logStr = ('{}\t{}\n'.format(hashref['disp_m5C']))
                    logInfo = ','.join(map(str, msp_starts)) + '\n'
                    logStr  =  ('{}\t{}\t'.format(name, strand)) + '#msp_starts:\t'+ logInfo
                    out_log.write(logStr)
                    logInfo = ','.join(map(str, msp_lengths)) + '\n\n'
                    logStr  =  ('{}\t{}\t'.format(name, strand)) + '#msp_lengths:\t'+ logInfo
                    out_log.write(logStr)


            disp_msp = disp_string[:]
            for n, x in enumerate(msp_starts):
                start = chromStart - highlightMin0 + x
                end = min(seg_len, start + msp_lengths[n])
                start = max(0, start)
                if start < seg_len and end > 0:
                    disp_msp[start: end] = [mspSymbol] * (end - start)

        mean = np.nanmean(sort_vector)
        hash = {
            'chrom': chrom,
            'strand': strand,   # AAG
            'name': name,
            'line': line,
            'vector': sort_vector,
            'disp_m6A': ''.join(disp_m6A),
            'disp_m5C': ''.join(disp_m5C),
            'disp_nuc': ''.join(disp_nuc),
            'disp_msp': ''.join(disp_msp),
            'start': chromStart,
            'end': chromEnd,
            'idx0': idx0,
            'idx1': idx1,
            'inputorder': cnt + 1,
            'meanmeth': mean
        }
        # print(hash)

        hashes.append(hash)

    # AAG: PART 3 OK, we are done going through the BAM file ###############################################
    print('Records processed: {}'.format(tot))

    makedirs(outputFolder, exist_ok=True)

    if not fullScore:   # TODO: add prefix to file names
        # Matrix of the m6A statuses within excerpt region
        fileMatrix = path.join(outputFolder, 'matrix_{}_{}_{}.tsv'.format(chromosome, highlightMin0, highlightMax1))
        out_matrix = open(fileMatrix, 'w')
        matrixHeader = '\t'.join(['ID', 'chrom', 'start', 'end', 'type', referenceString]) + '\n'
        out_matrix.write(matrixHeader)

        # Original m6A BED12 lines, just sorted (should I redo them with the header?)
        fileSorted = path.join(outputFolder, 'included_m6A_bed_tracks_sorted_by_methylation.txt')
        out_sorted = open(fileSorted, 'w')

        # Same order for outputs
        # AAG: everything lives in hashes
        methsorted = sorted(hashes, key=lambda x: x.get('meanmeth'), reverse=True)  # AAG: do we want sorted by methylation? does it slow things down?
        for hashref in methsorted:
            # print(hashref['meanmeth'])
            out_sorted.write('{}\n'.format('\t'.join([str(x) for x in hashref['line']])))

            out_matrix.write(  # AAG: add break
                '{}\t{}\t{}\t{}\t{}\n\t\t\t\tm6A\t{}\n'.format(hashref['name'], hashref['chrom'], hashref['start'], hashref['end'], hashref['strand'],  # AAG add strand
                                              hashref['disp_m6A']))
            if show_m5C:
                out_matrix.write('\t\t\t\tm5C\t{}\n'.format(hashref['disp_m5C']))
            if show_nuc:
                out_matrix.write('\t\t\t\tnuc\t{}\n'.format(hashref['disp_nuc']))
            if show_msp:
                out_matrix.write('\t\t\t\tmsp\t{}\n'.format(hashref['disp_msp']))

        print('Record count (min, max m6A): {} ({}, {})'.format(len(methsorted),
                                    methsorted[-1]['meanmeth'] if len(methsorted) else 0.0,
                                    methsorted[0]['meanmeth'] if len(methsorted) else 0.0))

        out_sorted.close()
        out_matrix.close()

    print('BAM ingestion : {:.1f} sec'.format((timer() - start_time)))

    ###########################################################################################

    # AAG: PART 4 now calculate summary stats across the records     DDDDDDDDDDDDDDDDDDDDDDDDDD
    start_time = timer()
    # open out file
    out_meth.write("# {}:{}-{}    assembly used:{}\n".format(chromosome, highlightMin0, highlightMax1, fa_file))
    out_meth.write("# nt\tmA\tA\tunkA\tmC\tC\tunkC\n")

    # Convert list of dictionaries to DataFrame
    df = pd.DataFrame(hashes)
    if df.empty:
        print("Warning: empty df at {}:{}-{}".format(chromosome, highlightMin0, highlightMax1))
        out_meth.write("## EMPTY RECORDS")

        exit()
    print('To df : {:.1f} sec'.format((timer() - start_time)))
    start_time = timer()


    # DOESNT WORK !! much faster: about 0.1 sec on 100K chunk for two data types
    # AAG TODO: check that actual conversions are correct; Also will need to change values for nuc and msp so that the same function can be used

    def generate_vectors(column):
        # Convert the column to numeric representation
        numeric_representation = column.apply(list).apply(pd.Series)

        # Define the conditions
        cond_yes = numeric_representation == '9'
        cond_no = numeric_representation == '0'
        cond_unknown_or_unmarked = numeric_representation.isin([str(i) for i in range(1, 8)] + [unMarked])

        # Count the conditions across the records
        vector_1 = cond_yes.sum()
        vector_2 = cond_no.sum()
        vector_3 = cond_unknown_or_unmarked.sum()

        return vector_1, vector_2, vector_3


    def generate_vectors_2(column):
        # Convert the column to numeric representation
        numeric_representation = column.apply(list).values
        numeric_representation = np.vstack(numeric_representation)

        # Count the conditions across the records
        vector_1 = np.sum(numeric_representation == '9', axis=0)
        vector_2 = np.sum(numeric_representation == '0', axis=0)
        # vector_3 = np.sum((numeric_representation >= '1') & (numeric_representation <= '8'), axis=0)
        vector_3 = np.sum((numeric_representation >= '1') & (numeric_representation <= '8') | (numeric_representation == unMarked), axis=0)


        return vector_1, vector_2, vector_3

    # Generate vectors for both columns

    m6A_vectors = generate_vectors_2(df['disp_m6A'])
    m5C_vectors = generate_vectors_2(df['disp_m5C'])   #

    # print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX:")
    # print (m6A_vectors)
    # print("cond_0:", m6A_vectors[1])

    print('generated vectors : {:.1f} sec'.format((timer() - start_time)))
    start_time = timer()


    # Convert the vectors to numpy arrays
    m6A_np_arrays = [np.array(vec) for vec in m6A_vectors]
    m5C_np_arrays = [np.array(vec) for vec in m5C_vectors]

    # Vertically stack the arrays
    stacked_arrays = np.vstack((*m6A_np_arrays, *m5C_np_arrays))


    # Iterate over rows for printing
    # print("Meth output:\n")
    for i in range(stacked_arrays.shape[1]):  # assuming all vectors have the same length
        # print("\t".join(str(stacked_arrays[j, i]) for j in range(3)))
        # print("\t".join(str(stacked_arrays[j+3, i]) for j in range(3)))
        out_meth.write(referenceString[i]+'\t')
        out_meth.write("\t".join(str(stacked_arrays[j, i]) for j in range(3)) +'\t')
        out_meth.write("\t".join(str(stacked_arrays[j+3, i]) for j in range(3)))
        out_meth.write("\n")
    out_meth.close()

    print('numpy zip and print: {:.1f} sec'.format((timer() - start_time)))

if __name__ == '__main__':
    fibseq_bam()
    # cProfile.run('fibseq_bam()', sort='tottime')

