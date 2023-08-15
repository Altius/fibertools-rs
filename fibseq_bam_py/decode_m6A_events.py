#!/usr/bin/env python

import sys
import pysam
import numpy as np
import copy

chrom = sys.argv[1]
in_fn = sys.argv[2]

'''
§1.7 https://samtools.github.io/hts-specs/SAMtags.pdf

"Base modifications, including base methylation, are represented as a series 
of edits from the primary unmodified sequence as originally reported by the 
sequencing instrument. This potentially differs to the sequence stored in the 
main SAM SEQ field if the latter has been reverse complemented, in which case 
SAM FLAG 0x10 must be set. This means modification positions are also recorded 
against the original orientation (i.e. starting at the 5’ end), and count the 
original base types."

Calling x.get_forward_sequence() via pysam will retrieve the forward strand if 
the segment/read/molecule is forward-oriented. If reverse-oriented, this API 
call claims it will return the reverse complement of the sequence, which would 
give us the forward sequence. However, the documentation is wrong or worded
incorrectly and we need to take the reverse complement. Either way, we apply 
deltas to the sequence and treat unmodified bases as canonical to retrieve 
offsets from the reference 5' end.

Once we have these raw offsets, we use a processed form of operations in the 
CIGAR string to modify and filter offsets. For instance, offsets that align to
mismatch operations ('X') should be removed for purposes of later display or
aggregation. Insertion, deletion, and soft and hard padding operations also 
move events relative to the assembly reference sequence.
'''

ASSEMBLY = 'hg38'

ASSEMBLY_BOUNDS = {
  'hg38':{
    'chr1':{'ub':248956422},
    'chr10':{'ub':133797422},
    'chr11':{'ub':135086622},
    'chr12':{'ub':133275309},
    'chr13':{'ub':114364328},
    'chr14':{'ub':107043718},
    'chr15':{'ub':101991189}, 
    'chr16':{'ub':90338345},
    'chr17':{'ub':83257441},
    'chr18':{'ub':80373285},
    'chr19':{'ub':58617616},
    'chr2':{'ub':242193529},
    'chr20':{'ub':64444167},
    'chr21':{'ub':46709983},
    'chr22':{'ub':50818468},
    'chr3':{'ub':198295559},
    'chr4':{'ub':190214555},
    'chr5':{'ub':181538259},
    'chr6':{'ub':170805979},
    'chr7':{'ub':159345973},
    'chr8':{'ub':145138636},
    'chr9':{'ub':138394717},
    'chrX':{'ub':156040895},
    'chrY':{'ub':57227415},
  },
}

INCLUDE_CLIPPING_OPS = True

MO_SKELETON = {
  'unmodified_base' : '',
  'code' : '',
  'strand' : '',
  'offsets' : None,
  'probabilities' : None,
}

'''
correction for changing position indexing from zero-based 
to one-based coordinate system
'''
ONE_BASED_INDEX_CORRECTION = 1

'''
fast pure-Python reverse complement, courtesy of Devon Ryan
ref. https://bioinformatics.stackexchange.com/a/3585/776
'''
DNA_TABLE = str.maketrans("ACTG", "TGAC")
def reverse_complement(seq):
    return seq.translate(DNA_TABLE)[::-1]

def get_all_5p_indices(arr, val):
  indices = []
  for i in range(len(arr)):
    if (arr[i] == val):
      indices.append(i)
  return indices

def decode_raw_mm_ml_tags(mm_tag, ml_tag, seq):
  base_modifications = mm_tag.split(';')
  base_probabilities = ml_tag
  offset_objs = []
  current_offset_cnt = 0
  for bm in base_modifications:
    if len(bm) == 0: continue
    mo = copy.deepcopy(MO_SKELETON)
    elems = bm.split(',')
    n_modifications = len(elems) - 1
    mo['unmodified_base'] = elems[0][0]
    mo['strand'] = elems[0][1]
    mo['code'] = elems[0][2]
    mo['offsets'] = np.zeros((n_modifications,), dtype=np.uint32)
    mo['probabilities'] = np.zeros((n_modifications,), dtype=np.ubyte)
    base_indices = get_all_5p_indices(seq, mo['unmodified_base'])
    deltas = [int(x) for x in elems[1:len(elems)]]
    o = 0
    offset_idx = 0
    for i, d in enumerate(deltas):
      o += d
      try:
        base_offset = base_indices[o]
        base_probability = base_probabilities[i + current_offset_cnt]
        mo['offsets'][offset_idx] = base_offset
        mo['probabilities'][offset_idx] = base_probability
        o += 1
        offset_idx += 1
      except IndexError as err:
        break
    mo['offsets'].resize(offset_idx, refcheck=True)
    mo['probabilities'].resize(offset_idx, refcheck=True)
    offset_objs.append(mo)
    current_offset_cnt += n_modifications
  assert(current_offset_cnt == len(base_probabilities))
  return offset_objs

def convert_cigar_ops_to_substitution_objs(cigar_ops, include_clipping_ops):
  substitution_objs = []
  current_position = 0
  for (op, op_length) in cigar_ops:
    if op == pysam.CSOFT_CLIP and include_clipping_ops: # S
      substitution_objs.append({
        'pos': current_position,
        'length': op_length,
        'type': op,
        'letter': 'S',
      })
    elif op == pysam.CHARD_CLIP and include_clipping_ops: # H
      substitution_objs.append({
        'pos': current_position,
        'length': op_length,
        'type': op,
        'letter': 'H',
      })
    elif op == pysam.CDIFF: # X
      substitution_objs.append({
        'pos': current_position,
        'length': op_length,
        'type': op,
        'letter': 'X',
      })
      current_position += op_length
    elif op == pysam.CINS: # I
      substitution_objs.append({
        'pos': current_position,
        'length': op_length,
        'type': op,
        'letter': 'I',
      })
    elif op == pysam.CDEL: # D
      substitution_objs.append({
        'pos': current_position,
        'length': op_length,
        'type': op,
        'letter': 'D',
      })
      current_position += op_length
    elif op == pysam.CREF_SKIP: # N
      substitution_objs.append({
        'pos': current_position,
        'length': op_length,
        'type': op,
        'letter': 'N',
      })
      current_position += op_length
    elif op == pysam.CEQUAL: # =
      current_position += op_length
    elif op == pysam.CMATCH: # M
      current_position += op_length
    else: # skipping any other operation
      pass 
  return substitution_objs

def modify_raw_moos_with_substitution_objs(raw_moos, sos):
  modified_moos = []
  for mo in raw_moos:
    modified_mo = copy.deepcopy(MO_SKELETON)
    modified_mo['unmodified_base'] = mo['unmodified_base']
    modified_mo['strand'] = mo['strand']
    modified_mo['code'] = mo['code']
    n_offsets = mo['offsets'].size
    n_probabilities = n_offsets
    modified_mo['offsets'] = np.zeros((n_offsets,), dtype=np.uint32)
    modified_mo['probabilities'] = np.zeros((n_probabilities,), dtype=np.ubyte)
    offsets = mo['offsets']
    probabilities = mo['probabilities']
    offset_idx = 0
    offset_modifier = 0
    modified_mo_idx = 0
    for so in sos:
      if so['type'] == pysam.CSOFT_CLIP or so['type'] == pysam.CHARD_CLIP: # starting with 'S' or 'H'
        offset_modifier -= so['length']
        continue
      if offset_idx < n_offsets:
        while (offsets[offset_idx] + offset_modifier) < so['pos']:
          modified_mo['offsets'][modified_mo_idx] = offsets[offset_idx] + offset_modifier
          modified_mo['probabilities'][modified_mo_idx] = probabilities[offset_idx]
          modified_mo_idx += 1
          offset_idx += 1
          if offset_idx == n_offsets:
            break
      if offset_idx < n_offsets:
        if so['type'] == pysam.CDIFF and (offsets[offset_idx] + offset_modifier) == so['pos']: # 'X'
          offset_idx += 1
        elif so['type'] == pysam.CDEL: # 'D'
          offset_modifier += so['length']
        elif so['type'] == pysam.CINS: # 'I'
          offset_modifier -= so['length']
          offset_idx += 1
        elif so['type'] == pysam.CREF_SKIP: # 'N'
          offset_modifier += so['length']
        elif so['type'] == pysam.CSOFT_CLIP or so['type'] == pysam.CHARD_CLIP: # ending with 'S' or 'H'
          offset_modifier += so['length']
    modified_mo['offsets'].resize(modified_mo_idx, refcheck=True)
    modified_mo['probabilities'].resize(modified_mo_idx, refcheck=True)
    modified_moos.append(modified_mo)
  return modified_moos

m6A_posn_events = {}
m6A_posn_events[chrom] = {}
m6A_posn_molecules = {}

'''
API calls
https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment
'''
with pysam.AlignmentFile(in_fn, "rb") as in_fh:
  bam_start = 0
  bam_stop = ASSEMBLY_BOUNDS[ASSEMBLY][chrom]['ub']
  m6A_posn_molecules[chrom] = np.zeros((bam_stop,), dtype=np.uint16)
  iter = in_fh.fetch(chrom, bam_start, bam_stop)
  line_idx = 0
  for x in iter:
    start = x.reference_start # 0-based leftmost coordinate
    stop = x.reference_end
    for posn in range(start, stop):
      m6A_posn_molecules[chrom][posn] += 1
    mm_tag = x.get_tag('MM')
    ml_tag = x.get_tag('ML')
    seq = x.get_forward_sequence() if x.is_forward else reverse_complement(x.get_forward_sequence())
    raw_mm_offset_objs = decode_raw_mm_ml_tags(mm_tag, ml_tag, seq)
    cigar_ops = x.cigartuples
    substitution_objs = convert_cigar_ops_to_substitution_objs(cigar_ops, INCLUDE_CLIPPING_OPS)
    modified_mm_offset_objs = modify_raw_moos_with_substitution_objs(raw_mm_offset_objs, substitution_objs)
    for mmoo in modified_mm_offset_objs:
      if (mmoo['unmodified_base'] == 'A' and mmoo['strand'] == '+' and mmoo['code'] == 'a') \
        or (mmoo['unmodified_base'] == 'T' and mmoo['strand'] == '-' and mmoo['code'] == 'a'):
        for o in mmoo['offsets']:
          posn = start + o
          if posn not in m6A_posn_events[chrom]:
            m6A_posn_events[chrom][posn] = 0
          try:
            if posn < start or posn >= stop:
              m6A_posn_molecules[chrom][posn] += 1
          except IndexError as err: # ignore clipped bases
            pass
          m6A_posn_events[chrom][posn] += 1
    line_idx += 1

'''
molecule-normalized m6A event means
'''
for posn, events in sorted(m6A_posn_events[chrom].items(), key=lambda item: int(item[0])):
  try:
    mean_val = float(events) / float(m6A_posn_molecules[chrom][posn])
    start = int(posn) + ONE_BASED_INDEX_CORRECTION
    stop = start + 1
    out_line = '{}\t{}\t{}\t{:.6f}'.format(chrom, start, stop, mean_val)
    sys.stdout.write('{}\n'.format(out_line))
  except (KeyError, ZeroDivisionError, ValueError, IndexError) as err: # ignore clipped bases
    pass
