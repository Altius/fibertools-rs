---
---
`fibertools-rs/fiber_bam_py`
==============

`fibertools-rs` a CLI tool for creating and interacting with fiberseq bam files.
[fiberseq-rs readme](../README.md)

`fibertools-rs/fiber_bam_py` is a set of python scripts that provide stand-alone python for
interacting with fiberseq bam files or can use build `ft extract` in `fibertools-rs`
to provide identical capabilities with a ~4x performance boost.

While these scripts alone may provide useful output, it's expected that they will be starting
point for extending their analysis and output capabilities.  If you are modifying these scripts
directly instead of cloning for a new tool, please clearly note what you have added.

Available scripts
-----------------

### fibseq_bed.py

This is a python version of Eric Haugen's perl script to process fiberseq data from bed files
(`region_matrix_from_m6a_bed_2023apr19.pl`) with some performance enhancements.  It should
provide nearly identical output to that script with the same inputs.
This is here mainly for reference.

### fibseq_bam.py

This is a stand-alone python script to read fiberseq data from bam files and accurately
capture m6a and 5mC quality data.  It uses `pysam` for all the bam heavy lifting and then 
reads the appropriate tags to align methylation data with the reference genome.  There
are relatively simple outputs available currently.

Necessary python modules are listed in `requirements.txt`

slurm/hpc execution is available via `run_fiberseq_bam.sh`

### fibseq_rust.py

This is identical in functionality to `fibseq_bam.py` but calls `ft extract --all` to carry out the
reading and preliminary analysis of the bam file rather than python.  The output is then
extracted from that output.

This requires the version of `fibertools-rs` which is in this Altius fork.  It adds support for:
 - 5mC_qual values
 - region filtering
 - writing output to stdout

slurm/hpc execution is available via `run_fiberseq_rust.sh`



