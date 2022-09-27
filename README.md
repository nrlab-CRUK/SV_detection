# Genomic rearrangement detection by fuzzy matching of junction flanking sequences

[`find_junction_spanning_sequences.R`](find_junction_spanning_sequences.R)
is an R script that identifies genomic sequencing reads containing junctions
from structural rearrangements (also known as structural variants or SVs) by
fuzzy matching of sequences that flank either side of the junction. It was
written to support a study exploring different approaches for detecting known
structural variants in tumours from breast cancer patients, especially to
monitor disease progression or recurrence over time.

TODO: reference and link to online paper to go here

This was originally written to find SV-supporting reads from a multiplex PCR
assay in which the PCR primer pairs are the flanking sequences for which the
fuzzy matching search is run. The script searches for the primer pairs for the
specified SV junctions within all sequence reads in a FASTQ file, allowing for
a specified number of mismatches in each primer sequence. It reports all matches
and, if provided with the complete amplicon or junction sequence will also
compute the total number of mismatches (edit distance).

It has also been used on whole genome sequencing (WGS) data to find the very
small subset of reads that support known or expected structural variants. In
this case the sequences flanking the junction need to be provided in the same
way that the primer pair sequences were used for the multiplex PCR assay. The
fuzzy matching operations are computationally expensive and running these on all
reads within a high-depth WGS dataset is not practical. Instead the approach we
have taken in the breast cancer study is to extract soft-clipped reads from the
sequence data aligned to the human reference genome as junction-spanning reads
will align to the two different parts of the genome and will have clipped
alignments. This produces a subset of reads on which the R script can be run
using significantly less computing resources. This approach is implemented as
a workflow written using the [Nextflow](https://www.nextflow.io) tool for
data-driven computational pipelines. The workflow is provided here in the
[`junction_detection.nf`](junction_detection.nf) file and instructions on how
to use this are given below.

## Running the R script on multiplex PCR sequencing data


