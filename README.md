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

## Running the `find_junction_spanning_sequences.R` script on multiplex PCR sequencing data

### Dependencies

The `find_junction_spanning_sequences.R` script is an R script and requires
both R and the following packages to be installed: `readr`, `dplyr`, `purrr`,
`stringr` and `stringi`. The script has been run using R version 4.2.0 as part
of the breast cancer study described above. The dependent packages are part of
the [`tidyverse`](https://www.tidyverse.org) and can be installed from the R
command prompt with the following command:

```
install.packages("tidyverse")
```

### Inputs

The script operates on a single FASTQ file and iterates over sequence records
contained within looking for matches to junction/flanking sequences given in a
CSV file.

The flanking sequences CSV file must contain `ID`, `LeftFlankingSequence` and
`RightFlankingSequence` columns and may optionally contain a `JunctionSequence`
column that contains the entire junction-spanning sequence that includes both
flanking sequences. For multiplexed PCR data in which the sequence reads are
from multiple polymerase chain reactions each targeting a different junction
amplicon, the left and right flanking sequences will be the PCR primer pair.

An example of the flanking sequences CSV file used in the breast cancer study,
truncated to show just the first 2 SV junctions, is given below:

```
ID,LeftFlankingSequence,RightFlankingSequence,JunctionSequence
PBCP_001_R002,ATCAGGGTAAAGCCAGAGCT,CCTTCTTATATAGACACCAG,TGTGGCCAGTCACTATGGAAACATCAAGCTGGTGAAGTTTCTGCTGCAGCACCAGGCAGATGTCAATGCCAAGACCAAGGTACAGGGGTGCCCCAGCCCCGACTCCTGCACTGACCCTTCTCCATGCCACATCAGGGTAAAGCCAGAGCTGCTCCCCTTCTTATATAGACACCAGTCATATAGAATTAGGGCCCCACCCTTATGATCGTAATGTTTGAACCTTACTTTTCTTGTGATAGTTAATTTTGTATGTCAACTTGACTGGACCCAGAAATTTGGTCAAACATTCTGGGTACATTTGTGAG
PBCP_001_R003,TGATGGCAAAATTCCTTTAT,AAAGAAGGGAGCAAAGAGAA,TTAAAGAGTCATGACAAAGAACTAAATATTCTTATAAGTAGCATGTTTTGAATAGTCATACCATTAAATTTTAAAAATAGGATTTACTGAGTTGTTAATTCTGTGAAAATTGGAAGTGCCTACCTTAAATTGATGGCAAAATTCCTTTATCCTGGAAAGAAGGGAGCAAAGAGAAACATCTGAGCATGTGACCAAACAGAAGAAACCCCAGTGTAGGGGTTTGAACAGTAACCCCCCTGAAAGGATATGTCCACCCAGAACCTGTGGGCATGACCTTATTTTGAAAAAGGGTCTTTGCAGATGTA
```

### Running the script

The script can be run as follows:

```
Rscript find_junction_spanning_sequences.R \
	--id=Patient_4_T4 \
	--fastq=Patient_4_T4.fq.gz \
	--flanking-sequences=flanking_sequences.csv \
	--output=Patient_4_T4.tsv \
	--max-distance=2 \
	--umi-length=10
```

Several flags have been specified, details of which are given in the table
below.

In this example, the maximum edit distance for each of the flanking
sequences is 2, i.e. up to 2 mismatches are allowed in each of the flanking
sequences. Each mismatch increases the edit distance by 1, as does an insertion
or deletion of one base. Unique molecular identifiers have been incorporated as
part of the library preparation in this example to avoid overcounting of reads
originating from the same DNA fragment/molecule. The `umi-length` parameter
allows for the first few bases at the beginning of each read to be excluded
from the search. 

Option/flag        | Description
-------------------|------------
id                 | Identifier for the dataset (useful when collating results from several FASTQ files)
fastq              | FASTQ file containing sequences to be searched for breakpoint junctions
flanking-sequences | CSV file containing flanking sequences on either side of the breakpoint junction; this is expected to contain ID, LeftFlankingSequence, RightFlankingSequence and JunctionSequence columns (default: flanking_sequences.csv)
output             | Tab-delimited output file containing sequences that match breakpoint junctions (default: junction_sequence_matches.txt)
max-distance       | Maximum edit distance for matches to each of the flanking sequences (default: 2)
umi-length         | Number of UMI bases at beginning of read to omit from search (default: 0)
