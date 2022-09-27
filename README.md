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

### Output

The output file is a tab-delimited file (TSV) containing the columns described
in the following table. Each row corresponds to a match for one of the FASTQ
sequence records to one of the junctions specified in the flanking sequences CSV
file.

Column                   | Description
-------------------------|---------------------------
ID                       | The dataset identifier as specified using the `id` flag
SequenceID               | The sequence name or ID from a FASTQ record that matches one of the SV junctions
UMI                      | The UMI tag excluded from the search if specified by the `umi-length` flag
Sequence                 | The sequence from the FASTQ record (excluding the trimmed UMI tag if present)
SVID                     | The matching SV identifier (the ID from the flanking sequences CSV file)
LeftFlankingSequence     | The left flanking sequence for the SV that the read sequence matches
RightFlankingSequence    | The right flanking sequence for the SV that the sequence record matches
Direction                | The direction of the match, either 'Forward' or 'Reverse'
MatchingPosition         | The position of the matching segment within the sequence
MatchingSequence         | The matching segment within the sequence that contains the left and right flanking sequences
LeftDistance             | The edit distance for the matching left flanking sequence
RightDistance            | The edit distance for the matching right flanking sequence
JunctionPosition         | The position of the matching junction sequence or NA if there is no match or the junction sequence is not specified in the flanking sequences CSV file
JunctionMatchingSequence | The portion of the sequence that matches that junction sequence if specified
JunctionDistance         | The edit distance for the matching junction sequence
JunctionOverlap          | The minimum overlap of the matching sequence with the junction sequence assuming that the mid-point of the junction sequence is the actual junction location


## Running the `junction_detection.nf` workflow on WGS data using Nextflow

The `junction_detection.nf` workflow can be used to run the
`find_junction_spanning_sequences.R` script on very large sequencing datasets
such as would be obtained from deep whole genome sequencing (WGS). Unlike
amplicon sequencing in which most reads are expected to match a rearrangement
junction since they result from PCR amplification of the genomic rearrangement,
only a very small fraction of reads in a WGS dataset will match any of the
expected SV junctions. These reads are likely to have split alignments when
aligned to the reference genome with an aligner such as
[bwa-mem](https://github.com/lh3/bwa).

Split alignments are represented in
the [BAM](https://en.wikipedia.org/wiki/Binary_Alignment_Map) files produced by
the aligner as a primary alignment with soft-clipping of the part of the read
that maps to the other side of the junction and a supplementary alignmnet for
that other part. The workflow extracts all soft-clipped reads which have a
configurable minimum number of bases that have been clipped. It then splits the
resulting set of sequences into chunks so that the computational work of running
the `find_junction_spanning_sequences.R` script can be distributed over multiple
processors, e.g. it can be parallelized by running on a multi-core server or a
high-performance compute cluster.

### Prerequisites

In addition to the R package dependencies described above, the junction
detection pipeline requires Nextflow and a Java 11 runtime to be installed.
Assuming that Java is already installed, installing Nextflow is very
straightforward with the following command:

```
curl -s https://get.nextflow.io | bash 
```

This will create a file called `nextflow`. This is an executable file that is
used to run Nextflow pipelines. The `nextflow` file can be moved to the home
directory or a `bin` subdirectory or another directory that is added to the
PATH.

The pipeline also requires [`samtools`](http://www.htslib.org) which is used
to extract soft-clipped reads from the input BAM files.

### Inputs

The input files for the workflow are BAM files resulting from aligning the
sequence data against a reference genome with an aligner, such as bwa, capable
of producing clipped alignments.

In addition to the flanking sequences CSV file required by the
`find_junction_spanning_sequences.R` script and described above, the workflow is
configured using a configuration file, an example of which is provided in this
repository (see [`nextflow.config`](nextflow.config)).

Also required is a sample sheet CSV file, named `sample_sheet.csv` by default,
that contains ID and BAM columns, where the ID is the dataset ID provided to the
R script using the `id` flag (see the table listing the options given above) and
the BAM column contains the name or path of the BAM file for that dataset.

### Configuring the workflow

Create a file called `junction_detection.config` in the run directory and make changes to
parameter settings as required. An example is given below:

```
// junction_detection.config
params {
    sample_sheet       = "sample_sheet.csv"
    flanking_sequences = "flanking_sequences.csv"
    results_dir        = "results"
    min_soft_clipped   = 20
    chunk_size         = 100000
    max_distance       = 2
    umi_length         = 10
    samtools           = "samtools"
    rscript            = "Rscript"
}
```

The parameters are essentially the same as those described above for the R
script. In addition it is possible to specify the chunk size, i.e. the number of
reads within each chunk to be processed as a separate job, and also paths to the
`Rscript` and `samtools` executables used by the pipeline (note that by default
those available on the PATH are used).

Run the workflow as follows:

```
nextflow run nrlab-CRUK/SV_detection -config junction_detection.config
```

Nextflow will download the workflow from the GitHub repository automatically the
first time it is run.

The pipeline is configured with 3 different run profiles. The standard profile
used by default will use up to 4 CPU cores and 8GB memory. Also available is a
profile called `bigserver` that can use up to 40 CPU cores and 128GB memory;
this can be selected using the `-profile` option.

```
nextflow run nrlab-CRUK/SV_detection -config junction_detection.config -profile bigserver
```

The third profile can be used for submitting jobs to a compute cluster with the
SLURM scheduler. More details on configuring and specifying run profiles can be
found in the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html).

### Outputs

The output files are a single TSV file for each dataset with the same columns as
described above for the `find_junction_spanning_sequences.R` script.
