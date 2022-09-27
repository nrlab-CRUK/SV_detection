#!/usr/bin/env nextflow
 
nextflow.enable.dsl = 2

process extract_soft_clipped_reads {
    memory "2 GB"
    time "24 hour"

    input:
        tuple val(id), path(bam)

    output:
        tuple val(id), path(fastq)

    script:
        fastq = "${id}.fq.gz"
        """
        set -o pipefail
        ( ${params.samtools} view -H ${bam}; ${params.samtools} view -F2048 ${bam} | gawk '( match(\$6, /^([0-9]+)S/, s) && s[1] >= ${params.min_soft_clipped} ) || ( match(\$6, /([0-9]+)S\$/, s) && s[1] >= ${params.min_soft_clipped} )' ) | ${params.samtools} fastq -N -o ${fastq}
        """
}


process find_junction_spanning_reads {
    memory "2 GB"
    time "24 hour"

    input:
        tuple val(id), path(fastq), path(flanking_sequences)

    output:
        path matches

    script:
        matches = "${id}.tsv"
        """
        ${params.rscript} ${projectDir}/find_junction_spanning_sequences.R --id=${id} --fastq=${fastq} --flanking-sequences=${flanking_sequences} --output=${matches} --max-distance=${params.max_distance} --umi-length=${params.umi_length}
        """
}


workflow {

    flanking_sequences = Channel.fromPath(params.flanking_sequences, checkIfExists: true)

    sample_sheet = Channel.fromPath(params.sample_sheet, checkIfExists: true)

    bam = sample_sheet
        .splitCsv(header: true)
        .map { row -> tuple(row.ID, file(row.BAM, checkIfExists: true)) }

    bam \
      | extract_soft_clipped_reads \
      | splitFastq(by: params.chunk_size, file: true, compress: true) \
      | combine(flanking_sequences) \
      | find_junction_spanning_reads \
      | collectFile(storeDir: params.results_dir, keepHeader: true)
}

