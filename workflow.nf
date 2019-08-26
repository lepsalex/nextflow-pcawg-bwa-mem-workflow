#!/usr/bin/env nextflow

// Script parameters
params.outputFilePrefix = "processed"
params.threads = 1
params.sortMemMb = 1024

// Inputs (TODO: should check for vals first via ifEmpty() - see nf-core examples)
Channel.fromPath('input/*.bam').into { bams_rh; bams_cr }

Channel.fromPath(params.reference_gz).into { reference_gz_align; reference_gz_reads }
Channel.fromPath(params.reference_gz_fai).set { reference_gz_fai_ch }
Channel.fromPath(params.reference_gz_amb).set { reference_gz_amb_ch }
Channel.fromPath(params.reference_gz_ann).set { reference_gz_ann_ch }
Channel.fromPath(params.reference_gz_bwt).set { reference_gz_bwt_ch }
Channel.fromPath(params.reference_gz_pac).set { reference_gz_pac_ch }
Channel.fromPath(params.reference_gz_sa).set { reference_gz_sa_ch }


process readHeader {
    container "quay.io/pancancer/pcawg-bwa-mem"

    input:
    file bam from bams_rh

    output:
    set val("${bam.getBaseName()}"), file(bam), file("${bam.getBaseName()}_header.txt") into headers

    """
    samtools view -H ${bam} | \\
    perl -nae 'next unless /^\\@RG/; s/\\tPI:\\t/\\t/; s/\\tPI:\\s*\\t/\\t/; s/\\tPI:\\s*\\z/\\n/; s/\\t/\\\\t/g; print' > \\
    ${bam.getBaseName()}_header.txt
    """
}

process countReads {
    container "quay.io/pancancer/pcawg-bwa-mem"

    input:
    file bam from bams_cr

    output:
    set val("${bam.getBaseName()}"), file("${bam.getBaseName()}_read_count.txt") into counts

    """
    samtools view ${bam} | \\
    wc -l > ${bam.getBaseName()}_read_count.txt
    """
}

process align {
    container "quay.io/pancancer/pcawg-bwa-mem"

    input:
    val reference_gz from reference_gz_align
    // Map headers to extract header text then join with counts
    set val(bamName), file(bam), file(bamHeader), val(headerText), file(readCount) from headers.map {
        [it[0], it[1], it[2], file(it[2]).text]
    }.join(counts)

    output:
    set val(bamName), file("${bamName}_aligned.bam"), file(bamHeader), file(readCount) into aligned_bams

    """
    bamtofastq exlcude=QCFAIL,SECONDARY,SUPPLEMENTARY T=${bamName}.t S=${bamName}.s O=${bamName}.o O2=${bamName}.o2 collate=1 tryoq=1 filename=${bam} | \\
    bwa mem -p -t ${params.threads} -T 0 -R "${headerText.trim()}" $reference_gz - | \\
    bamsort blockmb=${params.sortMemMb} inputformat=sam level=1 outputthreads=2 calmdnm=1 calmdnmrecompindetonly=1 calmdnmreference=$reference_gz tmpfile=${bamName}.sorttmp O=${bamName}_aligned.bam
    """
}

process bam_stats_qc {
    container "quay.io/pancancer/pcawg-bwa-mem"

    input:
    set val(bamName), file(bam), file(bamHeader), file(readCount) from aligned_bams

    output:
    file "${bamName}.bas" into bam_stats
    file bam into verified_bams

    """
    bam_stats -i ${bam} -o ${bamName}.bas \
    && \
    verify_read_groups.pl --header-file ${bamHeader} --bas-file ${bamName + ".bas"} --input-read-count-file ${readCount}
    """
}

process merge {
    container "quay.io/pancancer/pcawg-bwa-mem"

    input:
    // Collecting in channel, good thread on exactly this use-case:
    // https://groups.google.com/d/msg/nextflow/jt77_-uApMs/2X_74ireBQAJ
    file("verified_bam_") from verified_bams.collect()

    output:
    file "${params.outputFilePrefix}.bam" into merged_bams
    file "${params.outputFilePrefix}.metrics" into merged_bam_metrics   

    """
    bammarkduplicates \
    I=${verified_bam_.join(" I=")} \
    O=${params.outputFilePrefix}.bam \
    M=${params.outputFilePrefix}.metrics \
    tmpfile=${params.outputFilePrefix}.biormdup \
    markthreads=${params.threads} \
    rewritebam=1 \
    rewritebamlevel=1 \
    index=1 \
    md5=1
    """
}

process extract_unaligned_reads {
    container "quay.io/pancancer/pcawg-bwa-mem"

    input:
    val reference_gz from reference_gz_reads
    // combine into two channels to be processed with different f values
    set file(merged_bam), val(f) from merged_bams.combine(Channel.from(4, 8))

    output:
    file "${params.outputFilePrefix}_unmappedReads_f${f}.bam" into unmapped

    // samtools: -f option = Only output alignments with all bits set in INT present in the FLAG field.
    """
    samtools view -h -f ${f} ${merged_bam} | \\
    remove_both_ends_unmapped_reads.pl | \\
    bamsort blockmb=${params.sortMemMb} inputformat=sam level=1 outputthreads=2 calmdnm=1 calmdnmrecompindetonly=1 calmdnmreference=$reference_gz tmpfile=${params.outputFilePrefix}.sorttmp O=${params.outputFilePrefix}_unmappedReads_f${f}.bam
    """
}
