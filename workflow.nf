#!/usr/bin/env nextflow

// Script parameters
params.outputFilePrefix = "aligned_"
params.threads = 1
params.sortMemMb = 512

// Inputs
bams = Channel.fromPath('input/*.bam')
bams.into { bams_rh; bams_cr }

reference_gz = file(params.reference_gz)
reference_gz_fai = file(params.reference_gz_fai)
reference_gz_amb = file(params.reference_gz_amb)
reference_gz_ann = file(params.reference_gz_ann)
reference_gz_bwt = file(params.reference_gz_bwt)
reference_gz_pac = file(params.reference_gz_pac)
reference_gz_sa = file(params.reference_gz_sa)


process readHeader {
    container "quay.io/pancancer/pcawg-bwa-mem"

    input:
    file(bam) from bams_rh

    output:
    set val("${bam.getBaseName()}"), file(bam), file("${bam.getBaseName()}_header.txt") into headers

    """
    samtools view -H ${bam} | \
    perl -nae 'next unless /^\\@RG/; s/\\tPI:\\t/\\t/; s/\\tPI:\\s*\\t/\\t/; s/\\tPI:\\s*\\z/\\n/; s/\\t/\\\t/g; print' > \
    ${bam.getBaseName()}_header.txt
    """
}

process countReads {
    container "quay.io/pancancer/pcawg-bwa-mem"

    input:
    file(bam) from bams_cr

    output:
    set val("${bam.getBaseName()}"), file("${bam.getBaseName()}_read_count.txt") into counts

    """
    samtools view ${bam} | \
    wc -l > "${bam.getBaseName()}_read_count.txt"
    """
}

// Join headers with counts
headers_with_counts = headers.join(counts)

process align {
    container "quay.io/pancancer/pcawg-bwa-mem"

    input:
    set val(bamName), file(bam), file(bamHeader), file(readCount) from headers_with_counts

    output:
    set val(bamName), file("${bamName}_aligned.bam"), file(bamHeader), file(readCount) into aligned_bams

    """
    bamtofastq exlcude=QCFAIL,SECONDARY,SUPPLEMENTARY T=${bamName}.t S=${bamName}.s O=${bamName}.o O2=${bamName}.o2 collate=1 tryoq=1 filename=${bam} | \
    bwa mem -p -t ${params.threads} -T 0 -R "${bamHeader}" ${reference_gz} - | \
    bamsort blockmb=${params.sortMemMb} inputformat=sam level=1 outputthreads=2 calmdnm=1 calmdnmrecompindetonly=1 calmdnmreference=${reference_gz} tmpfile=${bamName}.sorttmp O=${bamName}_aligned.bam
    """
}

process bam_stats_qc {
    container "quay.io/pancancer/pcawg-bwa-mem"

    publishDir "bas_out"

    input:
    set val(bamName), file(bam), file(bamHeader), file(readCount) from aligned_bams

    output:
    file "${bamName}.bas" into bam_stats

    """
    bam_stats -i ${bam} -o ${bamName}.bas
    """
}