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
    set file(bam), val("${bam.getBaseName()}"), file("${bam.getBaseName()}_header.txt") into headers

    script:
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

    script:
    """
    samtools view ${bam} | \
    wc -l > "${bam.getBaseName()}_read_count.txt"
    """
}

process align {
    container "quay.io/pancancer/pcawg-bwa-mem"

    input:
    set file(bam), val(bamName), file(bamHeader) from headers

    output:
    file "${bamName}_aligned.bam" into aligned_bams

    script:
    """
    bamtofastq exlcude=QCFAIL,SECONDARY,SUPPLEMENTARY T=${bamName + ".t"} S=${bamName + ".s"} O=${bamName + ".o"} O2=${bamName + ".o2"} collate=1 tryoq=1 filename=${bam} | \
    bwa mem -p -t ${params.threads} -T 0 -R "${bamHeader}" ${reference_gz} - | \
    bamsort blockmb=${params.sortMemMb} inputformat=sam level=1 outputthreads=2 calmdnm=1 calmdnmrecompindetonly=1 calmdnmreference=${reference_gz} tmpfile=${bamName + ".sorttmp"} O=${bamName + "_aligned.bam"}
    """
}