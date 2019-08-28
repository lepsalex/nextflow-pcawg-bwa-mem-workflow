#!/usr/bin/env nextflow

// Script parameter defaults
params.outputDir = "./out"
params.outputFilePrefix = "example"
params.inputDir = "./input"
params.threads = 1
params.sortMemMb = 1024

// Inputs

// TODO: should check for existence first via ifEmpty() - see any nf-core example
// ----- Also applicable to the file() syntax below - checkIfExists: true
// ----- https://github.com/nf-core/chipseq/blob/master/main.nf
Channel.fromPath("${params.inputDir}/*.bam").into { bams_rh; bams_cr }

// Nextflow docs are a little conflicting, they say that a single file channel can
// be reused multiple times but when trying to do it results in an error so for now
// going to use direct file declaration like this, based on example on nextflow.oi:
// https://www.nextflow.io/example4.html
reference_gz_ch = file(params.reference_gz)
reference_gz_fai_ch = file(params.reference_gz_fai)
reference_gz_amb_ch = file(params.reference_gz_amb)
reference_gz_ann_ch = file(params.reference_gz_ann)
reference_gz_bwt_ch = file(params.reference_gz_bwt)
reference_gz_pac_ch = file(params.reference_gz_pac)
reference_gz_sa_ch = file(params.reference_gz_sa)


process readHeader {

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

    input:
    file reference_gz from reference_gz_ch
    // While not explicitly used, files are implicitly used by
    // convention in bwa mem (ie. explicit file_name.fa.gz => file_name.fa.gz.fai) 
    file reference_gz_fai from reference_gz_fai_ch
    file reference_gz_amb from reference_gz_amb_ch
    file reference_gz_ann from reference_gz_ann_ch
    file reference_gz_bwt from reference_gz_bwt_ch
    file reference_gz_pac from reference_gz_pac_ch
    file reference_gz_sa from reference_gz_sa_ch
    // Map headers to extract header text then join with counts
    set val(bamName), file(bam), file(bamHeader), val(headerText), file(readCount) from headers.map {
        [it[0], it[1], it[2], file(it[2]).text]
    }.join(counts)

    output:
    set val(bamName), file("${bamName}_aligned.bam"), file(bamHeader), file(readCount) into aligned_bams

    """
    bamtofastq exclude=QCFAIL,SECONDARY,SUPPLEMENTARY T=${bamName}.t S=${bamName}.s O=${bamName}.o O2=${bamName}.o2 collate=1 tryoq=1 filename=${bam} | \\
    bwa mem -p -t ${params.threads} -T 0 -R "${headerText.trim()}" ${reference_gz} - | \\
    bamsort blockmb=${params.sortMemMb} inputformat=sam level=1 outputthreads=2 calmdnm=1 calmdnmrecompindetonly=1 calmdnmreference=${reference_gz} tmpfile=${bamName}.sorttmp O=${bamName}_aligned.bam
    """
}

process bam_stats_qc {

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

process merge_aligned {

    publishDir "${params.outputDir}", mode: 'copy', overwrite: true

    input:
    // Collecting in channel, good thread on exactly this use-case:
    // https://groups.google.com/d/msg/nextflow/jt77_-uApMs/2X_74ireBQAJ
    file("verified_bam_") from verified_bams.collect()

    output:
    file "${params.outputFilePrefix}_aligned.bam" into mb_for_extract_ur, mb_for_extract_bru
    file "${params.outputFilePrefix}_aligned.metrics" into merged_bam_metrics   

    """
    bammarkduplicates \
    I=${verified_bam_.join(" I=")} \
    O=${params.outputFilePrefix}_aligned.bam \
    M=${params.outputFilePrefix}_aligned.metrics \
    tmpfile=${params.outputFilePrefix}_aligned.biormdup \
    markthreads=${params.threads} \
    rewritebam=1 \
    rewritebamlevel=1 \
    index=1 \
    md5=1
    """
}

process extract_unaligned_reads {
    // will create one file for each f value
    extract_flags = [4, 8]

    input:
    file reference_gz from reference_gz_ch
    file reference_gz_fai from reference_gz_fai_ch
    file merged_bam from mb_for_extract_ur
    each f from extract_flags

    output:
    file "${params.outputFilePrefix}_unmappedReads_f${f}.bam" into extract_ur_unmapped

    // samtools: -f option = Only output alignments with all bits set in INT present in the FLAG field.
    """
    samtools view -h -f ${f} ${merged_bam} | \\
    remove_both_ends_unmapped_reads.pl | \\
    bamsort blockmb=${params.sortMemMb} inputformat=sam level=1 outputthreads=2 calmdnm=1 calmdnmrecompindetonly=1 calmdnmreference=${reference_gz} tmpfile=${params.outputFilePrefix}.sorttmp O=${params.outputFilePrefix}_unmappedReads_f${f}.bam
    """
}

process extract_both_reads_unaligned {

    input:
    file merged_bam from mb_for_extract_bru

    output:
    file "${params.outputFilePrefix}_unmappedReads_f12.bam" into extract_bru_unmapped

    """
    samtools view -h -b -f 12 ${merged_bam} > "${params.outputFilePrefix}_unmappedReads_f12.bam"
    """
}

process merge_unmappedReads {

    publishDir "${params.outputDir}", mode: 'copy', overwrite: true
    
    input:
    // Collecting in channel, good thread on exactly this use-case:
    // https://groups.google.com/d/msg/nextflow/jt77_-uApMs/2X_74ireBQAJ
    file("unmapped_") from extract_ur_unmapped.concat(extract_bru_unmapped).collect()

    output:
    file "${params.outputFilePrefix}_unmappedRead.bam"
    file "${params.outputFilePrefix}_unmappedRead.metrics" 

    """
    bammarkduplicates \
    I=${unmapped_.join(" I=")} \
    O=${params.outputFilePrefix}_unmappedRead.bam \
    M=${params.outputFilePrefix}_unmappedRead.metrics \
    tmpfile=${params.outputFilePrefix}_unmappedRead.biormdup \
    markthreads=${params.threads} \
    rewritebam=1 \
    rewritebamlevel=1 \
    index=1 \
    md5=1
    """
}