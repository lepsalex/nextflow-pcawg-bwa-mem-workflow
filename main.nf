#!/usr/bin/env nextflow

// Script parameter defaults
params.outputFilePrefix = "example"

// normal processes
params.cpus = 1
params.threads = 2
params.mem = 1024

// bigger processes (align)
params.cpus_mid = 2
params.threads_mid = 4
params.mem_mid = 2048

// largest processes (align)
params.cpus_big = 4
params.threads_big = 8
params.mem_big = 4096


// Get Inputs from Minio or other S3 compatible bucket
process fetchFiles {
    container 'minio-mc-bash:latest'

    cpus params.cpus
    memory "${params.mem} MB"

    output:
    // TODO: should check for existence first via ifEmpty() - see any nf-core example
    // ----- Also applicable to the file() syntax below - checkIfExists: true
    // ----- https://github.com/nf-core/chipseq/blob/master/main.nf
    file "${params.inputsBucket}/*.bam" into bams_rh, bams_cr

    // // Nextflow docs are a little conflicting, they say that a single file channel can
    // // be reused multiple times but when trying to do it results in an error so for now
    // // going to use direct file declaration like this, based on example on nextflow.oi:
    // // https://www.nextflow.io/example4.html
    file "${params.referenceBucket}/${params.reference_gz}" into reference_gz_ch
    file "${params.referenceBucket}/${params.reference_gz_fai}" into reference_gz_fai_ch
    file "${params.referenceBucket}/${params.reference_gz_amb}" into reference_gz_amb_ch
    file "${params.referenceBucket}/${params.reference_gz_ann}" into reference_gz_ann_ch
    file "${params.referenceBucket}/${params.reference_gz_bwt}" into reference_gz_bwt_ch
    file "${params.referenceBucket}/${params.reference_gz_pac}" into reference_gz_pac_ch
    file "${params.referenceBucket}/${params.reference_gz_sa}" into reference_gz_sa_ch

    """
    mc config host add store ${params.minioURI} ${params.minioAccessKey} ${params.minioSecretKey}
    mc cp -r store/${params.referenceBucket} .
    mc cp -r store/${params.inputsBucket} .
    """
}

process readHeader {

    cpus params.cpus
    memory "${params.mem} MB"

    input:
    file(bam) from bams_rh.flatMap()

    output:
    set val("${bam.baseName}"), file(bam), file("${bam.baseName}_header.txt") into headers

    """
    samtools view -H ${bam} | \\
    perl -nae 'next unless /^\\@RG/; s/\\tPI:\\t/\\t/; s/\\tPI:\\s*\\t/\\t/; s/\\tPI:\\s*\\z/\\n/; s/\\t/\\\\t/g; print' > \\
    ${bam.baseName}_header.txt
    """
}

process countReads {

    cpus params.cpus
    memory "${params.mem} MB"

    input:
    file bam from bams_cr.flatMap()

    output:
    set val("${bam.baseName}"), file("${bam.baseName}_read_count.txt") into counts

    """
    samtools view ${bam} | \\
    wc -l > ${bam.baseName}_read_count.txt
    """
}

process align {

    cpus params.cpus_big
    memory "${params.mem_big} MB"

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
    bamsort blockmb=${params.mem} inputformat=sam level=1 outputthreads=2 calmdnm=1 calmdnmrecompindetonly=1 calmdnmreference=${reference_gz} tmpfile=${bamName}.sorttmp O=${bamName}_aligned.bam
    """
}

process bam_stats_qc {

    cpus params.cpus
    memory "${params.mem} MB"

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

    cpus params.cpus_mid
    memory "${params.mem_mid} MB"

    input:
    // Collecting in channel, good thread on exactly this use-case:
    // https://groups.google.com/d/msg/nextflow/jt77_-uApMs/2X_74ireBQAJ
    file("verified_bam_") from verified_bams.collect()

    output:
    file "${params.outputFilePrefix}_aligned.bam" into mb_for_extract_ur, mb_for_extract_bru, mappedReadsOutput
    file "${params.outputFilePrefix}_aligned.metrics" into merged_bam_metrics, mappedReadsMetricsOutput

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
    cpus params.cpus
    memory "${params.mem} MB"

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
    bamsort blockmb=${params.mem} inputformat=sam level=1 outputthreads=2 calmdnm=1 calmdnmrecompindetonly=1 calmdnmreference=${reference_gz} tmpfile=${params.outputFilePrefix}.sorttmp O=${params.outputFilePrefix}_unmappedReads_f${f}.bam
    """
}

process extract_both_reads_unaligned {
    cpus params.cpus
    memory "${params.mem} MB"

    input:
    file merged_bam from mb_for_extract_bru

    output:
    file "${params.outputFilePrefix}_unmappedReads_f12.bam" into extract_bru_unmapped

    """
    samtools view -h -b -f 12 ${merged_bam} > "${params.outputFilePrefix}_unmappedReads_f12.bam"
    """
}

process merge_unmappedReads {

    cpus params.cpus_mid
    memory "${params.mem_mid} MB"
    
    input:
    // Collecting in channel, good thread on exactly this use-case:
    // https://groups.google.com/d/msg/nextflow/jt77_-uApMs/2X_74ireBQAJ
    file("unmapped_") from extract_ur_unmapped.concat(extract_bru_unmapped).collect()

    output:
    file "${params.outputFilePrefix}_unmappedRead.bam" into unmappedReadsOutput
    file "${params.outputFilePrefix}_unmappedRead.metrics" into unmappedReadsMetricsOutput

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

process outputMappedResults {
    
    cpus params.cpus
    memory "${params.mem} MB"

    container 'minio-mc-bash:latest'

    input:
    file output from mappedReadsOutput

    """
    mc config host add store ${params.minioURI} ${params.minioAccessKey} ${params.minioSecretKey}
    mc cp ${output} store/${params.outputBucket}
    """
}

process outputMappedMetricsResults {
    
    cpus params.cpus
    memory "${params.mem} MB"

    container 'minio-mc-bash:latest'

    input:
    file output from mappedReadsMetricsOutput

    """
    mc config host add store ${params.minioURI} ${params.minioAccessKey} ${params.minioSecretKey}
    mc cp ${output} store/${params.outputBucket}
    """
}

process outputUnmappedResults {
    
    cpus params.cpus
    memory "${params.mem} MB"

    container 'minio-mc-bash:latest'

    input:
    file output from unmappedReadsOutput

    """
    mc config host add store ${params.minioURI} ${params.minioAccessKey} ${params.minioSecretKey}
    mc cp ${output} store/${params.outputBucket}
    """
}

process outputUnmappedMetricsResults {
    
    cpus params.cpus
    memory "${params.mem} MB"

    container 'minio-mc-bash:latest'

    input:
    file output from unmappedReadsMetricsOutput

    """
    mc config host add store ${params.minioURI} ${params.minioAccessKey} ${params.minioSecretKey}
    mc cp ${output} store/${params.outputBucket}
    """
}
