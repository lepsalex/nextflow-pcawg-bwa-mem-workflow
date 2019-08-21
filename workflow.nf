#!/usr/bin/env nextflow

process getBasename {

    container "quay.io/pancancer/pcawg-bwa-mem"

    // we don't need filenames in the future (I think)
    input:
    val bamPath from params.unalignedBams

    output:
    set val(bamPath), stdout into namedBams

    """
    basename ${bamPath}
    """
}

// process readHeader {

//     container "quay.io/pancancer/pcawg-bwa-mem"

//     input:
//     set file(bamPath), bamName from namedBams

//     output:

//     """
//     samtools view -H ${bamPath} | \
//     perl -nae 'next unless /^\@RG/; s/\tPI:\t/\t/; s/\tPI:\s*\t/\t/; s/\tPI:\s*\z/\n/; s/\t/\\t/g; print' > "${bamName}_header.txt"
//     """
// }

bamNames.println { it }