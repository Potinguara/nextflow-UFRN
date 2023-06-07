#!/usr/bin/env nextflow

params.vcfpath = "$projectDir/data/*recode.vcf"
params.indlist = "$projectDir/data/inds_freq.txt"
params.outdir = "result"


process VCFTOGENOT {
    container 'lfreitasl/vcfprocess:latest'
    publishDir(
        path: "${params.outdir}/",
        mode: 'copy'
    )
    input:
    val x
    val z

    output:
    path 'matrix_genotype_*'

    """
    vcfconverter.R $x $z
    """
}

workflow {
    matrices_ch = VCFTOGENOT(Channel.fromPath(params.vcfpath),params.indlist)
    matrices_ch.view{ it }
}
