#!/usr/bin/env nextflow

params.vcfpath = "$projectDir/data/*recode.vcf"
params.indlist = "$projectDir/data/inds_freq.txt"
params.outdir = "result"


process VCFTOGENOT {
    container 'lfreitasl/vcfprocess:latest'
    
    input:
    path x
    path z

    output:
    path '*.txt'

    """
    vcfconverter.R $x $z
    """
}

process HIERARCHCLUST {
    container 'lfreitasl/lalgorithms:latest'
    publishDir(
        path: "${params.outdir}",
        mode: 'copy'
    )
    input:
    path genotype

    output:
    path '*.pdf'

    """
    hierarchclust.py $genotype
    """
}

workflow {
    matrices_ch = VCFTOGENOT(Channel.fromPath(params.vcfpath),params.indlist)
    clusters_ch = HIERARCHCLUST(matrices_ch)
    clusters_ch.view{ it }
}
