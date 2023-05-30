#!/usr/bin/env nextflow

params.greeting = 'Hello world! This is a bit longer test, because I wanna see if tower is really that funny as everyone is saying!'
greeting_ch = Channel.of(params.greeting)

process SPLITLETTERS {
    input:
    val x

    output:
    path 'letras_*'

    """
    printf '$x' | split -b 2 - letras_
    """
}

process CONVERTTOUPPER {
    input:
    path y

    output:
    stdout

    """
    cat $y | tr '[a-z]' '[A-Z]' 
    """
}

workflow {
    letters_ch = SPLITLETTERS(greeting_ch)
    results_ch = CONVERTTOUPPER(letters_ch.flatten())
    results_ch.view{ it }
}

