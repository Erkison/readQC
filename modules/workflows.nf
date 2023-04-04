include { FASTQC; MULTIQC } from '../modules/processes.nf' 

workflow QC {
    take:
        reads_ch

    main: 
        FASTQC(reads_ch)

        MULTIQC(FASTQC.out.fastqc_directories.collect())

    emit:
        FASTQC.out.fastqc_directories.collect()
}