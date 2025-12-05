/************************** 
* CALL
**************************/

// (nextflow_env) $ 
// nextflow run main.nf \
// --samplesheet samplesheet.txt
// --krakendb kraken2_20240112/ \

/************************** 
* WORKFLOW
**************************/

// read in csv file (sep=',' make sure to fix, if you used MS Excel)
// sample,raw_read,raw_folder

include { ONTQC } from './subworkflows/local/ontqc'

// define workflow
workflow {
    ch_infiles = Channel
                .fromPath(params.samplesheet)
                .splitCsv(header: true, sep: ',')
                .map {row -> tuple(
                                    [id: row.sample, single_end:true], //meta map
                                    file(row.raw_folder + "/" + row.raw_read)
                                    )}


    ONTQC(ch_infiles, params.krakendb, params.outdir, params.multiqc_config)

}
