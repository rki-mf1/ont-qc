/************************** 
* DEFINE INCLUDES
**************************/
// Taxonomic classification
include { KRAKEN2_KRAKEN2 } from '../../modules/nf-core/kraken2/kraken2/main'
include { KRAKEN2_KRAKEN2 as KRAKEN2_CHOPPED } from '../../modules/nf-core/kraken2/kraken2/main'
include { KRAKENTOOLS_KREPORT2KRONA } from '../../modules/nf-core/krakentools/kreport2krona/main'
include { KRONA_KTIMPORTTEXT } from '../../modules/nf-core/krona/ktimporttext/main'
//Read QC
include { CHOPPER } from '../../modules/nf-core/chopper/main'
include { FASTQC } from '../../modules/nf-core/fastqc/main'
include { FASTQC as FASTQC_CHOPPED } from '../../modules/nf-core/fastqc/main'
include { NANOPLOT } from '../../modules/nf-core/nanoplot/main'
include { NANOPLOT as NANOPLOT_CHOPPED} from '../../modules/nf-core/nanoplot/main'

include { MULTIQC } from '../../modules/nf-core/multiqc/main'


//
// Get software versions for pipeline
//
def myProcessVersionsFromYAML(yaml_file) {
    def yaml = new org.yaml.snakeyaml.Yaml()
    def loadedYaml = yaml.load(yaml_file)
    // Ensure loadedYaml is a Map and flatten nested maps
    def versions = (loadedYaml instanceof Map ? loadedYaml : [:])
        .collectEntries { k, v ->
            if (v instanceof Map) {
                // If value is a Map, extract its key-value pairs
                v.collectEntries { nestedKey, nestedValue ->
                    [nestedKey, nestedValue instanceof String ? nestedValue.trim() : nestedValue]
                }
            } else if (v instanceof String && v.trim()) {
                // For non-nested strings, keep them if non-empty
                [k.tokenize(':')[-1], v.trim()]
            } else {
                // Skip other cases
                [:]
            }
        }

    // use indentation for correct formatting
    List<String> indented_versions = []
    for (String tool : versions.keySet()){
        indented_tool = "  " + tool
        indented_versions.add(indented_tool + ": \"" + versions.get(tool)+"\"")
    }
    return indented_versions.join("\n")
}

workflow ONTQC {
    take:
        ch_infiles
        ch_krakendb
        ch_outdir
        ch_multiqc_config

    main:
    ch_versions = Channel.empty()

    // read trimming chopper
    // INPUT: raw reads
    // OUTPUT: chopped reads
    ch_trimmed = CHOPPER(ch_infiles)
    ch_versions = ch_versions.mix(ch_trimmed.versions.first())

    // read classification
    // INPUT: reads
    // OUTPUT: classification report
    // CALL: KRAKEN2_KRAKEN2(tuple val(meta), path(reads), krakendb, save_output_fastqs, save_reads_assignment)
    ch_kraken = KRAKEN2_KRAKEN2(ch_infiles, ch_krakendb, false, false)
    ch_kraken_chopped = KRAKEN2_CHOPPED(ch_trimmed.fastq, ch_krakendb, false, false)
    ch_versions = ch_versions.mix(ch_kraken.versions.first())

    // Krona plots of classification
    // INPUT: Kraken report
    // OUTPUT: Krona report
    ch_krona_rprt = KRAKENTOOLS_KREPORT2KRONA(ch_kraken.report)
    ch_versions = ch_versions.mix(ch_krona_rprt.versions.first())
    // INPUT: Krona report
    // OUTPUT: HTML report
    ch_krona = KRONA_KTIMPORTTEXT(ch_krona_rprt.txt)
    ch_versions = ch_versions.mix(ch_krona.versions.first())

    // read QC fastqc
    // INPUT: reads
    // OUTPUT: QC report
    ch_fastqc = FASTQC(ch_infiles)
    ch_fastqc_chopped = FASTQC_CHOPPED(ch_trimmed.fastq)
    ch_versions = ch_versions.mix(ch_fastqc.versions.first())

    // read QC NanoPlot
    // INPUT: raw and chopped reads
    // OUTPUT: QC report
    ch_nanoplot_raw = NANOPLOT(ch_infiles)
    ch_nanoplot_chopped = NANOPLOT_CHOPPED(ch_trimmed.fastq)

    ch_collated_versions = Channel.value('software_versions:')
                                .concat(ch_versions.unique().map { version -> myProcessVersionsFromYAML(version) })
                                .collectFile(name: 'software_versions.yaml', storeDir: "${ch_outdir}/pipeline_info", newLine: true, sort: 'index')

    
    // multiQC
    // INPUT: QC and classification reports
    // MULTIQC(ch_multiqc_files, multiqc_config_yaml, extra_multiqc_config, multiqc_logo)
    // OUTPUT: multiQC report
    // collecting input
    ch_multiqc_files = Channel.empty()
    // FastQC raw and chopped
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_fastqc.zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_fastqc_chopped.zip.collect{it[1]}.ifEmpty([]))
    // Kraken
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_kraken.report.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_kraken_chopped.report.collect{it[1]}.ifEmpty([]))
    // Nanoplot raw and chopped
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_nanoplot_raw.txt.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_nanoplot_chopped.txt.collect{it[1]}.ifEmpty([]))

    // mix versions into channel
    // https://training.nextflow.io/hello_nextflow/09_hello_nf-core/#adding-parameters-to-your-pipeline_1
    // https://github.com/nextflow-io/nextflow/discussions/4925
    //ch_multiqc_files = ch_multiqc_files.mix(
     //   ch_collated_versions)

    ch_multiqc = MULTIQC(ch_multiqc_files.collect(), ch_multiqc_config, ch_collated_versions, [], [], [])

}