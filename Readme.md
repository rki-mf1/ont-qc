# Nanopore Quality Control
## Scope

Quality control (QC) of ONT FASTQ data with read trimming and preliminary taxonomic classification
## Background and purpose

![ont_qc_rki_pipeline](https://github.com/user-attachments/assets/f200368a-987f-4aee-b4b7-522065e4a00b)

MF1 regularly uses Nanopore (ONT, MinION) sequencing methods for all kinds of projects within the RKI. Normally, the raw data (POD5 data) is already base called on the sequencer into FASTQ data, and also demultiplexed. Hence, this SOP only covers quality control steps including trimming and a preliminary taxonomic classification of FASTQ data of one or multiple samples, i.e. multiple read files. 

- FastQC and NanoPlot are applied to measure quality metrics before and after read trimming
- Chopper is used for quality and length trimming and filtering
- Kraken2 provides a preliminary species assignment to estimate sample purity and composition
- MultiQC collects all QC informationen of all processed samples in one QC report


# Quick start

The pipeline needs 

- samplesheet
- config file
- output directory
- Kraken2 DB (default see nextflow.config)

Download your Kraken2 DB from https://benlangmead.github.io/aws-indexes/k2
For example, Standard-8, the Standard with DB (Refeq archaea, bacteria, viral, plasmid, human1, UniVec_Core) capped at 8 GB.

## Samplesheet format

The samplesheet must be a comma separated file with headers 

```
sample,raw_read,raw_folder
```

# Install

```
[${YOUR_PROJECT}/scripts]
git clone https://github.com/rki-mf1/ont-qc.git
```

## Example call

```{bash}
[${YOUR_PROJECT}]
$ nextflow run scripts/ont-qc/main.nf \
--samplesheet analyses/samplesheet.csv \
--outdir ${YOUR_PROJECT}/analyses \
-profile rki_slurm,rki_singularity \
-c scripts/ont-qc/nextflow.config \

###
# optional parameters for the pipeline
# given examples are defaults loaded from nextflow.config
#
# optional --krakendb kraken2_20240112 \
# optional --chopper_args " -q 20 -l 500 "

###
# optional parameters to Nextflow
# optional -with-report logs/test.$(date +%T).report.html \
# optional -resume \
# optional -stub

```



# All parameters

`--samplesheet` - comma separated input file with [sample,fwd_read,rev_read,raw_folder]

`--krakendb` - define kraken database

`--outdir` - define output folder

`--chopper_args` - set params for Chopper

additional paramters transferred to Nextflow

`-c` - config file for Nextflow run

`-profile` - recommended options are rki_slurm,rki_mamba

`-stub` - run test run

`-resume` - continue process from last interrupt

`-with-report` - emit a report on run time efficiency
