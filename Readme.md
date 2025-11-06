# Nanopore Quality Control
## Scope

Quality control (QC) of ONT FASTQ data with read trimming and preliminary taxonomic classification

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

Retrieve code via 

```
[${YOUR_PROJECT}/scripts]
git clone https://github.com/rki-mf1/ont-qc.git
```

## Example call

```{bash}
[${YOUR_PROJECT}]
$ nextflow run scripts/ont-qc/main.nf \
--samplesheet samplesheet.csv \
--outdir ${YOUR_PROJECT}/analyses \
-profile slurm,singularity \
-c scripts/ont-qc/nextflow.config \

###
# optional parameters for the pipeline
# given examples are defaults loaded from nextflow.config
#
# optional --krakendb /path/to/kraken2_DB \
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

additional parameters transferred to Nextflow

`-c` - config file for Nextflow run

`-profile` - choose e.g. singularity vs. mamba (if you use mama, add the --with-conda statement to the call)

`-stub` - run test run

`-resume` - continue process from last interrupt


`-with-report` - emit a report on run time efficiency
