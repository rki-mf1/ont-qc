# Quick start

The pipeline needs 

- samplesheet
- config file
- output directory
- Kraken2 DB (default see nextflow.config)

## Samplesheet format

The samplesheet must be a comma separated file with headers 

```
sample,raw_read,raw_folder
```

# Install

Retrieve code via 

```
hpc-login02[/scratch/projekte/MF1_BI-Support/${YOUR_PROJECT}/scripts]
git clone https://bitbucket.rki.local/scm/bbap/i008_trappek_ontqc.git
```

## Example call

```{bash}
hpc-login02[/scratch/projekte/MF1_BI-Support/${YOUR_PROJECT}]
$ nextflow run scripts/i008_trappek_ontqc/main.nf \
--samplesheet analyses/samplesheet.csv \
--outdir /scratch/projekte/MF1_BI-Support/${YOUR_PROJECT}/analyses \
--with-conda \
-profile rki_slurm,rki_mamba \
-c scripts/i008_trappek_ontqc/nextflow.config \

###
# optional parameters for the pipeline
# given examples are defaults loaded from nextflow.config
#
# optional --krakendb /scratch/databases/kraken2_20240112/ \
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