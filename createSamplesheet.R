## CALL: Rscript createSamplesheet.R --folder /path/to/fastq_data --outfile samplesheet_pipeline_benchmark01.csv
packages <- c("optparse")

## INSTALL PACKAGES ##
# Check if the directory './Rlibs' exists, if not, create it
if (!dir.exists("Rlibs")) {
  dir.create("Rlibs", showWarnings = FALSE, recursive = TRUE)
}

# Add './Rlibs' to the library path
.libPaths(c("Rlibs", .libPaths()))

# Check which packages are not installed
installed_packages <- packages %in% rownames(installed.packages())

# Install missing packages into './Rlibs'
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages], lib="Rlibs")
} else {
  cat("All packages are already installed.\n")
}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

## ARGUMENT PARSER ##
# Parse arguments
option_list <- list(
  make_option(c("-f", "--folder"), type="character", default=NULL,
              help="directory path to the raw data", metavar="character"),
  make_option(c("-o", "--outfile"), type="character", default=NULL,
              help="name of the output file [example: samplesheet_pipeline_benchmark01.csv]", metavar="character")
)


# Create an OptionParser object
opt_parser <- OptionParser(option_list=option_list)

# Parse the command line arguments
opt <- parse_args(opt_parser)

# Print the parsed arguments
cat("Folder: ", opt$folder, "\n")
cat("Outfile: ", opt$outfile, "\n")

# Hand over 
folder <- opt$folder
header <- c("sample","raw_read","raw_folder")


## CREATE SAMPLESHEET ##
# Extract read name (with ending)
raw_read=unlist(list.files(folder, pattern=".fastq.gz"),use.names = F)

# Extract sample name (without ending)
samples <- gsub("^(.+)\\.fastq.gz","\\1", raw_read)

# Save in data frame
df <- data.frame(samples,raw_read,folder)

# Set column names
colnames(df) <- header
outfile <- opt$outfile

# Write to output
write.table(df, outfile, row.names = F, col.names = T, sep = ",", quote = F)

