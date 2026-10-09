#!/usr/bin/env Rscript
#### created by Diego G. Salas
#
# Stage 04 simulates paired-end fragments with polyester for every sample in the info file.
# The design fixes how many fragments each transcript gets, the sample sets read length and seed.
# The truth is then read back from the read names, so it counts what was written and not what was asked.

##### Getopts #####
#-----------------------------------------------------------
suppressPackageStartupMessages(library(optparse))
option.list = list(
  make_option("--fasta",      type = "character", default = NULL, help = "transcripts.fa from stage 02"),
  make_option("--design",     type = "character", default = NULL, help = "design.tsv from stage 03"),
  make_option("--info",       type = "character", default = NULL, help = "samples.tsv with sample name readLength seed source"),
  make_option("--out",        type = "character", default = NULL, help = "folder for one subfolder per sample"),
  make_option("--frag-mean",  type = "integer",   default = 250L),
  make_option("--frag-sd",    type = "integer",   default = 25L),
  make_option("--error-rate", type = "double",    default = 0.005),
  make_option("--log",        type = "character", default = NULL, help = "run log, TSV"),
  make_option("--stamp",      type = "character", default = "",   help = "timestamp written to the log")
)
opt = parse_args(OptionParser(option_list = option.list))

error.flag = FALSE
for (name in c("fasta", "design", "info", "out")) {
  if (is.null(opt[[name]])) { cat("/!\\ Error : --", name, " is required\n", sep = "", file = stderr()); error.flag = TRUE }
}
if (error.flag) { print_help(OptionParser(option_list = option.list)); quit(status = 1) }
for (f in c(opt$fasta, opt$design, opt$info)) {
  if (!file.exists(f)) { cat("/!\\ Error : ", f, " does not exist, run the earlier stages first\n", sep = "", file = stderr()); quit(status = 1) }
}

##### Libraries #####
#-----------------------------------------------------------
t = Sys.time()
suppressPackageStartupMessages({
  library(polyester)
  library(Biostrings)
})
cat("Libraries\tdone in ", round(difftime(Sys.time(), t, units = "secs")), "s\tpolyester ",
    as.character(packageVersion("polyester")), "\n", sep = "")

##### Functions #####
#-----------------------------------------------------------
truth.from.reads = function(reads.file) {
  read.names = names(readDNAStringSet(reads.file))
  units = sub(";.*$", "", sub("^[^/]*/", "", read.names))
  counts.df = as.data.frame(table(transcript_id = units), stringsAsFactors = FALSE)
  names(counts.df)[2] = "fragments_simulated"
  counts.df
}

##### Data files #####
#-----------------------------------------------------------
t = Sys.time()
tx.dna = readDNAStringSet(opt$fasta)
design.df = read.delim(opt$design, stringsAsFactors = FALSE)
if (!identical(names(tx.dna), design.df$transcript_id)) {
  cat("/!\\ Error : the design and the FASTA list different transcripts or in a different order\n", file = stderr())
  quit(status = 1)
}
info.df = read.delim(opt$info, stringsAsFactors = FALSE)
cat("Data files\tdone in ", round(difftime(Sys.time(), t, units = "secs")), "s\t",
    length(tx.dna), " transcripts, ", sum(design.df$fragments), " fragments per sample\n", sep = "")

##### Simulation #####
#-----------------------------------------------------------
keep = design.df$fragments > 0
keep.fasta.file = file.path(opt$out, "transcripts_with_fragments.fa")
dir.create(opt$out, showWarnings = FALSE, recursive = TRUE)
writeXStringSet(tx.dna[keep], keep.fasta.file)

for (i in seq_len(nrow(info.df))) {
  name = info.df$name[i]
  out.dir = file.path(opt$out, name)
  reads.file = file.path(out.dir, "sample_01_1.fasta.gz")
  truth.file = file.path(out.dir, "truth.tsv")
  if (file.exists(truth.file)) { cat("[04_simulate] ", truth.file, " exists, skipping ", name, "\n", sep = "", file = stderr()); next }

  cmd = sprintf("simulate_experiment_countmat(readlen = %d, seed = %d, fraglen = %d, fragsd = %d, error_rate = %g, paired = TRUE, gzip = TRUE)",
                info.df$readLength[i], info.df$seed[i], opt$`frag-mean`, opt$`frag-sd`, opt$`error-rate`)
  cat("[04_simulate] ", name, "\t", cmd, "\n", sep = "")
  if (!is.null(opt$log)) cat(name, "\t04_simulate\t", cmd, "\t", opt$stamp, "\n", sep = "", file = opt$log, append = TRUE)
  t = Sys.time()
  dir.create(out.dir, showWarnings = FALSE, recursive = TRUE)
  simulate_experiment_countmat(fasta = keep.fasta.file, readmat = matrix(design.df$fragments[keep], ncol = 1),
                               outdir = out.dir, readlen = info.df$readLength[i], seed = info.df$seed[i],
                               fraglen = opt$`frag-mean`, fragsd = opt$`frag-sd`, error_rate = opt$`error-rate`,
                               paired = TRUE, gzip = TRUE)
  truth.df = truth.from.reads(reads.file)
  write.table(truth.df, truth.file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat(name, "\t", sum(truth.df$fragments_simulated), " fragments written\tdone in ",
      round(difftime(Sys.time(), t, units = "secs")), "s\n", sep = "")
}
