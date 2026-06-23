#!/usr/bin/env Rscript

# ============================================================
#  EucsPhylogenomics
#
#  Usage: Rscript run_pipeline.R --config config.yaml
#         Rscript run_pipeline.R --config config.yaml --redo
# ============================================================

# --- Load libraries and function ----------------------------
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(yaml))

# create a function to retrieve parameter value with default
f_get_param <- function(value, default) {
  if (is.null(value) || identical(value, "")) {
    default
  } else {
    value
  }
}

# --- Argument parsing ----------------------------------------
option_list <- list(
  make_option(c("-c", "--config"), type="character", default=NULL,
              help="Path to YAML config file [required]", metavar="FILE"),
  make_option(c("-r", "--redo"), action="store_true", default=FALSE,
              help="Re-run all analyses and override previous results")
)

# parse the arguments
opt <- parse_args(OptionParser(option_list=option_list))

# stop the code if config file is invalid
if (is.null(opt$config)) {
  stop("-c/--config is required.\n
        Usage: Rscript run_pipeline.R --config config.yaml")
}

if (!file.exists(opt$config)) {
  stop(paste("Config file not found:", opt$config))
}

# --- Load and validate config --------------------------------
cfg <- yaml::read_yaml(opt$config)

# set required parameters
required_fields <- c("codedir", "outdir", "fn_captus_sample_metadata", "fn_captus_target_metadata", "dir_fastq", "fn_eucs_metadata")
missing <- setdiff(required_fields, names(cfg))
if (length(missing) > 0) {
  stop(paste("Missing required config fields:", paste(missing, collapse=", ")))
}

# check if metadata files are invalid
if (!file.exists(path.expand(cfg$fn_captus_sample_metadata))) {
  stop(paste("fn_captus_sample_metadata file not found:", cfg$fn_captus_sample_metadata))
}

if (!file.exists(path.expand(cfg$fn_captus_target_metadata))) {
  stop(paste("fn_captus_target_metadata file not found:", cfg$fn_captus_target_metadata))
}

if (!file.exists(path.expand(cfg$fn_eucs_metadata))) {
  stop(paste("fn_eucs_metadata file not found:", cfg$fn_eucs_metadata))
}

# check if FASTQ directory is invalid
if (!dir.exists(path.expand(cfg$dir_fastq))) {
  stop(paste("dir_fastq directory not found:", cfg$dir_fastq))
}

# check if prefix is set, otherwise use the captus sample metadata filename
if (is.null(cfg$prefix) || cfg$prefix == "") {
  cfg$prefix <- "EucsPhylogenomics_output"
}

# --- Apply CLI overrides -------------------------------------
if (opt$redo) {
  message("Note: --redo flag set via CLI, overriding config.")
  cfg$redo <- TRUE
}

# --- Map config to rmarkdown params --------------------------
render_params <- list(
  codedir               = cfg$codedir,
  prefix                = cfg$prefix,
  outdir                = cfg$outdir,
  thread                = as.integer(f_get_param(cfg$thread, 1)),
  redo                  = as.logical(f_get_param(cfg$redo, FALSE)),

  exe_amas              = f_get_param(cfg$exe_amas, "AMAS.py"),
  exe_captus            = f_get_param(cfg$exe_captus, "captus"),
  exe_mafft             = f_get_param(cfg$exe_mafft, "mafft"),
  exe_fasttree          = f_get_param(cfg$exe_fasttree, "fasttree"),
  exe_astral4           = f_get_param(cfg$exe_astral4, "astral4"),
  exe_iqtree            = f_get_param(cfg$exe_iqtree, "iqtree3"),

  fn_captus_sample_metadata = cfg$fn_captus_sample_metadata,
  fn_captus_target_metadata = cfg$fn_captus_target_metadata,
  dir_fastq                 = cfg$dir_fastq,

  fn_species_tree     = f_get_param(cfg$fn_species_tree, ""),
  dir_locus_alignment = f_get_param(cfg$dir_locus_alignment, ""),

  run_phasing   = as.logical(f_get_param(cfg$run_phasing, FALSE)),
  exe_bwa       = f_get_param(cfg$exe_bwa, "bwa"),
  exe_samtools  = f_get_param(cfg$exe_samtools, "samtools"),
  exe_bcftools  = f_get_param(cfg$exe_bcftools, "bcftools"),
  exe_whatshap  = f_get_param(cfg$exe_whatshap, "whatshap"),

  fn_eucs_metadata = cfg$fn_eucs_metadata
)

# --- Run EucsPhylogenomics ------------------------------------
rmd_path <- file.path(path.expand(render_params$codedir), "codes", "1_main.Rmd")
if (!file.exists(rmd_path)) {
  stop(paste("1_main.Rmd not found:", rmd_path))
}

message("Starting EucsPhylogenomics pipeline...")
message("  Config:   ", opt$config)
message("  Prefix:   ", render_params$prefix)
message("  Output:   ", render_params$outdir)
message("  Samples:  ", render_params$fn_captus_sample_metadata)
message("  Threads:  ", render_params$thread)

# render the Rmarkdown file
rmarkdown::render(
  input       = rmd_path,
  params      = render_params,
  output_file = paste0(render_params$prefix, "_report.html"),
  output_dir  = file.path(path.expand(render_params$outdir), render_params$prefix),
  quiet       = FALSE
)

message("Done. Report: ",
        file.path(path.expand(render_params$outdir), render_params$prefix,
                  paste0(render_params$prefix, "_report.html")))
