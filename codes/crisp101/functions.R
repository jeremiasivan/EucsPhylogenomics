# functions for codes/crisp101

# run HybPiper
f_hybpiper <- function(sample, dir_sample, target_file, data_type, thread, outdir, exe_hybpiper) {
  cmd_hybpiper <- paste(exe_hybpiper, "assemble")
  
  # check data type
  if (data_type == "dna") {
    cmd_hybpiper <- paste(cmd_hybpiper, "-t_dna", target_file)
  } else if (data_type == "aa") {
    cmd_hybpiper <- paste(cmd_hybpiper, "-t_aa", target_file)
  } else {
    return(paste0("Invalid data type for ", sample, ". Skipped."))
  }
  
  # check if FASTQ files exist
  ls_fastq <- list.files(dir_sample, pattern=paste0(sample, "_R[12].fastq.gz"))
  if (length(ls_fastq) == 0) {
    return(paste0("No FASTQ files found for ", sample, ". Skipped."))
  }

  # run HybPiper
  cmd_hybpiper <- paste(cmd_hybpiper,
                        "-r", paste0(dir_sample, "/", sample, "_R*.fastq.gz"),
                        "--prefix", sample,
                        "--hybpiper_output", outdir,
                        "--bwa --cpu", thread)
  system(cmd_hybpiper)

  return(paste0("HybPiper run finished for ", sample, "."))
}

# check HybPiper paralogs
f_check_hybpiper_paralogs <- function(sample, outdir) {
  # initiate variable
  ls_paralogs_final <- c()

  # read files containing paralogs
  fn_paralogs <- paste0(outdir, "/", sample, "_genes_with_long_paralog_warnings.txt")
  if (file.exists(fn_paralogs)) {
      ls_paralogs <- readLines(fn_paralogs)
      ls_paralogs <- sapply(ls_paralogs, function(x) { unlist(strsplit(x, "-"))[2] })
      ls_paralogs_final <- c(ls_paralogs_final, ls_paralogs)
  }

  fn_paralogs_depth <- paste0(outdir, "/", sample, "_genes_with_paralog_warnings_by_contig_depth.csv")
  if (file.exists(fn_paralogs_depth) && file.size(fn_paralogs_depth)>0) {
      ls_paralogs <- data.table::fread(fn_paralogs_depth)
      ls_paralogs <- ls_paralogs$V2[ls_paralogs$V4]
      ls_paralogs <- gsub("gene ", "", ls_paralogs)
      ls_paralogs_final <- c(ls_paralogs_final, ls_paralogs)
  }

  return(unique(ls_paralogs_final))
}