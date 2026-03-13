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
