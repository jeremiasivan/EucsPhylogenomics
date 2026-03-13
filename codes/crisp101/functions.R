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

# run MAFFT
f_mafft <- function(fn_input, fn_out, exe_mafft) {
  cmd_mafft <- paste(exe_mafft, "--retree 2", fn_input, ">", fn_out)
  system(cmd_mafft)
}

# run MAFFT --add
f_mafft_add <- function(fn_ref, fn_sample, fn_out, exe_mafft) {
  cmd_mafft <- paste(exe_mafft, "--addfull", fn_sample, "--keeplength", fn_ref, ">", fn_out)
  system(cmd_mafft)
}

# function: delete sequence with >=threshold gaps
f_remove_seq <- function(fn_fasta, fn_output, threshold) {
    # read the DNA alignment
    seq <- Biostrings::readBStringSet(fn_fasta, format="fasta")

    # iterate over sequences
    pl <- c()
    for (i in 1:length(seq)) {
        seq_chr <- as.character(seq[i])
        
        # count the proportion of gaps
        n_gaps <- stringr::str_count(seq_chr, "-")
        n_total <- stringr::str_count(seq_chr)
        prop_gaps <- n_gaps / n_total

        if (prop_gaps >= threshold) {
            pl <- c(pl, i)
        }
    }
   
    # remove sequences with >=50% gaps
    if (length(pl)>0) {
      seq <- seq[-pl]
    }

    # save the new DNA alignment
    Biostrings::writeXStringSet(seq, filepath=fn_output)
}

# function: run IQ-TREE 2
f_iqtree2 <- function(fn_input, exe_iqtree2) {
    cmd_iqtree2 <- paste(exe_iqtree2,
                         "-s", fn_input,
                         "-bb 1000",
                         "-T 1 --quiet -redo")
    system(cmd_iqtree2)
}

# function: run ASTRAL-III 
f_astral <- function(fn_input, fn_output, fn_log, exe_astral) {
    cmd_astral <- paste("java -jar", exe_astral,
                    "-i", fn_input,
                    "-o", fn_output,
                    "-t 2 2>", fn_log)
    system(cmd_astral)
}