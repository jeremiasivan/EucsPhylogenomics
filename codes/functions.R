# functions for running EucsPhylogenomics

# function: run CAPTUS to clean FASTQ files
f_run_captus_cleaning <- function(exe_captus, input_dir, output_dir, thread, is_redo) {
  captus_cmd <- paste(exe_captus, "clean -r", input_dir, "-o", output_dir, "--threads", thread)

  if (is_redo) {
    captus_cmd <- paste(captus_cmd, "--overwrite")
  }

  system(captus_cmd)
}

# function: run CAPTUS to assemble short reads
f_run_captus_assembly <- function(exe_captus, input_dir, output_dir, thread, is_redo) {
  captus_cmd <- paste(exe_captus, "assemble -r", input_dir, "-o", output_dir, "--threads", thread)

  if (is_redo) {
    captus_cmd <- paste(captus_cmd, "--overwrite")
  }
  
  system(captus_cmd)
}

# function: run CAPTUS to extract target loci
f_run_captus_extraction <- function(exe_captus, fn_target_loci, input_dir, output_dir, thread, is_redo) {
  captus_cmd <- paste(exe_captus, "extract -f", fn_target_loci, "-a", input_dir, "-o", output_dir, "--threads", thread)

  if (is_redo) {
    captus_cmd <- paste(captus_cmd, "--overwrite")
  }

  system(captus_cmd)
}

# function: run CAPTUS to align extracted loci
f_run_captus_alignment <- function(exe_captus, input_dir, output_dir, thread, is_redo) {
  captus_cmd <- paste(exe_captus, "align -e", input_dir, "-o", output_dir, "--threads", thread)

  if (is_redo) {
    captus_cmd <- paste(captus_cmd, "--overwrite")
  }

  system(captus_cmd)
}

# function: add new sequences from CAPTUS to existing alignment
f_add_captus_seq <- function(fn_fasta, fn_captus, fn_output) {
    # read the existing DNA alignment
    seq <- Biostrings::readBStringSet(fn_fasta, format="fasta")

    # read the new sequences from CAPTUS
    seq_captus <- Biostrings::readBStringSet(fn_captus, format="fasta")

    # combine the two sets of sequences
    seq_combined <- c(seq, seq_captus)

    # save the new DNA alignment
    Biostrings::writeXStringSet(seq_combined, filepath=fn_output)
}

# run MAFFT --add
f_mafft_add <- function(fn_ref, fn_sample, fn_out, exe_mafft) {
  cmd_mafft <- paste(exe_mafft, "--auto --add", fn_sample, "--keeplength", fn_ref, ">", fn_out)
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
        if (stringr::str_count(seq_chr,"-") >= threshold*stringr::str_count(seq_chr)) {
            pl <- c(pl, i)
        }
    };
   
    # remove sequences with >=50% gaps
    if (length(pl)>0) {
      seq <- seq[-pl]
    }

    # save the new DNA alignment
    Biostrings::writeXStringSet(seq, filepath=fn_output)
}

# function: delete columns with >=threshold gaps
f_remove_col <- function(fn_fasta, fn_output, threshold) {
    # read the multiple sequence alignment
    seq <- Biostrings::readBStringSet(fn_fasta, format="fasta")
    seq <- as.character(seq)

    # convert sequence into matrix
    seq_matrix <- do.call(rbind, strsplit(seq, ""))

    # count the number of sequences
    n_seqs <- nrow(seq_matrix)

    # function to check gap proportion per column
    is_valid_site <- function(column) {
        gap_count <- sum(column == "-")
        return((gap_count / n_seqs) < threshold)
    }

    # extract columns with <50% gaps
    valid_sites <- apply(seq_matrix, 2, is_valid_site)
    seq_matrix_filtered <- seq_matrix[, valid_sites]
    seq_matrix_filtered <- Biostrings::DNAStringSet(apply(seq_matrix_filtered, 1, paste0, collapse=""))

    # save the new DNA alignment
    Biostrings::writeXStringSet(seq_matrix_filtered, filepath=fn_output)
}

# function: trim CAPTUS labels
f_trim_captus_label <- function(fn_input, locus, fn_output) {
    # read the DNA alignment
    seq <- Biostrings::readBStringSet(fn_input)

    # extract CAPTUS headers (source: Claude)
    df_output <- rbindlist(lapply(names(seq), function(header) {
      # extract sequence ID
      seq_id <- trimws(str_extract(header, "^\\S+"))
      
      # extract all [key=value] pairs
      matches <- str_match_all(header, "\\[([^=]+)=([^]]+)\\]")[[1]]
      if (nrow(matches) == 0) {
        return(NULL)
      }
      
      # matches is a matrix: col1=full match, col2=key, col3=value
      pairs <- setNames(as.list(trimws(matches[, 3])), trimws(matches[, 2]))
      
      c(list(seq_id = seq_id, locus = locus), pairs)
    }))

    # update sequence headers
    names(seq) <- sapply(names(seq), function(x) { unlist(strsplit(x, split=" "))[1] })

    # save the new DNA alignment
    Biostrings::writeXStringSet(seq, filepath=fn_output)

    return(df_output)
}

# function: run AMAS
f_amas_dna <- function(input_regex, fn_output_aln, fn_output_part, exe_amas) {
    cmd_amas <- paste("python3", exe_amas, "concat -f fasta -d dna",
                      "-i", input_regex,
                      "--part-format nexus",
                      "-t", fn_output_aln,
                      "-p", fn_output_part)
    system(cmd_amas)
}

# function: convert "N" and "?" to gaps
f_unknown2gap <- function(fn_fasta, fn_output) {
    # read the DNA alignment
    seq <- Biostrings::readBStringSet(fn_fasta, format="fasta")

    # convert "N" to gaps
    for (i in 1:length(seq)) {
        seq[i] <- gsub("N", "-", seq[i])
        seq[i] <- gsub("\\?", "-", seq[i])
    }

    # save the new DNA alignment
    Biostrings::writeXStringSet(seq, filepath=fn_output)
}

# function: run FastTree
f_fasttree <- function(fn_fasta, fn_output, exe_fasttree) {
    cmd_fasttree <- paste(exe_fasttree, "-gtr", fn_fasta, ">", fn_output)
    system(cmd_fasttree)
}

# function: run IQ-Tree 2
f_iqtree2 <- function(fn_input, fn_tree, fn_partition, prefix, thread, is_redo, exe_iqtree2) {
    cmd_iqtree2 <- paste(exe_iqtree2,
                         "-s", fn_input,
                         "-t", fn_tree,
                         "-p", fn_partition,
                         "--prefix", prefix,
                         "-T", thread, "--quiet")
                         
    if (is_redo) {
        cmd_iqtree2 <- paste(cmd_iqtree2, "-redo")
    }

    system(cmd_iqtree2)
}

# function: run IQ-Tree 2 (-T 1)
f_iqtree2_singlethread <- function(fn_input, prefix, is_redo, exe_iqtree2) {
    cmd_iqtree2 <- paste(exe_iqtree2,
                         "-s", fn_input,
                         "--prefix", prefix,
                         "-T 1 --quiet")
    
    if (is_redo) {
        cmd_iqtree2 <- paste(cmd_iqtree2, "-redo")
    }

    system(cmd_iqtree2)
}

# function: run ASTRAL-IV
f_astral4 <- function(fn_input, fn_tree, fn_output, fn_log, thread, exe_astral) {
    cmd_astral <- paste(exe_astral,
                        "-u 2",
                        "-i", fn_input,
                        "-c", fn_tree,
                        "-o", fn_output,
                        "-t", thread,
                        "2>", fn_log)
    system(cmd_astral)
}

# function: calculate sCF and gCF
f_calculate_cf <- function(fn_all_trees, fn_sp_tree, dir_fasta, prefix, thread, is_redo, exe_iqtree2) {
    cmd_cf <- paste(exe_iqtree2,
                    "-t", fn_sp_tree,
                    "--gcf", fn_all_trees,
                    "-p", dir_fasta,
                    "--scf 100",
                    "-T", thread,
                    "--prefix", prefix)

    if (is_redo) {
        cmd_cf <- paste(cmd_cf, "-redo")
    }
    
    system(cmd_cf)
}
