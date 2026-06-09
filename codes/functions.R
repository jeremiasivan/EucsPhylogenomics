# functions for running EucsPhylogenomics

# function: run CAPTUS to clean FASTQ files
f_run_captus_cleaning <- function(exe_captus, input_dir, output_dir, thread) {
  captus_cmd <- paste(exe_captus, "clean -r", input_dir, "-o", output_dir, "--threads", thread, "--overwrite")
  system(captus_cmd)
}

# function: run CAPTUS to assemble short reads
f_run_captus_assembly <- function(exe_captus, input_dir, output_dir, thread) {
  captus_cmd <- paste(exe_captus, "assemble -r", input_dir, "-o", output_dir, "--threads", thread, "--overwrite")
  system(captus_cmd)
}

# function: run CAPTUS to extract target loci
f_run_captus_extraction <- function(exe_captus, fn_target_loci, input_dir, output_dir, thread) {
  captus_cmd <- paste(exe_captus, "extract -f", fn_target_loci, "-a", input_dir, "-o", output_dir, "--threads", thread, "--overwrite")
  system(captus_cmd)
}

# function: run CAPTUS to align extracted loci
f_run_captus_alignment <- function(exe_captus, input_dir, output_dir, thread) {
  captus_cmd <- paste(exe_captus, "align -e", input_dir, "-o", output_dir, "--threads", thread, "--overwrite")
  system(captus_cmd)
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

# function: run AMAS
f_amas_dna <- function(input_regex, fn_output_aln, fn_output_part, exe_amas) {
    cmd_amas <- paste("python3", exe_amas, "concat -f fasta -d dna",
                      "-i", input_regex,
                      "--part-format nexus",
                      "-t", fn_output_aln,
                      "-p", fn_output_part)
    system(cmd_amas)
}

# function: run FastTree
f_fasttree <- function(fn_fasta, fn_output, exe_fasttree) {
    cmd_fasttree <- paste(exe_fasttree, "-gtr", fn_fasta, ">", fn_output)
    system(cmd_fasttree)
}

# function: run MAFFT --keeplength
f_mafft_keeplen <- function(fn_input, fn_alignment, fn_output, exe_mafft) {
    cmd_mafft <- paste(exe_mafft, "--add",
                       fn_input, "--keeplength", fn_alignment,
                       ">", fn_output)
    system(cmd_mafft)
}

# function: run IQ-Tree 2
f_iqtree2 <- function(fn_input, thread, exe_iqtree2) {
    cmd_iqtree2 <- paste(exe_iqtree2,
                         "-s", fn_input,
                         "-B 1000",
                         "-T", thread, "--quiet -redo")
    system(cmd_iqtree2)
}

# function: run ASTRAL-IV
f_astral4 <- function(fn_input, fn_output, fn_log, exe_astral) {
    cmd_astral <- paste(exe_astral,
                        "-u 2",
                        "-i", fn_input,
                        "-o", fn_output,
                        "2>", fn_log)
    system(cmd_astral)
}