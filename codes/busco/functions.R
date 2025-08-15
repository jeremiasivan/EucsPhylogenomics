# functions for codes/busco596

# function: combine individual FASTA as MSA
f_fasta2msa <- function(fn_input, header, fn_out) {
    # initiate variable
    first_sequence <- TRUE

    # open the FASTA file
    con <- file(fn_input, "r")

    # iterate over lines
    while (length(line <- readLines(con, n = 1)) > 0) {
        if (grepl("^>+", line)) {
            if (first_sequence) {
                write.table(paste0(">", header), file=fn_out, quote=F, row.names=F, col.names=F, append=T)
                first_sequence <- FALSE
            }
        } else {
            write.table(line, file=fn_out, quote=F, row.names=F, col.names=F, append=T)
        }
    }
    
    # close the file connection
    close(con)
}

# function: extract BUSCO
f_extract_busco <- function(ls_species, dir_busco, lineage, dir_output, is_redo) {
    # create output directory
    dir_output_all <- paste0(dir_output, "/all/")
    if (!dir.exists(dir_output_all)) {
        dir.create(dir_output_all, recursive=T)
    }

    # iterate over reference for Angophora and Corymbia
    for (sp in ls_species) {
        # BUSCO output directory
        dir_busco_fna <- paste0(dir_busco, "/", sp, "/run_", lineage, "/busco_sequences/single_copy_busco_sequences/")

        # list BUSCO
        ls_busco <- list.files(dir_busco_fna, pattern="*.fna$", recursive=F, full.names=F)
        ls_busco <- gsub(".fna", "", ls_busco)

        # iterate over BUSCO
        for (busco in ls_busco) {
            fn_input_fasta <- paste0(dir_busco_fna, busco, ".fna")
            fn_output_fasta <- paste0(dir_output_all, busco, ".fna")
        
            # combine BUSCO into one file
            f_fasta2msa(fn_input_fasta, sp, fn_output_fasta)
        }
    }
}

# function: check BUSCO sequences
f_check_busco <- function(eucs_min_sp, non_eucs_min_sp, std_error, dir_output, thread, is_redo, exe_mafft) {
    # create output directory
    dir_output_filtered <- paste0(dir_output, "/filtered/")
    if (!dir.exists(dir_output_filtered)) {
        dir.create(dir_output_filtered, recursive=T)
    }

    # list BUSCO
    dir_output_all <- paste0(dir_output, "/all/")
    ls_busco <- list.files(dir_output_all, pattern="*.fna$", recursive=F, full.names=F)
    ls_busco <- gsub(".fna", "", ls_busco)

    # create doSNOW cluster
    nwcl <- makeCluster(thread)
    doSNOW::registerDoSNOW(nwcl)

    # iterate over BUSCO
    foreach (busco = ls_busco) %dopar% {
        fn_input <- paste0(dir_output_all, busco, ".fna")
        fn_output <- paste0(dir_output_filtered, "step1_", busco, ".fna")
        fn_output_aligned <- paste0(dir_output_filtered, "step1_", busco, "_aligned.fna")
        if (all(file.exists(fn_output, fn_output_aligned)) && !is_redo) {
            return(NULL)
        }

        # delete the output file
        unlink(fn_output)
        unlink(fn_output_aligned)

        # open the sequence
        seq <- seqinr::read.fasta(fn_input)
        ls_len <- sapply(seq, function(x) { length(x) })

        median_len <- median(ls_len)
        max_len <- median_len + std_error*median_len
        min_len <- median_len - std_error*median_len
        
        # filter out sequences that are too short or long
        ls_idx_filtered <- c()
        for (i in 1:length(seq)) {
            if (length(seq[[i]]) >= min_len && length(seq[[i]]) <= max_len) {
                ls_idx_filtered <- c(ls_idx_filtered, i)
            }
        }

        # extract sequences
        filtered_seq <- seq[ls_idx_filtered]

        ls_filtered_sp <- names(filtered_seq)
        ls_filtered_eucs <- ls_filtered_sp[grepl("^E_",ls_filtered_sp)]
        ls_filtered_non_eucs <- ls_filtered_sp[!ls_filtered_sp%in%ls_filtered_eucs]

        # check if the number of species suffices
        if (length(ls_filtered_eucs) >= eucs_min_sp && length(ls_filtered_non_eucs) >= non_eucs_min_sp) {
            seqinr::write.fasta(filtered_seq, names=names(filtered_seq), file.out=fn_output)

            # run MAFFT
            system(paste(exe_mafft, "--retree 2", fn_output, ">", fn_output_aligned))
        }
    }

    stopCluster(nwcl)
}

# function: run TAPER
f_taper_dna <- function(fn_fasta, fn_output, exe_taper) {
    cmd_taper <- paste("julia", exe_taper, "-m -", fn_fasta, ">", fn_output)
    system(cmd_taper)
}

# function: convert "N" to gaps
f_n2gap <- function(fn_fasta, fn_output) {
    # read the DNA alignment
    seq <- Biostrings::readDNAStringSet(fn_fasta, format="fasta")

    # convert "N" to gaps
    for (i in 1:length(seq)) {
        seq[i] <- gsub("N", "-", seq[i])
    }

    # save the new DNA alignment
    Biostrings::writeXStringSet(seq, filepath=fn_output)
}

# function: delete sequence with >=threshold gaps
f_remove_seq <- function(fn_fasta, fn_output, threshold) {
    # read the DNA alignment
    seq <- Biostrings::readDNAStringSet(fn_fasta, format="fasta")

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
    seq <- Biostrings::readDNAStringSet(fn_fasta, format="fasta")
    seq <- as.character(seq)

    # convert sequence into matrix
    seq_matrix <- do.call(rbind, strsplit(seq, ""))

    # count the number of sequences
    n_seqs <- nrow(seq_matrix)

    # function to check gap proportion per column
    is_valid_site <- function(column) {
        gap_count <- sum(column == "-")
        return((gap_count / n_seqs) <= threshold)
    }

    # extract columns with <50% gaps
    valid_sites <- apply(seq_matrix, 2, is_valid_site)
    seq_matrix_filtered <- seq_matrix[, valid_sites]
    seq_matrix_filtered <- Biostrings::DNAStringSet(apply(seq_matrix_filtered, 1, paste0, collapse=""))

    # save the new DNA alignment
    Biostrings::writeXStringSet(seq_matrix_filtered, filepath=fn_output)
}

# function: run FastTree
f_fasttree <- function(fn_fasta, fn_output, exe_fasttree) {
    cmd_fasttree <- paste(exe_fasttree, fn_fasta, ">", fn_output)
    system(cmd_fasttree)
}

# function: run TreeShrink
f_treeshrink <- function(fn_input, prefix, dir_output, exe_treeshrink) {
    cmd_treeshrink <- paste("python", exe_treeshrink,
                            "-t", fn_input,
                            "-O", prefix,
                            "-o", dir_output)
    system(cmd_treeshrink)
}

# function: run MAFFT
f_mafft <- function(fn_input, fn_output, params_mafft, exe_mafft) {
    cmd_mafft <- paste(exe_mafft, params_mafft,
                       fn_input, ">", fn_output)
    system(cmd_mafft)
}

# function: run IQ-Tree 2
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