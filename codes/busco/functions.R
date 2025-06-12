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
            
            if (!file.exists(fn_output_fasta) || is_redo) {
                # combine BUSCO into one file
                unlink(fn_output_fasta)
                f_fasta2msa(fn_input_fasta, sp, fn_output_fasta)
            }
        }
    }
}

# function: check BUSCO sequences
f_check_busco <- function(eucs_min_sp, non_eucs_min_sp, std_error, dir_output, thread, is_redo) {
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
        fn_output <- paste0(dir_output_filtered, busco, ".fna")
        if (file.exists(fn_output) && !is_redo) {
            return(NULL)
        }

        # delete the output file
        unlink(fn_output)

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
            seqinr::write.fasta(seq, names=names(seq), file.out=fn_output)
        }
    }

    stopCluster(nwcl)
}