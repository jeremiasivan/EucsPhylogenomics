# function: add FASTA sequence into one file
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

# function: run Easy353 to build database
f_easy353_build_db <- function(reftaxonomy, thread, dir_output, exe_build_db) {
    cmd_build_db <- paste(exe_build_db,
                          "-o", dir_output,
                          "-c", reftaxonomy,
                          "-t", thread,
                          "-generate")
    system(cmd_build_db)
}

# function: run Easy353
f_easy353 <- function(ls_fastq, dir_refdb, dir_output, exe_easy353) {
    cmd_easy353 <- paste(exe_easy353, "-1", ls_fastq[1])

    # add the second FASTQ file (if any)
    if (length(ls_fastq)>1) {
        cmd_easy353 <- paste(cmd_easy353, "-2", ls_fastq[2])
    }

    cmd_easy353 <- paste(cmd_easy353,
                         "-r", dir_refdb,
                         "-o", dir_output)
    system(cmd_easy353)
}

# function: run Hybpiper
f_hybpiper <- function(fn_target_gene, fn_fastq, prefix, dir_output, thread, exe_hybpiper) {
    hybpiper_cmd <- paste(exe_hybpiper, "assemble",
                          "-t_dna", fn_target_gene, 
                          "-r", fn_fastq,
                          "--prefix", prefix,
                          "-o", dir_output,
                          "--bwa",
                          "--cpu", thread)
    system(hybpiper_cmd)
}

# function: run Orthofinder
f_orthofinder <- function(dir_fasta, dir_output, thread, exe_orthofinder) {
    cmd_orthofinder <- paste(exe_orthofinder,
                         "-f", dir_fasta,
                         "-o", dir_output,
                         "-t", thread)
    system(cmd_orthofinder)
}

# function: run ROADIES
f_roadies <- function(fn_config, thread, exe_roadies) {
    cmd_roadies <- paste("python", exe_roadies,
                         "--config", fn_config,
                         "--cores", thread, "--converge")
    system(cmd_roadies)
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