# function: add entry to logfile
f_write_log <- function(fn_log, msg) {
    write.table(msg, file=fn_log, quote=F, row.names=F, col.names=F, append=T)
}

# function: create Mash sketch for all reference alignments
f_mash_sketch <- function(exe_mash, str_refseq, thread, prefix) {
    cmd <- paste(exe_mash, "sketch",
                "-o", prefix,
                "-p", thread,
                str_refseq)
    system(cmd)
}

# function: run Mash screen on short reads
f_mash_screen <- function(exe_mash, file_msh, file_input, file_output) {
    cmd_mash <- paste(exe_mash, "screen -w", file_msh, file_input, ">", file_output)
    system(cmd_mash)
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
f_iqtree2 <- function(fn_input, prefix, exe_iqtree2) {
    cmd_iqtree2 <- paste(exe_iqtree2,
                         "-s", fn_input,
                         "--prefix", prefix,
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

# function: quality control using AdapterRemoval
f_qc_adapterremoval <- function(fastq_one, fastq_two, fn_adapters, prefix, min_quality, thread, exe_adapterremoval) {
    cmd_qc <- paste(exe_adapterremoval, "--file1", fastq_one)
    
    # check the reverse fastq file
    if (!is.null(fastq_two)) {
        cmd_qc <- paste(cmd_qc, "--file2", fastq_two)
    }

    # update the adapters
    if (fn_adapters != "" && file.exists(fn_adapters)) {
        cmd_qc <- paste(cmd_qc, "--adapter-list", fn_adapters)
    }

    # set the minimum quality score and merge overlapping reads
    cmd_qc <- paste(cmd_qc, "--basename", prefix,
                    "--trimqualities --trimns --minquality", min_quality)
    
    # run AdapterRemoval to get prefix.collapsed.truncated
    system(cmd_qc)
}

# function: quality control using BBTools
f_qc_bbtools <- function(fastq, dir_output, dir_rqcfilterdata, exe_rqcfilter2) {
    cmd_qc <- paste0(exe_rqcfilter2,
                     " -Xmx=101077m jni=t",
                     " in=",fastq,
                     " path=",dir_output,
                     " barcodefilter=f clumpify=t dedupe=t kapa=f khist=t maxns=1 minlen=49 mlf=0.33 phix=f pigz=t pjet=f qtrim=r removecat=f removedog=f removehuman=f removemicrobes=f removemouse=f rna=f",
                     " rqcfilterdata=",dir_rqcfilterdata, 
                     " sketch skipfilter=t trimfragadapter=t trimpolyg=5 trimq=6 unpigz=t usejni=f")
    system(cmd_qc)
}

# function: convert the partitions from number to species
f_part_num2chr <- function(input) {
    # extract labels
    ls_label <- attr(input, "labels")

    # iterate over the partitions
    ls_partition <- lapply(input, function(x) {
        # convert the numbers to characters
        ls_species <- sapply(x, function(y) {
            ls_label[y]
        })
    })

    # add frequency
    attr(ls_partition, "number") <- attr(input, "number")

    return(ls_partition)
}

# function: compare two partitions
f_compare_parts <- function(input, ref_input) {
    # extract partitions
    input_names <- unlist(lapply(input, function(x) { paste0(sort(x), collapse="-") }))
    ref_names <- unlist(lapply(ref_input, function(x) { paste0(sort(x), collapse="-") }))

    return(ref_names[!ref_names%in%input_names])
}

# function: run IQ-TREE2 with contrained topology
f_iqtree2_constrained <- function(fn_input, fn_topology, prefix, exe_iqtree2) {
    cmd_iqtree2 <- paste(exe_iqtree2,
                         "-s", fn_input,
                         "-g", fn_topology,
                         "--prefix", prefix,
                         "-T 1 --quiet -redo")
    system(cmd_iqtree2)
}

# function: run AU test
f_iqtree2_au <- function(fn_input, fn_topology, prefix, exe_iqtree2) {
    cmd_iqtree2 <- paste(exe_iqtree2,
                         "-s", fn_input,
                         "-z", fn_topology,
                         "--prefix", prefix,
                         "-n 0 -zb 10000 -au",
                         "-T 1 --quiet -redo")
    system(cmd_iqtree2)
}