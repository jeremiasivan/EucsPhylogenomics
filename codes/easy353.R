# Run Easy353
dir_output <- "/data/jeremias/eucs/easy353/"
thread <- 10

exe_build_db <- "build_database.py"
exe_easy353 <- "easy353.py"
exe_mafft <- "mafft"
exe_iqtree2 <- "iqtree2"

dir_shortreads <- "/data/jeremias/eucs/shortreads/"

reftaxonomy <- "Eucalyptus"

################################################

library(doSNOW)

################################################

# function: run Easy353 to build database
f_build_db <- function(reftaxonomy, thread, dir_output, exe_build_db) {
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

################################################

# create output directories
dir_output_db <- paste0(dir_output, "/", tolower(reftaxonomy), "_db/")
if (!dir.exists(dir_output_db)) {
    dir.create(dir_output_db, recursive=T)
}

# build database
f_build_db(reftaxonomy, thread, dir_output_db, exe_build_db)

# extract reference sequences
ls_refseq <- list.files(dir_output_db, pattern="*.fasta$")
ls_refseq <- sapply(ls_refseq, function(x) {
    ls_element <- strsplit(x, split="\\.")
    ls_element[[1]][length(ls_element[[1]])-2]
})

# iterate over short-reads
dir_output_easy353 <- paste0(dir_output, "/easy353/")

ls_shortreads <- list.dirs(dir_shortreads, recursive=F, full.names=F)
for (shortread in ls_shortreads) {
    # extract FASTQ files
    ls_fastq <- c(paste0(dir_shortreads, "/", shortread, "_1.fastq"),
                  paste0(dir_shortreads, "/", shortread, "_2.fastq"))

    if (!file.exists(ls_fastq[1])) {
        next
    }

    # remove reverse FASTQ file if unavailable
    if (!file.exists(ls_fastq[2])) {
        ls_fastq <- ls_fastq[1]
    }

    # create output directory
    dir_output_easy353_sp <- paste0(dir_output_easy353, shortread)
    if (!dir.exists(dir_output_easy353_sp)) {
        dir.create(dir_output_easy353_sp, recursive=T)
    }

    # run Easy353
    f_easy353(ls_fastq, dir_output_db, dir_output_easy353_sp, exe_easy353)
}

# extract concatenated sequence
dir_output_tree <- paste0(dir_output, "/tree/")
if (!dir.exists(dir_output_tree)) {
    dir.create(dir_output_tree, recursive=T)
}

# create doSNOW cluster
nthread <- thread
nthread <- ifelse(thread > nthread, nthread, thread)

nwcl <- makeCluster(nthread)
doSNOW::registerDoSNOW(nwcl)

# iterate over short-reads
foreach (ref = ls_refseq) %dopar% {
    # create output directory
    dir_output_tree_sp <- paste0(dir_output_tree, ref, "/")
    if (!dir.exists(dir_output_tree_sp)) {
        dir.create(dir_output_tree_sp, recursive=T)
    }

    # output files
    fn_fasta_concat <- paste0(dir_output_tree_sp, "concat.fa")
    fn_fasta_concat_aligned <- paste0(dir_output_tree_sp, "concat_aligned.fa")

    for (shortread in ls_shortreads) {
        # input file
        dir_assembly <- paste0(dir_output_easy353, "/", shortread, "/assemble_out/")
        fn_fasta <- list.files(dir_assembly, pattern=paste0("*",ref,".a353.fasta"))

        # add sequence into one file
        f_fasta2msa(fn_fasta, shortread, fn_fasta_concat)
    }

    # run MAFFT using FFT-NS-2
    f_mafft(fn_fasta_concat, fn_fasta_concat_aligned, "--retree 2", exe_mafft)

    # run IQ-Tree 2
    f_iqtree2(fn_fasta_concat_aligned, exe_iqtree2)
}

stopCluster(nwcl)

################################################