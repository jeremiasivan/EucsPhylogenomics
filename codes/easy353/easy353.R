################################################
#############        EASY353       #############
################################################

dir_codes <- "/home/jeremias/EucsPhylogenomics/codes/"
dir_output <- "/data/jeremias/eucs/easy353/"
thread <- 10

# run Easy353
exe_build_db <- "build_database.py"
exe_easy353 <- "easy353.py"
exe_mafft <- "mafft"
exe_iqtree2 <- "iqtree2"

reftaxonomy <- "Eucalyptus"
dir_shortreads <- "/data/jeremias/eucs/shortreads/"

################################################

library(doSNOW)

source(paste0(dir_codes, "/functions.R"))

################################################

# create output directories
dir_output_db <- paste0(dir_output, "/", tolower(reftaxonomy), "_db/")
if (!dir.exists(dir_output_db)) {
    dir.create(dir_output_db, recursive=T)
}

# build database
if (!dir.exists(paste0(dir_output_db, "353gene/"))) {
    f_easy353_build_db(reftaxonomy, thread, dir_output_db, exe_build_db)
}

# extract reference sequences
ls_refseq <- list.files(dir_output_db, pattern="*.fasta$")
ls_refseq <- sapply(ls_refseq, function(x) {
    ls_element <- strsplit(x, split="\\.")
    ls_element[[1]][length(ls_element[[1]])-2]
})

################################################

# iterate over short-reads
dir_output_easy353 <- paste0(dir_output, "/easy353/")

ls_shortreads <- list.dirs(dir_shortreads, recursive=F, full.names=F)
for (shortread in ls_shortreads) {
    # extract FASTQ files
    ls_fastq <- c(paste0(dir_shortreads, "/", shortread, "/", shortread, ".1.fastq"),
                  paste0(dir_shortreads, "/", shortread, "/", shortread, ".2.fastq"))

    if (!file.exists(ls_fastq[1])) {
        next
    }

    # remove reverse FASTQ file if unavailable
    if (!file.exists(ls_fastq[2])) {
        ls_fastq <- ls_fastq[1]
    }

    # create output directory
    dir_output_easy353_sp <- paste0(dir_output_easy353, shortread, "/")
    if (!dir.exists(dir_output_easy353_sp)) {
        dir.create(dir_output_easy353_sp, recursive=T)
    }

    # check if log file exists
    fn_log <- paste0(dir_output_easy353_sp, "easy353.log")
    if (file.exists(fn_log)) {
        next
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
nthread <- ifelse(thread > length(ls_refseq), length(ls_refseq), thread)

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
    fn_fasta_concat_aligned_treefile <- paste0(dir_output_tree_sp, "concat_aligned.fa.treefile")

    # check if treefile exists
    if (file.exists(fn_fasta_concat_aligned_treefile)) {
        return(NULL)
    }

    for (shortread in ls_shortreads) {
        # input file
        dir_assembly <- paste0(dir_output_easy353, shortread, "/assemble_out/")
        fn_fasta <- list.files(dir_assembly, pattern=paste0("*",ref,".a353.fasta"), recursive=F, full.names=F)

        if (!file.exists(fn_fasta)) {
            next
        }

        # add sequence into one file
        f_fasta2msa(paste0(dir_assembly, fn_fasta), shortread, fn_fasta_concat)
    }

    # run MAFFT using FFT-NS-2
    f_mafft(fn_fasta_concat, fn_fasta_concat_aligned, "--retree 2", exe_mafft)

    # run IQ-Tree 2
    f_iqtree2(fn_fasta_concat_aligned, exe_iqtree2)
}

stopCluster(nwcl)

################################################