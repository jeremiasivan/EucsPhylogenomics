# Run Hybpiper
dir_codes <- "/home/jeremias/EucsPhylogenomics/codes/"
dir_output <- "/data/jeremias/eucs/hybpiper/"
thread <- 50

exe_hybpiper <- "hybpiper"
exe_mafft <- "mafft"
exe_iqtree2 <- "iqtree2"
exe_astral <- "astral"

fn_target_gene <- "target_gene.fa"

dir_shortreads <- "/data/jeremias/eucs/shortreads/"

################################################

library(doSNOW)

source(paste0(dir_codes, "/codes/functions.R"))

################################################

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

################################################

# create output directory
dir_output_hybpiper <- paste0(dir_output, "/hybpiper/")
if (!dir.exists(dir_output_hybpiper)) {
    dir.create(dir_output_hybpiper, recursive=T)
}

# iterate over short-reads
ls_shortreads <- list.dirs(dir_shortreads, recursive=F, full.names=F)
for (shortread in ls_shortreads) {
    # output file
    fn_output <- paste0(dir_output_hybpiper, shortread, "/target_tallies.txt")
    if (file.exists(fn_output)) {
        next
    }

    # extract FASTQ files
    ls_fastq <- paste0(dir_shortreads, "/", shortread, "/", shortread, ".*.fastq")

    # run Easy353
    f_hybpiper(fn_target_gene, ls_fastq, shortread, dir_output_hybpiper, thread, exe_hybpiper)
}

# extract concatenated sequence
dir_output_tree <- paste0(dir_output, "/tree/")
if (!dir.exists(dir_output_tree)) {
    dir.create(dir_output_tree, recursive=T)
}

# create doSNOW cluster
nwcl <- makeCluster(thread)
doSNOW::registerDoSNOW(nwcl)

# iterate over short-reads
for (shortread in ls_shortreads) {
    # output directory
    dir_output_hybpiper_sp <- paste0(dir_output_hybpiper, "/", shortread, "/")

    # extract genes
    ls_gene_sp <- list.dirs(dir_output_hybpiper_sp, recursive=F, full.names=F)

    foreach (gene = ls_gene_sp) %dopar% {
        # output directory
        dir_output_gene <- paste0(dir_output_tree, gene, "/")
        if (!dir.exists(dir_output_gene)) {
            dir.create(dir_output_gene, recursive=T)
        }

        # input file
        fn_fasta <- paste0(dir_output_hybpiper_sp, gene, "/", shortread, "/sequences/FNA/", gene, ".FNA")

        # output file
        fn_fasta_concat <- paste0(dir_output_gene, "concat.fa")
        fn_fasta_concat_aligned <- paste0(dir_output_tree, gene, "/concat_aligned.fa")
        
        # check if MSA exists
        if (file.exists(fn_fasta_concat_aligned)) {
            return(NULL)
        }

        # add sequence into one file
        f_fasta2msa(fn_fasta, shortread, fn_fasta_concat)
    }
}

# extract all genes
ls_genes <- list.dirs(dir_output_tree, recursive=F, full.names=F)

# iterate over genes
ls_trees <- foreach (gene = ls_genes, .combine='c') %dopar% {
    # output files
    fn_fasta_concat <- paste0(dir_output_tree, gene, "/concat.fa")
    fn_fasta_concat_aligned <- paste0(dir_output_tree, gene, "/concat_aligned.fa")
    fn_fasta_concat_aligned_treefile <- paste0(dir_output_tree, gene, "/concat_aligned.fa.treefile")

    # check if treefile exists
    if (file.exists(fn_fasta_concat_aligned_treefile)) {
        return(fn_tree)
    }

    # run MAFFT using FFT-NS-2
    f_mafft(fn_fasta_concat, fn_fasta_concat_aligned, "--retree 2", exe_mafft)

    # run IQ-Tree 2
    f_iqtree2(fn_fasta_concat_aligned, exe_iqtree2)

    # check if treefile exists
    if (!file.exists(fn_fasta_concat_aligned_treefile)) {
        return(NULL)
    }

    return(fn_tree)
}

stopCluster(nwcl)

# combine all BUSCO trees
ls_trees <- paste(ls_trees, collapse=" ")
fn_alltrees <- paste0(dir_output_tree, "alltrees.tre")
system(paste("cat", ls_trees, ">", fn_alltrees))

# run ASTRAL-III
fn_astral_outfile <- paste0(dir_output_tree, "alltrees.astral.tre")
fn_astral_logfile <- paste0(dir_output_tree, "alltrees.astral.log")
f_astral(fn_alltrees, fn_astral_outfile, fn_astral_logfile, exe_astral)

################################################