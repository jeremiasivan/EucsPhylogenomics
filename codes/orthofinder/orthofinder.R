################################################
#############      ORTHOFINDER     #############
################################################

dir_output <- "/data/jeremias/eucs/orthofinder/"
thread <- 10

# extract genes using GTF/GFF from assemblies
is_extract_gene <- FALSE
dir_fasta <- "/data/jeremias/eucs/genes/assemblies/"
dir_gtf <- "/data/jeremias/eucs/genes/gtf/"

exe_gffread <- "gffread"
fn_refseq_chr <- "eucs_refseq_chr.txt" # optional

# run Orthofinder
exe_orthofinder <- "Orthofinder"
exe_astral <- "astral"

dir_genes <- "/data/jeremias/eucs/genes/proteome/"

################################################

# require(data.table)
# require(tidyr)
library(doSNOW)

source(paste0(dir_codes, "/functions.R"))

# create output directory
if (!dir.exists(dir_output)) {
    dir.create(dir_output)
}

################################################

if (is_extract_gene) {
    # list all reference sequences
    ls_refseq <- list.files(dir_fasta, pattern="*.fa$", recursive=F, full.names=F)

    # create doSNOW cluster
    nwcl <- makeCluster(thread)
    doSNOW::registerDoSNOW(nwcl)

    # output directory
    dir_genes <- paste0(dir_output, "/genes/")
    if (!dir.exists(dir_genes)) {
        dir.create(dir_genes, recursive=T)
    }

    # iterate over reference sequences
    foreach (ref = ls_refseq) %dopar% {
        # update the input files
        fn_fasta <- paste0(dir_fasta, "/", ref, ".fa")
        fn_gtf <- paste0(dir_gtf, "/", ref, ".gtf")

        # output file
        fn_outfile <- paste0(dir_genes, "/", ref, ".faa")
        
        # read metadata file if provided
        if (fn_refseq_chr != "" && file.exists(fn_refseq_chr)) {
            # read GTF file
            df_gtf <- data.table::fread(fn_gtf)
            data.table::setnames(df_gtf, c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"))

            # open file
            df_refseq_chr <- data.table::fread(fn_refseq_chr)

            # convert GTF file to match naming convention
            df_refseq_chr_long <- tidyr::pivot_longer(df_refseq_chr, cols=-chr, names_to="species", values_to="assembly")
            df_refseq_chr_subset <- df_refseq_chr_long[df_refseq_chr_long$species==ref,]

            # iterate over chromosomes
            for (j in 1:nrow(df_refseq_chr_subset)) {
                # update the name of the chromosome following the assembly
                df_gtf$seqname <- gsub(df_refseq_chr_subset$chr[j], df_refseq_chr_subset$assembly[j], df_gtf$seqname)
            }

            # save the updated GTF file
            data.table::fwrite(df_gtf, file=fn_gtf, sep="\t", quote=F, col.names=F)
        }
        
        # run Gffread
        cmd_gffread <- paste(exe_gffread, "-g", fn_fasta, "-V -y", fn_outfile, fn_gtf)
        system(cmd_gffread)
    }

    stopCluster(nwcl)
}

################################################

# run Orthofinder
f_orthofinder(dir_genes, dir_output, thread, exe_orthofinder)

################################################

# extract the result folder
ls_dirs <- list.dirs(dir_output, recursive=F, full.names=T)
dir_output_run <- ls_dirs[grepl(".*Results*", ls_dirs)]
if (length(dir_output_run) != 1) {
    stop("Multiple Orthofinder output folders. Exited.")
}

# input files and directories
fn_orthogroup <- paste(dir_output_run, "/Orthogroups/Orthogroups.tsv")
fn_species_tree <- paste0(dir_output_run, "/Species_Tree/SpeciesTree_rooted.txt")

dir_output_run_alignment <- paste(dir_output_run, "/Single_Copy_Orthologue_Sequences/")
dir_output_run_genetrees <- paste(dir_output_run, "/Gene_Trees/")

# output directories
dir_sco <- paste0(dir_output, "/summary/single_copy_orthologues/")
dir_sco_fa <- paste0(dir_sco, "fasta/")
dir_sco_tre <- paste0(dir_sco, "trees/")
lapply(list(dir_sco_fa, dir_sco_tre), function(x){ if(!dir.exists(x)) dir.create(x, recursive=T) })

# open species tree
species_tre <- ape::read.tree(fn_species_tree)
ls_refseq <- species_tre$tip.label

# open data.table
df_orthogroup <- data.table::fread(fn_orthogroup)

################################################

# create doSNOW cluster
nwcl <- makeCluster(thread)
doSNOW::registerDoSNOW(nwcl)

# extract the list of FASTAs
ls_genes <- list.files(dir_output_run_alignment, full.names=F, recursive=F)

# iterate over FASTAs
ls_trees <- foreach (genealn = ls_genes, .combine='c') %dopar% {
    aln <- readLines(paste0(dir_output_run_alignment, "/", genealn))
    gene <- gsub(".fa", "", genealn)

    # iterate over reference sequences
    for (ref in ls_refseq) {
        ref_suffix <- df_orthogroup[[ref]][df_orthogroup$Orthogroup==gene]

        aln <- gsub(ref_suffix, ref, aln)
    }

    # save the alignment file
    writeLines(aln, con=paste0(dir_sco_fa, genealn))

    # copy the treefile
    fn_input_gene_tree <- paste0(dir_output_run_genetrees, gene, "_tree.txt")
    fn_output_gene_tree <- paste0(dir_sco_tre, gene, ".treefile")
    file.copy(fn_input_gene_tree, fn_output_gene_tree)

    return(fn_output_gene_tree)
}

stopCluster(nwcl)

# combine all gene trees
ls_trees <- paste(ls_trees, collapse=" ")
fn_alltrees <- paste0(dir_sco, "alltrees.tre")
system(paste("cat", ls_trees, ">", fn_alltrees))

# run ASTRAL-III
fn_astral_outfile <- paste0(dir_sco, "alltrees.astral.tre")
fn_astral_logfile <- paste0(dir_sco, "alltrees.astral.log")
f_astral(fn_alltrees, fn_astral_outfile, fn_astral_logfile, exe_astral)

################################################