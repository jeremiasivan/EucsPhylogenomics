# Generate target gene file for Hybpiper
dir_codes <- "/home/jeremias/EucsPhylogenomics/codes/"

fn_output <- "/data/jeremias/eucs/hybpiper/target_busco.fa"
dir_input_genes <- "/data/jeremias/eucs/busco/"

################################################

source(paste0(dir_codes, "/functions.R"))

################################################

ls_reference <- list.dirs(dir_input_genes, recursive=F, full.names=F)
ls_shared_gene <- c()

# iterate over references
for (ref in ls_reference) {
    # extract the list of genes
    ls_gene_sp <- list.files(paste0(dir_input_genes, "/", ref), pattern = "*.fa$", full.names = F, recursive = F)
    ls_gene_sp <- sapply(ls_gene_sp, function(x) { gsub(".fa", "", x) })
    
    # update the list of genes
    if (length(ls_shared_gene) == 0) {
        ls_shared_gene <- ls_gene_sp
    } else {
        ls_shared_gene <- intersect(ls_shared_gene, ls_gene_sp)
    }
}

# iterate over genes
for (gene in ls_shared_gene) {
    # iterate over references
    for (ref in ls_reference) {
        fn_fasta <- paste0(dir_input_genes, "/", ref, "/", gene, ".fa")

        # add sequence into one file
        f_fasta2msa(fn_fasta, paste0(ref, "-", gene), fn_output)
    }
}

################################################