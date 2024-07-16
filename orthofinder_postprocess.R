# Convert OrthoFinder ID to Species Names
thread <- 10

fn_refseq <- "/EucsPhylogenomics/eucs_refseq.txt"
fn_orthogroup <- "/EucsPhylogenomics/orthofinder/orthogroups/Orthogroups.tsv"
fn_species_tree <- "/EucsPhylogenomics/orthofinder/SpeciesTree_rooted.txt"

prefix_dir_fasta <- "/EucsPhylogenomics/orthofinder/single_copy_orthologue"
prefix_dir_tree <- "/EucsPhylogenomics/orthofinder/gene_tree"

################################################

# load libraries
# require(data.table)
library(doSNOW)

# open data.table
df_refseq <- data.table::fread(fn_refseq)
df_orthogroup <- data.table::fread(fn_orthogroup)

ls_refseq <- unique(df_refseq$id)

# ---------- Species Tree ----------
if (fn_species_tree != "") {
    tre <- ape::read.tree(fn_species_tree)

    for (tips in tre$tip.label) {
        tre$tip.label <- gsub(tips, df_refseq$species[df_refseq$id==tips], tre$tip.label)
    }

    ape::write.tree(tre, file=paste0(fn_species_tree, ".species"))
}

# ----------- Gene Trees -----------
if (prefix_dir_tree != "") {
    # create doSNOW cluster
    nwcl <- makeCluster(thread)
    doSNOW::registerDoSNOW(nwcl)

    # extract the list of gene trees
    ls_gene_trees <- list.files(prefix_dir_tree, full.names=F, recursive=F)

    # iterate over gene trees
    foreach (genetree=ls_gene_trees) %dopar% {
        tre <- ape::read.tree(paste0(prefix_dir_tree,"/",genetree))

        gene <- gsub("_tree.txt","",genetree)

        for (ref in ls_refseq) {
            ref_updated <- gsub("\\.","_",ref)
            ref_suffix <- df_orthogroup[[ref]][df_orthogroup$Orthogroup==gene]
            ref_updated <- paste0(ref_updated,"_",ref_suffix)

            tre$tip.label <- gsub(ref_updated, df_refseq$species[df_refseq$id==ref], tre$tip.label)
        }

        ape::write.tree(tre, file=paste0(prefix_dir_tree,"/",genetree,".species"))
    }

    stopCluster(nwcl)
}

# -------- FASTA Alignments --------
if (prefix_dir_fasta != "") {
    # create doSNOW cluster
    nwcl <- makeCluster(thread)
    doSNOW::registerDoSNOW(nwcl)

    # extract the list of FASTAs
    ls_genes <- list.files(prefix_dir_fasta, full.names=F, recursive=F)

    # iterate over FASTAs
    foreach (gene=ls_genes) %dopar% {
        aln <- ape::read.FASTA(paste0(prefix_dir_fasta,"/",gene))

        gene_updated <- gsub(".fa","",gene)

        for (ref in ls_refseq) {
            ref_suffix <- df_orthogroup[[ref]][df_orthogroup$Orthogroup==gene_updated]

            names(aln) <- gsub(ref_suffix, df_refseq$species[df_refseq$id==ref], names(aln))
        }

        ape::write.FASTA(tre, file=paste0(prefix_dir_fasta,"/",gene,".species"))
    }

    stopCluster(nwcl)
}
