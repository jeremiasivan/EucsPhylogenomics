################################################
#############     CALCULATE CF     #############
################################################

dir_output <- ""
thread <- 1

dir_fasta <- ""

exe_iqtree2 <- ""
exe_astral <- ""

################################################

library("ape")
library("doSNOW")

################################################

# extract all FASTA
ls_fasta <- list.files(dir_fasta, pattern="*.fa", recursive=F, full.names=F)
ls_locus <- sapply(ls_fasta, function(x) {
    gsub("*.fa", "", x)
})

# create doSNOW cluster
nwcl <- makeCluster(thread)
doSNOW::registerDoSNOW(nwcl)

# iterate over loci
foreach (locus = ls_locus) %dopar% {
    # set prefix
    fn_fasta <- paste0(dir_fasta, "/", locus, ".fa")
    prefix <- paste0(dir_output, "/", locus, "/", locus)

    # run IQ-TREE2
    f_iqtree2(fn_fasta, prefix, exe_iqtree2)
}

stopCluster(nwcl)

# combine all locus trees into one file
fn_all_trees <- paste0(dir_output, "/all_locus.treefile")
system(paste0("cat ", dir_output, "/*/*.treefile > ", fn_all_trees))

# generate species tree
fn_astral_log <- paste0(dir_output, "/astral.log")
fn_astral_tre <- paste0(dir_output, "/astral.treefile")
f_astral(fn_all_trees, fn_astral_tre, fn_astral_log, exe_astral)

# open the treefile
astral_tree <- ape::read.tree(fn_astral_tre)

# create doSNOW cluster
nwcl <- makeCluster(thread)
doSNOW::registerDoSNOW(nwcl)

# iterate over trees
foreach (locus = ls_locus) %dopar% {
    # input file
    fn_fasta <- paste0(dir_fasta, "/", locus, ".fa")

    # open locus treefile
    tree <- ape::read.tree(paste0(dir_output, "/", locus, "/", locus, ".treefile"))

    # subset the species tree
    astral_tree_sub <- ape::keep.tip(astral_tree, tip=tree$tip.label)

    # extract all branching events
    astral_tree_sub_branches <- ape::prop.part(astral_tree_sub)
    locus_tree_branches <- ape::prop.part(tree)

    # convert the number to species
    astral_tree_sub_branches <- f_part_num2chr(astral_tree_sub_branches)
    locus_tree_branches <- f_part_num2chr(locus_tree_branches)

    # extract incongruence
    incongruence_branches <- f_compare_parts(locus_tree_branches, astral_tree_sub_branches)

    # iterate over incongruence branches
    i <- 1
    for (br in incongruence_branches) {
        tree_tips <- strsplit(br, split="-")[[1]]
        if (length(tree_tips) < 4) {
            next
        }

        ref_topology <- ape::keep.tip(astral_tree, tip=tree_tips)
        ref_topology$node.label <- NULL
        ref_topology$edge.length <- NULL

        ref_topology <- ape::write.tree(ref_topology)

        # save the topology
        fn_output <- paste0(dir_output, "/", locus, "/constr", i, ".topfile")
        f_write_log(fn_output, msg=ref_topology)

        # run IQ-TREE2 with constrained topology
        prefix <- paste0(dir_output, "/", locus, "/constr", i)
        f_iqtree2_constrained(fn_fasta, fn_output, prefix, exe_iqtree2)

        # update i
        i <- i+1
    }

    # combine all trees
    fn_constr_alltrees <- paste0(dir_output, "/", locus, "/constr.alltrees")
    system(paste0("cat ", dir_output, "/", locus, "/constr*.treefile > ", fn_constr_alltrees))

    # run AU test
    prefix <- paste0(dir_output, "/", locus, "/au")
    f_iqtree2_au(fn_fasta, fn_constr_alltrees, prefix, exe_iqtree2)
}

stopCluster(nwcl)