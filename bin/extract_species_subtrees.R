#!/usr/bin/env Rscript

# this is made to make it run with other things that are not Docker and have R still work 
.libPaths(c("/usr/local/lib/R/site-library", "/usr/local/lib/R/library"))
print(.libPaths)

# Load necessary libraries
library(ape)
library(phangorn)
library(phytools)

# Function to parse command-line arguments
parse_args <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  if (length(args) < 5) {
    stop("Usage: script.R <BigTree file> <mapping file> <species file> <output file> <mode>", call. = FALSE)
  }
  
  mode <- args[5]
  if (!(mode %in% c("emp", "sim"))) {
    stop("<mode> must be either 'emp' or 'sim'", call. = FALSE)
  }
  
  list(
    bigtree_file = args[1],
    mapping_file = args[2],
    species_file = args[3],
    output_file = args[4],
    mode = mode
  )
}

# Main function
main <- function() {
  
  # Parse command-line arguments
  args <- parse_args()
  
  bigtree_file <- args$bigtree_file
  mapping_file <- args$mapping_file
  species_file <- args$species_file
  output_file <- args$output_file
  mode <- args$mode
  
  # Read the BigTree tree and unroot it
  full_aln_tree <- read.tree(bigtree_file)
  full_aln_tree <- unroot(full_aln_tree)
  
  # Read the sequence header mapping
  full_aln_codefile <- read.table(file = mapping_file, header = FALSE)
  
  # Read the names of the ortholog species
  orgs <- read.table(file = species_file, header = FALSE)
  
  # Prepare to collect the results
  all_species_raxml_subtrees_from_full_aln <- list()

  # Rename each tip label in the tree according to the mapping
  for (j in 1:length(full_aln_tree$tip.label)) { 
    full_aln_tree$tip.label[j] <- as.character(full_aln_codefile$V1[which(full_aln_codefile$V2 == full_aln_tree$tip.label[j])])
  }
  
  # Iterate over species and extract paralog subtrees
  for (species in orgs$V1) {
    
    # Get the tips that match the species name. two different ways of doing that based on the mode.
    # in sim mode the mapping is species_<num>. while in emp mode the mapping key is <gene_name>_species_<num>.
    # in one case is enought to rely on presence while in the other it needs to be more stringent.
    tips_to_keep <- NULL
    if (mode == "emp") {
      tips_to_keep <- grep(species, full_aln_tree$tip.label, perl = TRUE)
    } else {
      tips_to_keep <- grep(paste0("^", species), full_aln_tree$tip.label, perl = TRUE)
    }
    
    # Extract and unroot the subtree for this species
    test_subtree <- keep.tip(full_aln_tree, tips_to_keep)
    test_subtree <- unroot(test_subtree)
    
    print(sapply(strsplit(test_subtree$tip.label, split="_"), "[", 1))

    # Rename the tips. once again a different schema depending on the mode.
    # in sim mode the tip label will be in the form of Seq_<num>.
    # in emp mode the tip label will be in the form of Seq_<gene_name>.
    if (mode == "emp") {
      test_subtree$tip.label <- paste0("Seq_", sapply(strsplit(test_subtree$tip.label, split="_"), "[", 1))
    } else {
      test_subtree$tip.label <- gsub(paste0(species, "_"), "Seq_", test_subtree$tip.label)
    }
    
    # Change class of tree and add it to the collection
    test_subtree <- as.multiPhylo(test_subtree)
    all_species_raxml_subtrees_from_full_aln <- c(all_species_raxml_subtrees_from_full_aln, test_subtree)
  }
  
  # Set the class of the resulting collection
  class(all_species_raxml_subtrees_from_full_aln) <- "multiPhylo"
  
  # Write the subtrees to an output file in Newick format
  write.tree(all_species_raxml_subtrees_from_full_aln, file = output_file)
  
  cat("Subtrees successfully written to", output_file, "\n")
}

# Run the main function
main()