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
  
  if (length(args) < 3) {
    stop("Usage: script.R <BigTree file> <mapping file> <species file> <output file>", call. = FALSE)
  }
  
  list(
    bigtree_file = args[1],
    mapping_file = args[2],
    species_file = args[3],
    output_file = args[4]
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
    
    # Get the tips that match the species name
    tips_to_keep <- grep(paste0("^", species), full_aln_tree$tip.label, perl = TRUE)
    
    # Extract and unroot the subtree for this species
    test_subtree <- keep.tip(full_aln_tree, tips_to_keep)
    test_subtree <- unroot(test_subtree)
    
    # Rename the tips
    #test_subtree$tip.label <- gsub(paste0(species, "_"), "Seq_", test_subtree$tip.label)
    
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