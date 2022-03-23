#!/usr/bin/env Rscript

library("rlist")

argv = commandArgs(trailingOnly = TRUE)

if (length(argv)!=3)
{
  stop("Usage: Rscript intersect_orthologs.R <ORG_IDS> <FAMILY_ID> <NUM_OF_SPECIES> \n", call. = FALSE)
}

orgs<-read.table(file=argv[1])
mylist <- list()
for (org in levels(orgs$V1)) {
  filename <- paste(argv[2],".",org,"_domain_sequences_intersect_with_ref_after_OMA.fasta.headers",sep = "")
  if (file.size(filename) != 0) {
	myvar <- read.table(file = filename,header = FALSE,stringsAsFactors=FALSE)
  	mylist[[org]] = myvar$V1
  }
  
}

stored_intersect = c(0)
for (num_species in argv[3]:argv[3]) {
	all_comb = combn(names(mylist),num_species)
	for (i in 1:ncol(all_comb)) {
		print(paste0("Number of species: ",num_species,"  Combination ",i,"  out of ",ncol(all_comb)))
		test_list = list()
    	for (set in all_comb[,i]) {
      	  test_list = list.append(test_list,mylist[[set]])
    	}
    overlap_all <- Reduce(intersect, test_list)
    print(paste("Overlap length: ",length(overlap_all)))
    if (is.null(max(stored_intersect)) == FALSE) {
    if ( length(overlap_all) > max(stored_intersect) ) {
        save_species = all_comb[,i]
        save_genes = overlap_all
    }
    stored_intersect = c(stored_intersect,length(overlap_all))
    }
    }
}

write(save_genes,file = paste0(argv[2],".intersecting_genes"))
write(save_species,file = paste0(argv[2],".orthologs_org_ids_to_concatenate"))
