#!/bin/bash

# Usage map_genenames_to_fasta.bash <ORG_IDS> <PFAM_ID>

for file in `cat $1`; do
	filename=$2"."$file"_domain_sequences_intersect_with_ref.tab"
	fasta=$2"."$file"_domain_sequences.fasta"
	output=$2"."$file"_domain_sequences_intersect_with_ref_after_OMA.fasta"
	headers=$output".headers"
	perl /users/cn/abaltzis/projects/1_paralogs/paraconc-nf/paraconc-nf_v0.2/bin/map_genenames_to_fasta.pl $filename $fasta > $output
	grep ">" $output | tr -d ">" | tr "/" "\t" | cut -f1 > $headers 
done
