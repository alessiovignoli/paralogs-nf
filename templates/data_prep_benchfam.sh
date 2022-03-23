perl ${params.path_to_bin}/extract_fasta_for_org_from_BENCHFAM_fasta.pl ${fam}/${fam.baseName}"_pfam.fa" ${ortho_ids} ${fam}/${fam.baseName}"_ref.fa"
cd ${params.ref}
bash ${params.path_to_bin}/isolate_UNIPROT_ID.bash ${params.ref}"_"${fam.baseName}"_domain_sequences.fasta" > ${params.ref}"_"${fam.baseName}"_domain_sequences.fasta.headers"
perl ${params.path_to_bin}/find_1_to_1_orthologs_for_ref_org.pl ${params.ref}"_"${fam.baseName}"_domain_sequences.fasta.headers" | sort -V | uniq > ${params.ref}"_"${fam.baseName}"_domain_sequences.uniprot_oma-groups"
cut -f1 ${params.ref}"_"${fam.baseName}"_domain_sequences.uniprot_oma-groups" > ${params.ref}"_domain_sequences_after_oma.uniprot"
perl ${params.path_to_bin}/download_uniprot_acc_id_to_gene_name.pl ${params.ref}"_domain_sequences_after_oma.uniprot" > ${params.ref}"_domain_sequences_after_oma.uniprot_genenames"
perl ${params.path_to_bin}/map_genenames_to_fasta.pl ${params.ref}"_domain_sequences_after_oma.uniprot_genenames" ${params.ref}"_"${fam.baseName}"_domain_sequences.fasta" > ${params.ref}"_domain_sequences_intersect_with_ref_after_OMA.fasta"
grep ">" ${params.ref}"_domain_sequences_intersect_with_ref_after_OMA.fasta" | tr -d ">" | tr "/" "\t" | cut -f1 > ${params.ref}"_domain_sequences_intersect_with_ref_after_OMA.fasta.headers"
if [ `cat ${params.ref}"_domain_sequences_intersect_with_ref_after_OMA.fasta.headers" | wc -l` -gt 15 ]; then
	cd ../
        perl ${params.path_to_bin}/extract_1_to_1_orthologs_for_each_species_OMA.pl ${params.ref}/${params.ref}"_"${fam.baseName}"_domain_sequences.uniprot_oma-groups" ${ortho} ${params.ref}/${params.ref}"_domain_sequences_after_oma.uniprot_genenames"
        ${params.path_to_bin}/map_genenames_to_fasta.bash ${ortho} ${fam.baseName}
        Rscript ${params.path_to_bin}/intersect_orthologs_auto.R ${ortho_ids}
else
	cd ../
	rm -r *
fi

