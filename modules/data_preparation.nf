#!/usr/env nextflow

nextflow.enable.dsl=2

process data_preparation {
	tag "${fam}"
	publishDir "${fam}", mode: 'copy'	
	errorStrategy 'ignore'

	input:
	tuple val(fam), path(fasta)
        val ref
        val species_num

	output:
	path "*" optional true into myout
	tuple val(fam), path("${fam}.intersecting_genes"), path("${fam}.*_domain_sequences_intersect_with_ref_after_OMA.fasta") optional true into prep_out

	script:
	"""
        echo ${ref} > ref_org_id
	extract_fasta_for_org_from_PFAM_fasta.pl ${fasta} ref_org_id
	isolate_UNIPROT_ID.bash ${fam}'.'${ref}'_domain_sequences.fasta' > ${fam}'.'${ref}'_domain_sequences.fasta.headers'
	find_1_to_1_orthologs_for_ref_org.pl ${fam}"."${ref}"_domain_sequences.fasta.headers" | sort -V | uniq > ${fam}"."${ref}"_domain_sequences.uniprot_oma-groups"
        count_species_freq_in_oma_groups.bash ${fam}"."${ref}"_domain_sequences.uniprot_oma-groups" | head -30 | awk '{print \$2}' | sed "s|PIGXX|PIG|g" | sed "s|RATNO|RAT|g" > selected_orthologs
	extract_fasta_for_org_from_PFAM_fasta.pl ${fasta} selected_orthologs
        grep -v ${ref} selected_orthologs > selected_ortho
        cut -f1 ${fam}"."${ref}"_domain_sequences.uniprot_oma-groups" > ${fam}"."${ref}"_domain_sequences_after_oma.uniprot"
	awk '{print \$1"\t"\$1}' ${fam}"."${ref}"_domain_sequences_after_oma.uniprot" > ${fam}"."${ref}"_domain_sequences_after_oma.uniprot_map"
	map_genenames_to_fasta.pl ${fam}"."${ref}"_domain_sequences_after_oma.uniprot_map" ${fam}"."${ref}"_domain_sequences.fasta" > ${fam}"."${ref}"_domain_sequences_intersect_with_ref_after_OMA.fasta"
	grep ">" ${fam}"."${ref}"_domain_sequences_intersect_with_ref_after_OMA.fasta" | tr -d ">" | tr "/" "\t" | cut -f1 > ${fam}"."${ref}"_domain_sequences_intersect_with_ref_after_OMA.fasta.headers"
	if [ `cat ${fam}"."${ref}"_domain_sequences_intersect_with_ref_after_OMA.fasta.headers" | wc -l` -gt 10 ]; then
        	extract_1_to_1_orthologs_for_each_species_OMA.pl ${fam}"."${ref}"_domain_sequences.uniprot_oma-groups" selected_ortho ${fam}"."${ref}"_domain_sequences_after_oma.uniprot_map" ${fam} ${ref}
        	map_genenames_to_fasta.bash selected_ortho ${fam}
        	intersect_orthologs.R selected_orthologs ${fam} ${species_num}
        	if [[ ! -f "${fam}".intersecting_genes ]] || [[ `cat "${fam}".intersecting_genes | wc -l` -lt 10 ]]; then rm -r *; fi
	else
        	rm -r *
	fi
	"""
}
