#!/usr/env nextflow

nextflow.enable.dsl=2

process extract_fasta_aln_per_species_sim {
	label 'cpu'
	tag "${fam}"	
	publishDir "${fam}/results", mode: 'copy'

	input:
	tuple val(fam), path(init_aln), path(orthologs_ids)

	output:
	tuple val(fam),path("${fam}*.fasta_aln"), emit: all_aln_sim
	path "*"

	script:

	"""
	extract_seq_per_species_sim.pl ${orthologs_ids} ${init_aln}
	for fasta in `find *.fasta_aln -maxdepth 0 -type f | grep -v "Guinea_Pig"`; do cp "${fam}".ma_Guinea_Pig.fasta_aln \$fasta; done
	"""

}

process concatenate_alns_sim {
	label 'cpu'
	tag "${fam}"
	publishDir "${fam}/results", mode: 'copy'

	input:
	tuple val(fam), val(all_aln)
	val ortho

	output:
	path "*"

	shell:
	'''
	concatenate_aln.py -c -m B -f !{fam} !{all_aln.join(" ")}
	for sample_num in $(eval echo {1..10}); do
			for single_unit_num in $(eval echo {1..!{ortho}}); do
				eval "cat single_unit_{1..$single_unit_num}_sample_${sample_num}_concatenated_aln.phylip_fastme_tree.nwk > unit_${single_unit_num}_sample_${sample_num}_concatenated_aln.phylip_fastme_tree_for_supertree.nwk"
				/SuperFine/runSuperFine.py -r gmrp unit_${single_unit_num}_sample_${sample_num}_concatenated_aln.phylip_fastme_tree_for_supertree.nwk > unit_${single_unit_num}_sample_${sample_num}_concatenated_aln.phylip_fastme_superfine_tree.nwk
				eval "cat RAxML_bestTree.single_unit_{1..$single_unit_num}_sample_${sample_num}_concatenated_aln_raxml.nwk > unit_${single_unit_num}_sample_${sample_num}_concatenated_aln.raxml_tree_for_supertree.nwk"
				/SuperFine/runSuperFine.py -r gmrp unit_${single_unit_num}_sample_${sample_num}_concatenated_aln.raxml_tree_for_supertree.nwk > unit_${single_unit_num}_sample_${sample_num}_concatenated_aln.raxml_superfine_tree.nwk
			done
	done
	'''
}
