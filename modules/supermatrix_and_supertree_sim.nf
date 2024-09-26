#!/usr/env nextflow

nextflow.enable.dsl=2

process extract_fasta_aln_per_species_sim {
	label 'cpu'
	tag "${fam}"

	input:
	tuple val(fam), path(init_aln), path(orthologs_ids)

	output:
	tuple val(fam),path("${fam}*.fasta_aln"), emit: all_aln_sim
	path "*"

	script:

	"""
	extract_seq_per_species_sim.pl ${orthologs_ids} ${init_aln} -after
	"""

}

// takes a MSA and a list of species names, -> phylip format MSA with all sequence from the same species one after the other and header names changed into code values.
process extract_fasta_aln_per_species_sim_for_full_aln {
	label 'process_low'
	tag "${fam}"

	input:
	tuple val(fam), path(init_aln), path(orthologs_ids)

	output:
	tuple val(fam),path("${fam}_full_aln_coded.phylip"), emit: phylip_full_aln_sim
	//path "*"

	script:

	"""
	extract_seq_per_species_sim.pl ${orthologs_ids} ${init_aln} -prior
	cat *_for_full_aln.fasta_aln > ${fam}_full_aln.fasta_aln
	t_coffee -other_pg seq_reformat -in ${fam}_full_aln.fasta_aln -output code_name > full_aln.code_name
    t_coffee -other_pg seq_reformat -code full_aln.code_name -in ${fam}_full_aln.fasta_aln > ${fam}_full_aln_coded.fasta_aln
    t_coffee -other_pg seq_reformat -in ${fam}_full_aln_coded.fasta_aln -output phylip_aln > ${fam}_full_aln_coded.phylip
	"""

}

// computes the ML tree given an alignment in phylip form.
process run_phylo_ML_full_sim {
        label 'process_low'
        tag "${fam}"

        input:
        tuple val(fam), path(phylip_aln)

        output:
        path "*"

        script:
        raxml_output = phylip_aln.baseName + "_raxml.nwk"
        """
        time raxml -D -m PROTGAMMALG -s ${phylip_aln} -n ${raxml_output} -p 2233
        """
}

// computes the ME tree given an alignment in phylip form.
process run_phylo_ME_full_sim {
        label 'process_low'
        tag "${fam}"

        input:
        tuple val(fam), path(phylip_aln)

        output:
        path "*"

        script:
        """
        fastme -i ${phylip_aln} -p -g 1.0 -s -n
        """
}

process only_concatenate_aln_sim {
	label 'cpu'
	tag "${fam}"

	input:
	tuple val(fam), val(all_aln)

	output:
	tuple val(fam),path("${fam}_supermatrix.phylip"), emit: phylip_only_concatenate_aln_sim

	shell:
	'''
	concatenate_aln.py -f !{fam} -S !{all_aln.join(" ")}
	'''
}

process concatenate_alns_sim {
	label 'cpu'
	tag "${fam}"

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
				gotree compute consensus -i unit_${single_unit_num}_sample_${sample_num}_concatenated_aln.phylip_fastme_tree_for_supertree.nwk -o unit_${single_unit_num}_sample_${sample_num}_concatenated_aln.phylip_fastme_consensus_tree.nwk -f 0.5
				eval "cat RAxML_bestTree.single_unit_{1..$single_unit_num}_sample_${sample_num}_concatenated_aln_raxml.nwk > unit_${single_unit_num}_sample_${sample_num}_concatenated_aln.raxml_tree_for_supertree.nwk"
				/SuperFine/runSuperFine.py -r gmrp unit_${single_unit_num}_sample_${sample_num}_concatenated_aln.raxml_tree_for_supertree.nwk > unit_${single_unit_num}_sample_${sample_num}_concatenated_aln.raxml_superfine_tree.nwk
				gotree compute consensus -i unit_${single_unit_num}_sample_${sample_num}_concatenated_aln.raxml_tree_for_supertree.nwk -o unit_${single_unit_num}_sample_${sample_num}_concatenated_aln.raxml_consensus_tree.nwk -f 0.5
			done
	done
	'''
}
