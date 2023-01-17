#!/usr/env nextflow

nextflow.enable.dsl=2

process extract_fasta_per_species {
	label 'cpu'
        tag "${fam}"
        publishDir "${fam}/results_without_mouse_10_samples/${mode}", mode: 'copy'

        input:
	tuple val(fam), path(intersect), path(species), path(rest)
        val mode

        output:
        tuple val(fam), file("${fam}.*_domain_sequences_after_intersection.fasta")

	script:
        """

        extract_intersecting_seq.pl ${intersect} ${species} \$PWD -after ${fam}

        """
}

process run_alns {
	label 'cpu'
        tag "${fam} - ${mode}"
        publishDir "${fam}/results_without_mouse_10_samples/${mode}", mode: 'copy'

        input:
        tuple val(fam), path(fasta)
        val mode

        output:
        tuple val(fam), path(output_alns), emit: aln_output
        tuple val(fam), path("file.code_name"), emit: code_name

        script:
        output_alns = fasta.baseName + "_coded.fasta_aln"
        output_ph = fasta.baseName + "_coded.phylip"

        if (mode == 'tcoffee')
                template 'tcoffee.sh'
        else if (mode == 'psicoffee')
                template 'psicoffee.sh'
        else if (mode == 'ginsi')
                template 'ginsi.sh'
        else
                error "Invalid alignment mode: ${mode}"
}

process concatenate_alns {
	label 'cpu'
        tag "${fam}"
        publishDir "${fam}/results_without_mouse_10_samples/${mode}", mode: 'copy'

        input:
        tuple val (fam), val(all_aln)
        val mode
        val ortho

        output:
        path "*"

        shell:
	'''
	concatenate_aln.py -c -m B -f !{fam} !{all_aln.join(" ")}
	for sample_num in $(eval echo {1..10}); do
			for single_unit_num in $(eval echo {1..6}); do
				eval "cat single_unit_{1..$single_unit_num}_sample_${sample_num}_concatenated_aln.phylip_fastme_tree.nwk > unit_${single_unit_num}_sample_${sample_num}_concatenated_aln.phylip_fastme_tree_for_supertree.nwk"
				/SuperFine/runSuperFine.py -r gmrp unit_${single_unit_num}_sample_${sample_num}_concatenated_aln.phylip_fastme_tree_for_supertree.nwk > unit_${single_unit_num}_sample_${sample_num}_concatenated_aln.phylip_fastme_superfine_tree.nwk
				eval "cat RAxML_bestTree.single_unit_{1..$single_unit_num}_sample_${sample_num}_concatenated_aln_raxml.nwk > unit_${single_unit_num}_sample_${sample_num}_concatenated_aln.raxml_tree_for_supertree.nwk"
				/SuperFine/runSuperFine.py -r gmrp unit_${single_unit_num}_sample_${sample_num}_concatenated_aln.raxml_tree_for_supertree.nwk > unit_${single_unit_num}_sample_${sample_num}_concatenated_aln.raxml_superfine_tree.nwk
			done
	done
	'''
}
