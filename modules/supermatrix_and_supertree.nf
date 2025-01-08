#!/usr/env nextflow


// for each family tuple, it gets a list of species id and genes id. 
// Using these two information it will retain only the specified species files containing each only the speciefied gene sequences. 
process extract_fasta_per_species {
	label 'process_low'
        tag "${fam}"

        input:
	tuple val(fam), path(intersect), path(species), path(rest)

        output:
        tuple val(fam), file("${fam}.*_domain_sequences_after_intersection.fasta"), emit: selected_fasta

	script:
        """
        extract_intersecting_seq.pl ${intersect} ${species} \$PWD -after ${fam}
        """
}


// does exactlòy what the above does but puts everything into a single file and changes the header ids to be uniquea nd species specific.
process extract_fasta_per_species_for_full_aln {
	label 'process_low'
        tag "${fam}"

        input:
	tuple val(fam), path(intersect), path(species), path(rest)

        output:
        tuple val(fam), path("${intersect}"), path("${species}"), path("${fam}.domain_sequences_prior_after_intersection.fasta"), emit: selected_fasta

	script:
        """
        extract_intersecting_seq_prior.pl ${intersect} ${species} \$PWD -prior ${fam}
        cat ${fam}.*_domain_sequences_prior_after_intersection.fasta > ${fam}.domain_sequences_prior_after_intersection.fasta
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
/*
#split_seq_from_full_aln_based_on_species.pl ${intersect} ${species} ${fam}
        #for i in `find *_domain_sequences_prior_after_intersection.fasta_aln -maxdepth 0 -type f`; do
        #        t_coffee -other_pg seq_reformat -code ${code_file} -in \$i > \${i%".fasta_aln"}_coded.fasta_aln        
        #        t_coffee -other_pg seq_reformat -in \${i%".fasta_aln"}_coded.fasta_aln -output phylip_aln > \${i%".fasta_aln"}_coded.phylip 
        #done
*/
process run_full_alns {
	label 'cpu'
        tag "${fam} - ${mode}"
        publishDir "${fam}/results_without_mouse_10_samples/${mode}", mode: 'copy'

        input:
        tuple val(fam), path(intersect), path(species), path(fasta), path(code_file)
        val mode

        output:
        tuple val(fam), path("*.phylip"), emit: split_aln
        tuple val(fam), path("${output_alns_phy}"), emit: full_aln
        path ("full_aln.code_name")

        script:
        output_alns = fasta.baseName + "_full.fasta_aln"
        output_alns_coded = fasta.baseName + "_full_coded.fasta_aln"
        output_alns_phy = fasta.baseName + "_full_coded.phylip"

        """
        t_coffee ${fasta} -output fasta -outfile ${output_alns} -thread 2
        t_coffee -other_pg seq_reformat -in ${output_alns} -output code_name > full_aln.code_name
        t_coffee -other_pg seq_reformat -code full_aln.code_name -in ${output_alns} > ${output_alns_coded}
        t_coffee -other_pg seq_reformat -in ${output_alns_coded} -output phylip_aln > ${output_alns_phy}
        """
}

process run_phylo_full {
        label 'cpu'
        tag "${fam}"
        publishDir "${fam}/results_without_mouse_10_samples/${mode}", mode: 'copy'

        input:
        tuple val(fam), path(phylip_aln)
        val mode

        output:
        path "*"

        script:
        raxml_output = phylip_aln.baseName + "_raxml.nwk"
        """
        raxml -D -m PROTGAMMALG -s ${phylip_aln} -n ${raxml_output} -p 2233
        fastme -i ${phylip_aln} -p -g 1.0 -s -n
        """
}

process run_phylo_recon {
        label 'cpu'
        tag "${fam}"
        publishDir "${fam}/results_without_mouse_10_samples/${mode}", mode: 'copy'

        input:
        tuple val(fam), path(phylip_aln)
        val mode

        output:
        path "*"

        script:
        raxml_output = phylip_aln.baseName + "_raxml.nwk"
        """
        raxml -D -m PROTGAMMALG -s ${phylip_aln} -n ${raxml_output} -p 2233
        fastme -i ${phylip_aln} -p -g 1.0 -s -n
        """
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
