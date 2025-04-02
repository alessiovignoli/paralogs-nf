#!/usr/env nextflow


// for each family tuple, it gets a list of species id and genes id. 
// Using these two information it will retain only the specified species files containing each only the speciefied gene sequences. 
process extract_fasta_per_species {
	label 'process_low'
        tag "${fam}"
        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
                'docker://athbaltzis/paralogs:v0.10' :
                'athbaltzis/paralogs:v0.10' }"


        input:
	tuple val(fam), path(intersect), path(species), path(rest)

        output:
        tuple val(fam), file("${fam}.*_domain_sequences_after_intersection.fasta"), emit: selected_fasta

	script:
        """
        extract_intersecting_seq.pl ${intersect} ${species} \$PWD -after ${fam}
        """
}


// does exactlòy what the above does but puts everything into a single file and changes the header ids to be uniquea and species specific.
process extract_fasta_per_species_for_full_aln {
	label 'process_low'
        tag "${fam}"
        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
                'docker://athbaltzis/paralogs:v0.10' :
                'athbaltzis/paralogs:v0.10' }"

        input:
	tuple val(fam), path(intersect), path(species), path(rest)

        output:
        tuple val(fam), path("${fam}.domain_sequences_prior_after_intersection.fasta"), emit: selected_fasta

	script:
        """
        extract_intersecting_seq_prior.pl ${intersect} ${species} \$PWD -prior ${fam}
        cat ${fam}.*_domain_sequences_prior_after_intersection.fasta > ${fam}.domain_sequences_prior_after_intersection.fasta
        """
}

// takes a MSA and a list of species names, -> one msa in fasta format per species, each with all sequences in input fasta from the same species. Order of sequences is preserved.
process extract_fasta_aln_per_species_emp {
	label 'process_low'
	tag "${fam}"
        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
                'docker://athbaltzis/paralogs:v0.10' :
                'athbaltzis/paralogs:v0.10' }"

	input:
	tuple val(fam), path(init_aln), path(orthologs_ids)

	output:
	tuple val(fam), path("${fam}*.fasta_aln"), emit: species_aln

	script:
        prefix = init_aln.baseName
	"""
	# Loop through each ID in the ID file
	while read -r ID; do
		# Create the output file for the current ID
		OUTPUT_FILE="${prefix}_\${ID}.fasta_aln"

		# Use awk to extract the headers containing the ID and their sequences 
                # remove the id (species) from the header and prepend it with "Seq_"
		awk -v id="\$ID"  -v header_prefix="\$OUTPUT_FILE" '
		BEGIN { found=0; }
		/^>/ {
			if (index(\$0, id) > 0) { 
				found=1;
                                gsub(id, "", \$0);
                                gsub("__", "_", \$0);
                                gsub(">", ">Seq_", \$0);
				print \$0 > header_prefix;
			} else {
				found=0;
			}
		}
		!/^>/ {
			if (found) { 
				print \$0 > header_prefix;
			}
		}
		' "${init_aln}"

	done < ${orthologs_ids}
	"""
}

// depending on the alligner asked (aka mode). it executes a different template file present in the template dir.
// but all in all it alligns and renames the fasta in input. 
process run_alns {
	label 'process_low'
        tag "${fam} - ${mode}"
        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
                'docker://athbaltzis/paralogs:v0.10' :
                'athbaltzis/paralogs:v0.10' }"
        

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


// creates the bigtree phylip MSA, using all paralogs ortholog informations. It creates a header id code based om code_file.
process run_full_alns {
	label 'process_medium'
        tag "${fam} - ${params.mode}"
        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
                'docker://athbaltzis/paralogs:v0.10' :
                'athbaltzis/paralogs:v0.10' }"

        input:
        tuple val(fam), path(fasta)

        output:
        tuple val(fam), path("${output_alns_phy}"), emit: full_aln
        tuple val(fam), path("full_aln.code_name"), emit: msa_code_name
        tuple val(fam), path("${output_alns}"), emit: fasta_full_aln

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


// it simply removes the reference species keyword fromn the list of species, while returning the same tuple.
process remove_ref_species {
        label 'process_low'
        tag "${fam}"
        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
                'docker://athbaltzis/paralogs:v0.10' :
                'athbaltzis/paralogs:v0.10' }"

        input:
	tuple val(fam), path(intersect), path(species), path(rest)

        output:
        tuple val(fam), path(intersect), path("*_no_ref"), path(rest), emit: no_ref
        tuple val(fam), path("*_no_ref"), emit: modified_ortho

        script:
        """
        sed "/${params.ref}/d" ${species} > ${species}_no_ref
        """
}


// The following process will do the whole SuperMatrix and SuperTree computation. With 10 replication as well.
// concatenate_aln.py  will do the following:
// create all the MSA units and single units in phylip format. From the starting species alignments by shuffling and drowing columns.
//      unit = divide the concatenated shuffled MSA in 6 blocks, block1 is unit1, block1 + block2 is unit2 and so on.
//      single_unit = divide the concatenated shuffled MSA in 6 blocks, each block is is own single unit.
// this above process is repeated 10 times, so to generate statistical replicates called sample.
//      so will have unit1_sample1 that will be different from unit1_sample2 and so on.
// For all the above MSA (6 species/units 10 replicates/samples) -> 60 MSA unit, 60 MSa single_unit
// 	    for each unit MSA make raxml (ML) (60) and fastme (ME) (60 trees -> this are the SuperMatrix ML and Me final trees (120)
//      for each single_unit MSA make raxml (ML) (60) and fastme (ME) (60) trees
// the bash for loop will then take the (60) single_unit trees and append them in a file in blocks like: single_unit1 tree  (file1), single_unit1 tree  + single_unit2 tree (file2) ecc.. . this is done for both ME and ML trees respectively.
//      this (120) "appended" file are fed to superfine that will do the SuperTree approach, this are the final SuperTree trees.
process concatenate_alns {
	label 'process_high'
        tag "${fam}"
        container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
                'docker://athbaltzis/paralogs:v0.10' :
                'athbaltzis/paralogs:v0.10' }"

        input:
        tuple val (fam), val(all_aln)
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
