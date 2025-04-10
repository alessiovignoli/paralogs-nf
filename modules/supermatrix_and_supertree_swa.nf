#!/usr/env nextflow

// takes a set of fasta files from one family and selects randomly a given number of files. in each of those files swap a given random numbers of sequences between themselves.
process swap_seq_in_fam {
	label 'process_low'
	tag "${fam}"
	container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://athbaltzis/paralogs:v0.10' :
        'athbaltzis/paralogs:v0.10' }"


	input:
	tuple val(fam), val(num_species_swap), val(num_swap), path(all_aln)

	output:
	// also inputs file have to be outputted since the swaaped once are chosen at random as subset
	tuple val(fam), val(num_species_swap), val(num_swap), path("*fasta_aln*", includeInputs : true), emit: swapped_fastas

	script:
	"""
	# the simlynk to input file that has been modified needs to go 
	# so that the output can be a set of (25) files (number of species) 
	for i in \$(ls | sort -R | head -n ${num_species_swap}); do
		fasta_sequence_swapper.py -i \$i \\
			-o \$i.${num_species_swap}_${num_swap} \\
			-n ${params.paralog_num_sim} \\
			-s ${num_swap}
		rm \$i
	done
	"""
}

// transforms a fasta MSA into a phylip format MSA. does this at family level through a for loop.
process from_fasta_to_phylip_swa {
    label 'process_low'
	tag "${id}"
	container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/t-coffee:13.46.0.919e8c6b--hfc96bf3_0':
        'quay.io/biocontainers/t-coffee:13.46.0.919e8c6b--hfc96bf3_0' }"


	input:
	tuple val(id), path(all_msa_fasta)

	output:
	tuple val(id), path("*.phylip"), emit: phylip

	script:
	"""
	export TEMP='./'

	# output name is created so all extentions after fasta_aln are replced with phylip instead
	for i in \$(echo ${all_msa_fasta}); do
		t_coffee -other_pg seq_reformat -in \$i -output phylip_aln > \${i%fasta_aln*}phylip
	done
	"""
}


// computes the ML tree given an alignment in phylip form. does this at family level through a for loop.
process run_phylo_ML_supertree_aln_swa {
	label 'process_low'
	tag "${fam}"
	container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://athbaltzis/paralogs:v0.10' :
        'athbaltzis/paralogs:v0.10' }"


	input:
	tuple val(fam), path(all_phylip_aln)

	output:
	path "*"
	tuple val(fam), path("RAxML_bestTree*"), emit: ml_bestree

	script:
	"""
	# extention .phylip is replaced with _raxml.nwk
	for i in \$(echo ${all_phylip_aln}); do
		raxml -D -m PROTGAMMALG -s \$i -n \${i%.*}_raxml.nwk -p 2233
	done
	"""
}


// computes the ME tree given an alignment in phylip form. does this at family level through a for loop.
process run_phylo_ME_supertree_aln_swa {
	label 'process_low'
	tag "${fam}"
	container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://athbaltzis/paralogs:v0.10' :
        'athbaltzis/paralogs:v0.10' }"


	input:
	tuple val(fam), path(all_phylip_aln)

	output:
	path "*"
	tuple val(fam), path("*_fastme_tree.nwk"), emit: me_bestree

	script:
	"""
    # extention .phylip is replaced with 
	for i in \$(echo ${all_phylip_aln}); do
		fastme -i \$i -p -g 1.0 -s -n
	done
	"""
}