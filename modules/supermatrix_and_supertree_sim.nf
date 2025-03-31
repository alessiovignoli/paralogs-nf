#!/usr/env nextflow


// takes a MSA and a list of species names, -> one msa in fasta format per species, each with all sequences in input fasta from the same species. Order of sequences is preserved.
process extract_fasta_aln_per_species_sim {
	label 'process_low'
	tag "${fam}"

	input:
	tuple val(fam), path(init_aln), path(orthologs_ids)

	output:
	tuple val(fam), path("${fam}*.fasta_aln"), emit: all_aln_sim
	//path "*"

	script:
	"""
	extract_seq_per_species_sim.pl ${orthologs_ids} ${init_aln} -after
	"""
}


// takes a MSA and a list of species names, -> phylip format MSA with all sequence from the same species one after the other and header names changed into code values. Order of sequences is preserved.
process extract_fasta_aln_per_species_sim_for_full_aln {
	label 'process_low'
	tag "${fam}"

	input:
	tuple val(fam), path(init_aln), path(orthologs_ids)

	output:
	tuple val(fam), path("${fam}_full_aln_coded.phylip"), emit: phylip_full_aln_sim
	tuple val(fam), path("full_aln.code_name"), path(orthologs_ids), emit: msa_code_name

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
	tuple val(fam), path("RAxML_bestTree*"), emit: ml_bestree

	script:
	raxml_output = phylip_aln.baseName + "_raxml.nwk"
	"""
	raxml -D -m PROTGAMMALG -s ${phylip_aln} -n ${raxml_output} -p 2233
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
	tuple val(fam), path("*_fastme_tree.nwk"), emit: me_bestree

	script:
	"""
	fastme -i ${phylip_aln} -p -g 1.0 -s -n
	"""
}

// given a gene family tree extracts all the subfamilies trees and puts them into a file ready to be used by superfine 
process extract_species_submsa {
	label 'process_low'
	tag "${fam}"

	input:
	tuple val(fam), path(gene_family_tree), path(msa_code_name_file), path(orthologs_ids)

	output:
	tuple val(fam), path("*"), emit: species_subtrees

	script:
	output_name = gene_family_tree.baseName + ".species_subtrees.nwk"
	"""
	extract_species_subtrees.R ${gene_family_tree} ${msa_code_name_file} ${orthologs_ids}  ${output_name} ${params.data}
	"""

}

// mergees a set of trees using the program superfine. It can work with one file containg many trees or many files containing only one tree.
process superfine {
	label 'process_low'
	tag "${id}"

	input:
	tuple val(id), path(trees_tobe_merged)

	output:
	path output

	script:
	fam = id.split(" ")[0]
	output = ("${trees_tobe_merged}".split(" ")[0].split("\\."))[0..-3].join(".") + "_paralog_tree.nwk"
	"""
	cat ${trees_tobe_merged} > tmp
	/SuperFine/runSuperFine.py -r gmrp tmp > ${output}
	"""

}

// transforms a fasta MSA into a phylip format MSA.
process from_fasta_to_phylip {
	label 'process_low'
	tag "${msa_fasta}"

	input:
	tuple val(id), path(msa_fasta)

	output:
	tuple val(id), path(output), emit: phylip

	script:
	output = ("${msa_fasta}".split("\\."))[0..-2].join(".") + ".phylip"
	"""
	t_coffee -other_pg seq_reformat -in ${msa_fasta} -output phylip_aln > ${output}
	"""
}

// gets a list of MSA with same number of sequences in the same order, meaning that the first sequence in MSA 1 is expected to be related (ortholog) to sequence 1 in MSA2 3 4 .. ecc.. 
// it spits out a single MSA in phylip format with number of sequences the MSA1 num seq, but with len = num MSA * len aln in MSA. basucally puts all first sequences of all MSA into the same line (concatenating them), same for the second and so on. Order of concatenation is not the order of mSA filenames in the list. 
process only_concatenate_aln_sim {
	label 'process_low'
	tag "${fam}"

	input:
	tuple val(fam), val(all_aln)

	output:
	tuple val(fam), path("${fam}_supermatrix.phylip"), emit: phylip_only_concatenate_aln_sim

	shell:
	'''
	concatenate_aln.py -f !{fam} -S !{all_aln.join(" ")}
	'''
}


process swap_seq_in_fam {
	label 'process_low'
	tag "${fam}"

	input:
	tuple val(fam), val(num_species_swap), val(num_swap), path(all_aln)

	output:
	//tuple val(fam), path("${fam}_supermatrix.phylip"), emit: phylip_only_concatenate_aln_sim

	script:
	"""
	for i in \$(ls | sort -R | head -n ${num_species_swap}); do
		echo  ${num_swap} \$i 
	done
	"""
}


// The following process will do the whole SuperMatrix and SuperTree computation. With 10 replication as well.
// concatenate_aln.py  will do the following:
// create all the MSA units and single units in phylip format. From the starting species alignments by shuffling and drowing columns.
//      unit = divide the concatenated shuffled MSA in 25 blocks, block1 is unit1, block1 + block2 is unit2 and so on.
//      single_unit = divide the concatenated shuffled MSA in 25 blocks, each block is is own single unit.
// this above process is repeated 10 times, so to generate statistical replicates called sample.
//      so will have unit1_sample1 that will be different from unit1_sample2 and so on.
// For all the above MSA (25 species/units 10 replicates/samples) -> 250 MSA unit, 250 MSa single_unit
// 	    for each unit MSA make raxml (ML) (250) and fastme (ME) (250 trees -> this are the SuperMatrix ML and Me final trees (500)
//      for each single_unit MSA make raxml (ML) (250) and fastme (ME) (250) trees
// the bash for loop will then take the (250) single_unit trees and append them in a file in blocks like: single_unit1 tree  (file1), single_unit1 tree  + single_unit2 tree (file2) ecc.. . this is done for both ME and ML trees respectively.
//      this (500) "appended" file are fed to superfine that will do the SuperTree approach, this are the final SuperTree trees.
process concatenate_alns_sim {
	label 'process_high'
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
