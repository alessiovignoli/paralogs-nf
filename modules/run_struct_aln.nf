#!/usr/env nextflow

nextflow.enable.dsl=2

process run_struct_ref {
	label 'cpu'
	tag "${fam}"
	publishDir "${fam}/results/AF2", mode: 'copy'

	input:
	tuple val(fam), path(fasta), path(af2_models), path(code_file)
	val ref_species

	output:
	path ("*")

	shell:
	'''
	for mod in `find *.alphafold.pdb`; do extract_from_pdb -force -infile $mod > test.pdb; mv test.pdb $mod; done
	for i in `grep ">" !{fasta} | tr -d ">"`; do echo -e ">"$i "_P_" "$i"'.alphafold.pdb'; done > !{fam}.!{ref_species}_alphafold.template_list

	t_coffee !{fasta} -template_file !{fam}.!{ref_species}_alphafold.template_list -method TMalign_pair -output fasta_aln -outfile !{fam}_tmalign_alphafold.fa
	mv !{fam}_tmalign_alphafold.fa !{fam}_tmalign_alphafold.fa_original
	t_coffee -other_pg seq_reformat -code file.code_name -in !{fam}_tmalign_alphafold.fa_original > !{fam}_tmalign_alphafold.fa
	t_coffee -other_pg seq_reformat -in !{fam}_tmalign_alphafold.fa -output phylip_aln > !{fam}_tmalign_alphafold.ph  
	fastme -i !{fam}_tmalign_alphafold.ph -m BioNJ -p LG -g 1.0 -s -n -z 5
	raxml -D -m PROTGAMMALG -s !{fam}_tmalign_alphafold.ph -n !{fam}_tmalign_alphafold_raxml.nwk -p 2233
	'''
}
