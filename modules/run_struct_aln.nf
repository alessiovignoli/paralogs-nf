#!/usr/env nextflow

process run_struct_aln {
	label 'cpu'
	tag "${fasta}"
	publishDir "${fam}/results/3DCoffee", mode: 'copy'

	input:
	tuple val(fam), path(fasta), path(models), path(code_file)

	output:
	tuple val(fam), path("${fasta.baseName}_tmalign_colabfold_coded.fasta_aln"), emit: struct_aln_output

	shell:
	'''
	for mod in `find *.colabfold.pdb`; do extract_from_pdb -force -infile $mod > test.pdb; mv test.pdb $mod; done
	for i in `grep ">" !{fasta} | tr -d ">"`; do echo -e ">"$i "_P_" "$i"'.colabfold.pdb'; done > !{fasta.baseName}.colabfold.template_list
	t_coffee !{fasta} -template_file !{fasta.baseName}.colabfold.template_list -method TMalign_pair -output fasta_aln -outfile !{fasta.baseName}_tmalign_colabfold.fasta_aln
	t_coffee -other_pg seq_reformat -code file.code_name -in !{fasta.baseName}_tmalign_colabfold.fasta_aln > !{fasta.baseName}_tmalign_colabfold_coded.fasta_aln
	t_coffee -other_pg seq_reformat -in !{fasta.baseName}_tmalign_colabfold_coded.fasta_aln -output phylip_aln > !{fasta.baseName}_tmalign_colabfold_coded.phylip
	'''
}
