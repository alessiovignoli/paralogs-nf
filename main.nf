#!/usr/env nextflow

nextflow.enable.dsl=2

//include { data_preparation } from './modules/data_preparation.nf'
include { extract_fasta_per_species; run_alns; concatenate_alns} from './modules/supermatrix_and_supertree.nf'
include { extract_fasta_aln_per_species_sim; concatenate_alns_sim} from './modules/supermatrix_and_supertree_sim.nf'
include {split_multi_fasta; run_alphafold2} from './modules/run_af2.nf'
include {run_struct_ref} from './modules/run_struct_aln.nf'

params.data = 'pfam'

// data preparation inputs
/*
params.input_fasta = "*.fasta"
if (params.input_fasta) {
	Channel
	.fromPath(params.input_fasta)
	.map { it -> [it.baseName,it]}
	.set{in_fasta}
}*/

params.ref = "MOUSE"
params.species_num = 7
params.mode = 'tcoffee'
params.orthologs_ids = "PF*/PF*.orthologs_org_ids_to_concatenate"
params.fasta = "PF*/PF*_domain_sequences_intersect_with_ref_after_OMA.fasta"
params.intersecting_genes = "PF*/PF*.intersecting_genes"
params.AF2 = "PF*/results/AF2/*.pdb"


if (params.fasta) {
        Channel
	.fromPath(params.fasta)
	.map { it -> [it.baseName.toString().split("\\.")[0],it]}
	.groupTuple()
	.set{grouped_fasta}

        Channel
        .fromPath(params.fasta)
        .filter(~/.*${params.ref}_domain_sequences_intersect.*/)
        .map { it -> [it.baseName.toString().split("\\.")[0],it]}
        .set{fasta_for_af2}
}
if (params.intersecting_genes) {
        Channel
        .fromPath(params.intersecting_genes)
        .map { it -> [it.baseName.toString().split("\\.")[0],it]}
        .set{intersect}
}
if (params.orthologs_ids) {
        Channel
        .fromPath(params.orthologs_ids)
        .map { it -> [it.baseName.toString().split("\\.")[0],it]}
        .set{ortho}
}
intersect.combine(ortho,by:0).combine(grouped_fasta,by:0).set{input_fasta}
Channel.fromPath(params.AF2).map { it -> [it.toString().split('\\/')[-4],it]}.groupTuple().set{af2_models}



params.species_num_sim = 25
params.input_sim = "*.ma"
params.orthologs_ids_sim = "orthologs_org_ids_to_concatenate"
Channel.fromPath(params.input_sim).map{it -> [it.baseName,it]}.set {input_aln_sim}
Channel.fromPath(params.orthologs_ids_sim).set{orthologs_ids_sim}

if(params.data == 'pfam') {
        log.info """\
Paralogs Analysis - version 0.1
=====================================
Input paralogs datasets (FASTA)		        : ${params.fasta}
Number of species		                : ${params.species_num}
Input paralog genes		                : ${params.intersecting_genes}
Input species names		                : ${params.orthologs_ids}
MSA algorithm (TCoffee, PSICoffee, MAFFT-Ginsi) : ${params.mode}
"""
.stripIndent()
}
else {
        log.info """\
Paralogs Analysis - version 0.1
=====================================
Input simulated paralogs datasets (FASTA)       : ${params.input_sim}
Number of species		                : ${params.species_num_sim}
Input species names		                : ${params.orthologs_ids_sim}
"""
.stripIndent()
}

workflow pfam_data_with_AF2 {
        extract_fasta_per_species(input_fasta,params.mode)
        extract_fasta_per_species.out.transpose().set{aln_input}
        run_alns(aln_input,params.mode)
        run_alns.out.aln_output.groupTuple().set{alns_grouped}
        concatenate_alns(alns_grouped,params.mode,params.species_num)
        split_multi_fasta(fasta_for_af2)
        run_alphafold2(split_multi_fasta.out.transpose())
        extract_fasta_per_species.out.transpose().filter(~/.*${params.ref}_domain_sequences_after_intersection.*/).combine(run_alphafold2.out.af2_models.groupTuple(), by: 0).combine(run_alns.out.code_name.groupTuple().map{it -> [it[0], it[1][0]]}, by:0).set{struct_ref_input}
        run_struct_ref(struct_ref_input,params.ref)
}

workflow pfam_data_without_AF2 {
        extract_fasta_per_species(input_fasta,params.mode)
        extract_fasta_per_species.out.transpose().set{aln_input}
        run_alns(aln_input,params.mode)
        run_alns.out.aln_output.groupTuple().set{alns_grouped}
        concatenate_alns(alns_grouped,params.mode,params.species_num)
        extract_fasta_per_species.out.transpose().filter(~/.*${params.ref}_domain_sequences_after_intersection.*/).combine(af2_models, by: 0).combine(run_alns.out.code_name.groupTuple().map{it -> [it[0], it[1][0]]}, by:0).set{struct_ref_input}
        run_struct_ref(struct_ref_input,params.ref)
}

workflow simulated_data {
        extract_fasta_aln_per_species_sim(input_aln_sim.combine(orthologs_ids_sim))
        concatenate_alns_sim(extract_fasta_aln_per_species_sim.out.all_aln_sim,params.species_num_sim)
}

workflow {
        if (params.data == 'pfam')
                pfam_data_without_AF2()
        else
                simulated_data()
}

workflow.onComplete {
    println "Pipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
