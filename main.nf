#!/usr/env nextflow

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
/*
include { validateParameters; paramsHelp } from 'plugin/nf-validation'

// Print help message if needed
if (params.help) {
    log.info paramsHelp("nextflow run main.nf --csv my_data.csv --exp_conf exp.yaml --model my_model.py --tune_conf tune.yaml -profie your_profile ")
    System.exit(0)
}

// Validate input parameters
if (params.validate_params) {
    validateParameters()
}
*/
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    input handling
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
params.input_fasta = "*.fasta"
if (params.input_fasta) {
	Channel
	.fromPath(params.input_fasta)
	.map { it -> [it.baseName,it]}
	.set{in_fasta}
}

Channel.fromPath(params.AF2).map { it -> [it.toString().split('\\/')[-4],it]}.groupTuple().set{af2_models}
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { extract_fasta_per_species } from './modules/supermatrix_and_supertree.nf'
include { extract_fasta_per_species_for_full_aln } from './modules/supermatrix_and_supertree.nf'
include { run_alns } from './modules/supermatrix_and_supertree.nf'
include { run_full_alns } from './modules/supermatrix_and_supertree.nf'
include { run_phylo_full } from './modules/supermatrix_and_supertree.nf'
include { concatenate_alns } from './modules/supermatrix_and_supertree.nf'
include { concatenate_alns as concatenate_struct_alns } from './modules/supermatrix_and_supertree.nf'

include { extract_fasta_aln_per_species_sim } from './modules/supermatrix_and_supertree_sim.nf'
include { extract_fasta_aln_per_species_sim_for_full_aln } from './modules/supermatrix_and_supertree_sim.nf'
include { concatenate_alns_sim } from './modules/supermatrix_and_supertree_sim.nf'
include { run_phylo_ML_full_sim as run_phylo_ML_full_aln_sim } from './modules/supermatrix_and_supertree_sim.nf'
include { run_phylo_ML_full_sim as run_phylo_ML_supermatrix_aln_sim } from './modules/supermatrix_and_supertree_sim.nf'
include { run_phylo_ML_full_sim as run_phylo_ML_supertree_aln_sim } from './modules/supermatrix_and_supertree_sim.nf'
include { run_phylo_ME_full_sim as run_phylo_ME_full_aln_sim } from './modules/supermatrix_and_supertree_sim.nf'
include { run_phylo_ME_full_sim as run_phylo_ME_supermatrix_aln_sim } from './modules/supermatrix_and_supertree_sim.nf'
include { run_phylo_ME_full_sim as run_phylo_ME_supertree_aln_sim } from './modules/supermatrix_and_supertree_sim.nf'
include { only_concatenate_aln_sim } from './modules/supermatrix_and_supertree_sim.nf'
include { extract_species_submsa as extract_species_submsa_ML } from './modules/supermatrix_and_supertree_sim.nf'
include { extract_species_submsa as extract_species_submsa_ME } from './modules/supermatrix_and_supertree_sim.nf'
include { superfine } from './modules/supermatrix_and_supertree_sim.nf'
include { superfine as superfine_supertree } from './modules/supermatrix_and_supertree_sim.nf'
include { from_fasta_to_phylip } from './modules/supermatrix_and_supertree_sim.nf'

include { run_colabfold } from './modules/run_struct_model.nf'
include { split_multi_fasta } from './modules/run_struct_model.nf'
include { run_alphafold2 } from './modules/run_struct_model.nf'

include { run_struct_aln } from './modules/run_struct_aln.nf'


// TODO remove all the below

//include { data_preparation } from './modules/data_preparation.nf'
//include { extract_fasta_per_species; extract_fasta_per_species_for_full_aln; run_alns; run_full_alns; run_phylo_full; concatenate_alns; concatenate_alns as concatenate_struct_alns } from './modules/supermatrix_and_supertree.nf'
//include { extract_fasta_aln_per_species_sim; extract_fasta_aln_per_species_sim_for_full_aln; concatenate_alns_sim; run_phylo_ML_full_sim as run_phylo_ML_full_aln_sim; run_phylo_ML_full_sim as run_phylo_ML_supermatrix_aln_sim; run_phylo_ME_full_sim as run_phylo_ME_full_aln_sim; run_phylo_ME_full_sim as run_phylo_ME_supermatrix_aln_sim; only_concatenate_aln_sim} from './modules/supermatrix_and_supertree_sim.nf'
//include {run_colabfold; split_multi_fasta; run_alphafold2} from './modules/run_struct_model.nf'
//include {run_struct_aln} from './modules/run_struct_aln.nf'


//params.data = 'sim'
//params.db = "/users/cn/abaltzis/db/colabfolddb"
//params.ref = "MOUSE"
//params.species_num = 6
//params.mode = 'tcoffee'
//params.orthologs_ids = "PF*/PF*.orthologs_org_ids_to_concatenate"
//params.fasta = "PF*/PF*_domain_sequences_intersect_with_ref_after_OMA.fasta"
//params.intersecting_genes = "PF*/PF*.intersecting_genes"
//params.AF2 = "PF*/results/AF2/*.pdb"
//params.species_num_sim = 25
//params.input_sim = "*.ma"
//params.orthologs_ids_sim = "orthologs_org_ids_to_concatenate"


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow pfam_data_with_AF2 {

        take:
        fasta
        intersecting_genes
        orthologs_ids

        main:

        // first prepare and group together the inputs
        grouped_fasta = Channel.fromPath(fasta)
                        .map { it -> [it.baseName.toString().split("\\.")[0], it]}
	                .groupTuple()
        intersect     = Channel.fromPath(intersecting_genes)
                        .map { it -> [it.baseName.toString().split("\\.")[0], it]}
        ortho         = Channel.fromPath(orthologs_ids)
                        .map { it -> [it.baseName.toString().split("\\.")[0], it]}
        input_fasta   = intersect.combine(ortho, by:0).combine(grouped_fasta, by:0)
        
        // gathering the correct sequence inputs for each PFAM family
        extract_fasta_per_species(input_fasta)
        aln_input = extract_fasta_per_species.out.selected_fasta.transpose()
        extract_fasta_per_species_for_full_aln(input_fasta)

        //1   run_alns(aln_input, params.mode)
        //1   extract_fasta_per_species_for_full_aln.out.combine(run_alns.out.code_name.groupTuple().map{it -> [it[0], it[1][0]]}, by:0).set{run_full_alns_input}
        //1   run_full_alns(run_full_alns_input, params.mode)
        //1   run_phylo_full(run_full_alns.out.full_aln, params.mode)
        //run_phylo_recon(run_full_alns.out.transpose(),params.mode)
        //run_alns.out.aln_output.filter(~/^((?!MOUSE_domain).)*$/).groupTuple().set{alns_grouped}
        //concatenate_alns(alns_grouped,params.mode,params.species_num)
                //run_colabfold.out.colabfold_models.combine(run_alns.out.code_name.groupTuple().map{it -> [it[0], it[1][0]]}, by:0).set{struct_aln_input}
        //run_struct_aln(struct_aln_input)
        //run_struct_aln.out.struct_aln_output.groupTuple().set{struct_alns_grouped}
        //concatenate_struct_alns(struct_alns_grouped,'3DCoffee',params.species_num)
        
}

// TODO nout used so remove once sure
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

        take:
        input_sim
        orthologs_ids

        main:
        // first prepare the inputs
        input_aln_sim     = Channel.fromPath(input_sim)
                            .map{it -> [it.baseName, it]}
        orthologs_ids_sim = Channel.fromPath(orthologs_ids)
        
        // BigTree compotutation
        // takes a MSA and a list of species names, -> phylip format MSA with all sequence from the same species one after the other and header names changed into code values.
        extract_fasta_aln_per_species_sim_for_full_aln(input_aln_sim.combine(orthologs_ids_sim))
        // computes the ML tree from previous phylip alignment
        run_phylo_ML_full_aln_sim(extract_fasta_aln_per_species_sim_for_full_aln.out.phylip_full_aln_sim)
        // computes the ME tree previous phylip alignmen
        run_phylo_ME_full_aln_sim(extract_fasta_aln_per_species_sim_for_full_aln.out.phylip_full_aln_sim)
        // extracting the submsas per species from the ML trees connecting the tree with the codename file that is nedded to ranem the tree tips labels
        extract_inputs_ML = run_phylo_ML_full_aln_sim.out.ml_bestree.join(extract_fasta_aln_per_species_sim_for_full_aln.out.msa_code_name)
        extract_species_submsa_ML(extract_inputs_ML)
        tobe_merged_ML = extract_species_submsa_ML.out.species_subtrees.map{
                it -> [(it[0] + " - ML"), it[1]]
        }
        // extracting the submsas per species from the ME trees connecting the tree with the codename file that is nedded to ranem the tree tips labels
        extract_inputs_ME = run_phylo_ME_full_aln_sim.out.me_bestree.join(extract_fasta_aln_per_species_sim_for_full_aln.out.msa_code_name)
        extract_species_submsa_ME(extract_inputs_ME)
        tobe_merged_ME = extract_species_submsa_ME.out.species_subtrees.map{
                it -> [(it[0] + " - ME"), it[1]]
        }
        // put in the same channel ME and ML extracted trees
        ready_for_superfine = tobe_merged_ML.concat(tobe_merged_ME)
        // use the superfine program to merge all species_submsa into one paralog tree
        superfine(ready_for_superfine)

        
        // from the family MSA extraxct  one msa in fasta format per species, each with all sequences in input fasta from the same species. Order of sequences is preserved.
        extract_fasta_aln_per_species_sim(input_aln_sim.combine(orthologs_ids_sim))

        // SuperMatrix with all species and no shuffle.
        // gets the list of MSAs form each family extracted above (one family at the time) and concatenates all first sequences of all MSA in the list into the same line, same for the second and so on. Order of concatenation is not the order of mSA filenames in the list. The output is a single concatenated MSA in phylip format (columns are not shuffled yet). 
        only_concatenate_aln_sim(extract_fasta_aln_per_species_sim.out.all_aln_sim)
        run_phylo_ML_supermatrix_aln_sim(only_concatenate_aln_sim.out.phylip_only_concatenate_aln_sim)
        run_phylo_ME_supermatrix_aln_sim(only_concatenate_aln_sim.out.phylip_only_concatenate_aln_sim)
        
        // SuperTree with all species and no shuffle.
        // gets the list of MSAs form each family extracted above (one family at the time) and for each MSA it computes the (paralog) tree. Then thanks to superfine it merges all  this (25) trees into a singel (paralog) tree.
        // the creation of the species paralogs happens in parallel for each of them (not a for loop). So first the above list is divided in family identifier + single species MSA.
        msas_for_supertree = extract_fasta_aln_per_species_sim.out.all_aln_sim.transpose()
        // then each msa is made into phylip format and then the species paralog tree is made either using ME or ML.
        from_fasta_to_phylip(msas_for_supertree)
        run_phylo_ML_supertree_aln_sim(from_fasta_to_phylip.out.phylip)
        run_phylo_ME_supertree_aln_sim(from_fasta_to_phylip.out.phylip)
        // now put toghther all ML trees coming from the same family, do the same for ME. Then merge the two channles
        tobe_merged_supertree_ML = run_phylo_ML_supertree_aln_sim.out.ml_bestree.map{
                it -> [(it[0] + " - ML"), it[1]]
        }.groupTuple()
        tobe_merged_supertree_ME = run_phylo_ME_supertree_aln_sim.out.me_bestree.map{
                it -> [(it[0] + " - ME"), it[1]]
        }.groupTuple()
        ready_for_superfine_supertree = tobe_merged_supertree_ML.concat(tobe_merged_supertree_ME)
        // finally merge together using superfine all paralog trees from the same family.
        superfine_supertree(ready_for_superfine_supertree)


        // The following process will do the whole SuperMatrix and SuperTree computation. With 10 replication as well. it will generate all the relevant trees -> 500 (250 ST and 250 SM). 250 = 25 species/units * 10 replicates/samples. 
        concatenate_alns_sim(extract_fasta_aln_per_species_sim.out.all_aln_sim, params.species_num_sim)
        
}

workflow {

        if (params.data == 'pfam') {

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
                pfam_data_with_AF2(params.fasta, params.intersecting_genes, params.orthologs_ids)

        } else if (params.data == 'sim') {

                log.info """\
                Paralogs Analysis - version 0.1
                =====================================
                Input simulated paralogs datasets (FASTA)       : ${params.input_sim}
                Number of species		                : ${params.species_num_sim}
                Input species names		                : ${params.orthologs_ids_sim}
                """
                .stripIndent()
                simulated_data(params.input_sim, params.orthologs_ids_sim)

        } else {

                log.info "the data parameter has to be either 'pfam' or 'sim', pointing to the type of analysis and data to be used. \ngiven : ${params.data}"
                exit 1

        }

}

workflow.onComplete {
        println "Pipeline completed at: $workflow.complete"
        println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
