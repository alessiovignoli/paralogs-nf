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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


include { extract_fasta_per_species } from './modules/supermatrix_and_supertree.nf'
include { extract_fasta_per_species_for_full_aln } from './modules/supermatrix_and_supertree.nf'
include { extract_fasta_aln_per_species_emp } from './modules/supermatrix_and_supertree.nf'
include { run_alns } from './modules/supermatrix_and_supertree.nf'
include { run_full_alns } from './modules/supermatrix_and_supertree.nf'
include { concatenate_alns } from './modules/supermatrix_and_supertree.nf'
include { concatenate_alns as concatenate_struct_alns } from './modules/supermatrix_and_supertree.nf'
include { remove_ref_species } from './modules/supermatrix_and_supertree.nf'

include { extract_fasta_aln_per_species_sim } from './modules/supermatrix_and_supertree_sim.nf'
include { extract_fasta_aln_per_species_sim_for_full_aln } from './modules/supermatrix_and_supertree_sim.nf'
include { concatenate_alns_sim } from './modules/supermatrix_and_supertree_sim.nf'
include { run_phylo_ML_full_sim as run_phylo_ML_full_aln_emp } from './modules/supermatrix_and_supertree_sim.nf'
include { run_phylo_ML_full_sim as run_phylo_ML_supermatrix_aln_emp } from './modules/supermatrix_and_supertree_sim.nf'
include { run_phylo_ML_full_sim as run_phylo_ML_supertree_aln_emp } from './modules/supermatrix_and_supertree_sim.nf'
include { run_phylo_ML_full_sim as run_phylo_ML_full_aln_sim } from './modules/supermatrix_and_supertree_sim.nf'
include { run_phylo_ML_full_sim as run_phylo_ML_supermatrix_aln_sim } from './modules/supermatrix_and_supertree_sim.nf'
include { run_phylo_ML_full_sim as run_phylo_ML_supertree_aln_sim } from './modules/supermatrix_and_supertree_sim.nf'
include { run_phylo_ME_full_sim as run_phylo_ME_full_aln_emp } from './modules/supermatrix_and_supertree_sim.nf'
include { run_phylo_ME_full_sim as run_phylo_ME_supermatrix_aln_emp } from './modules/supermatrix_and_supertree_sim.nf'
include { run_phylo_ME_full_sim as run_phylo_ME_supertree_aln_emp } from './modules/supermatrix_and_supertree_sim.nf'
include { run_phylo_ME_full_sim as run_phylo_ME_full_aln_sim } from './modules/supermatrix_and_supertree_sim.nf'
include { run_phylo_ME_full_sim as run_phylo_ME_supermatrix_aln_sim } from './modules/supermatrix_and_supertree_sim.nf'
include { run_phylo_ME_full_sim as run_phylo_ME_supertree_aln_sim } from './modules/supermatrix_and_supertree_sim.nf'
include { only_concatenate_aln_sim } from './modules/supermatrix_and_supertree_sim.nf'
include { only_concatenate_aln_sim as only_concatenate_aln_emp } from './modules/supermatrix_and_supertree_sim.nf'
include { extract_species_submsa as extract_species_submsa_ML_emp } from './modules/supermatrix_and_supertree_sim.nf'
include { extract_species_submsa as extract_species_submsa_ML_sim } from './modules/supermatrix_and_supertree_sim.nf'
include { extract_species_submsa as extract_species_submsa_ME_emp } from './modules/supermatrix_and_supertree_sim.nf'
include { extract_species_submsa as extract_species_submsa_ME_sim } from './modules/supermatrix_and_supertree_sim.nf'
include { superfine as superfine_emp } from './modules/supermatrix_and_supertree_sim.nf'
include { superfine as superfine_sim } from './modules/supermatrix_and_supertree_sim.nf'
include { superfine as superfine_supertree_emp } from './modules/supermatrix_and_supertree_sim.nf'
include { superfine as superfine_supertree_sim } from './modules/supermatrix_and_supertree_sim.nf'
include { from_fasta_to_phylip as from_fasta_to_phylip_emp } from './modules/supermatrix_and_supertree_sim.nf'
include { from_fasta_to_phylip as from_fasta_to_phylip_sim } from './modules/supermatrix_and_supertree_sim.nf'

include { run_colabfold } from './modules/run_struct_model.nf'
include { split_multi_fasta } from './modules/run_struct_model.nf'
include { run_alphafold2 } from './modules/run_struct_model.nf'

include { run_struct_aln } from './modules/run_struct_aln.nf'


// this is separated from the rest because it might be temporary work.
// this are the alias for the BigTree without ref computation in the empirical step.
include { extract_fasta_per_species_for_full_aln as extract_fasta_per_species_for_full_aln_no_ref } from './modules/supermatrix_and_supertree.nf'
include { run_full_alns as run_full_alns_no_ref } from './modules/supermatrix_and_supertree.nf'
include { run_phylo_ML_full_sim as run_phylo_ML_full_aln_emp_no_ref } from './modules/supermatrix_and_supertree_sim.nf'
include { run_phylo_ME_full_sim as run_phylo_ME_full_aln_emp_no_ref } from './modules/supermatrix_and_supertree_sim.nf'
include { extract_species_submsa as extract_species_submsa_ML_emp_no_ref } from './modules/supermatrix_and_supertree_sim.nf'
include { extract_species_submsa as extract_species_submsa_ME_emp_no_ref } from './modules/supermatrix_and_supertree_sim.nf'
include { superfine as superfine_emp_no_ref } from './modules/supermatrix_and_supertree_sim.nf'

// this is separated from the rest because it might be temporary work.
// analysis of mislabeled paralogs
include { swap_seq_in_fam } from './modules/supermatrix_and_supertree_swa.nf'
include { from_fasta_to_phylip_swa } from './modules/supermatrix_and_supertree_swa.nf'
include { run_phylo_ML_supertree_aln_swa } from './modules/supermatrix_and_supertree_swa.nf'
include { run_phylo_ME_supertree_aln_swa } from './modules/supermatrix_and_supertree_swa.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow empirical_data {

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
        
        // BigTree compotutation
        // gathering all the correct sequence inputs for each PFAM family and putting them into a single file.
        extract_fasta_per_species_for_full_aln(input_fasta)
        // creates the bigtree phylip MSA, using all paralogs ortholog informations. It creates a header id code based om code_file.
        run_full_alns(extract_fasta_per_species_for_full_aln.out.selected_fasta)
        // computes the ML tree from previous phylip alignment
        run_phylo_ML_full_aln_emp(run_full_alns.out.full_aln)
        // computes the ME tree from previous phylip alignmen
        run_phylo_ME_full_aln_emp(run_full_alns.out.full_aln)
        // extracting the submsas per species from the ML trees connecting the tree with the codename file and species names that are nedded to raname the tree tips labels
        extract_inputs_ML = run_phylo_ML_full_aln_emp.out.ml_bestree.join(run_full_alns.out.msa_code_name).join(ortho)
        extract_species_submsa_ML_emp(extract_inputs_ML)
        tobe_merged_ML = extract_species_submsa_ML_emp.out.species_subtrees.map{
                it -> [(it[0] + " - ML"), it[1]]
        }
        // extracting the submsas per species from the ME trees connecting the tree with the codename file that is nedded to ranem the tree tips labels
        extract_inputs_ME = run_phylo_ME_full_aln_emp.out.me_bestree.join(run_full_alns.out.msa_code_name).join(ortho)
        extract_species_submsa_ME_emp(extract_inputs_ME)
        tobe_merged_ME = extract_species_submsa_ME_emp.out.species_subtrees.map{
                it -> [(it[0] + " - ME"), it[1]]
        }
        // put in the same channel ME and ML extracted trees
        ready_for_superfine = tobe_merged_ML.concat(tobe_merged_ME)
        // use the superfine program to merge all species_submsa into one paralog tree
        superfine_emp(ready_for_superfine)
        
        // BigTree compotutation but without the reference species sequences (aka params.ref).
        // process are aliased and repeated for simplicity for the computing running time in analysis step.
        // remove the keyword ref fromthe ortho file. Outputs directly what is needed by folowing processes.
        remove_ref_species(input_fasta)
        // repeat the same steps as in the BigTree approach as above from this point on. Basically do the same but without ref. 
        extract_fasta_per_species_for_full_aln_no_ref(remove_ref_species.out.no_ref)
        run_full_alns_no_ref(extract_fasta_per_species_for_full_aln_no_ref.out.selected_fasta)
        run_phylo_ML_full_aln_emp_no_ref(run_full_alns_no_ref.out.full_aln)
        run_phylo_ME_full_aln_emp_no_ref(run_full_alns_no_ref.out.full_aln)
        extract_inputs_ML_no_ref = run_phylo_ML_full_aln_emp_no_ref.out.ml_bestree.join(run_full_alns_no_ref.out.msa_code_name).join(remove_ref_species.out.modified_ortho)
        extract_species_submsa_ML_emp_no_ref(extract_inputs_ML_no_ref)
        tobe_merged_ML_no_ref = extract_species_submsa_ML_emp_no_ref.out.species_subtrees.map{
                it -> [(it[0] + " - ML"), it[1]]
        }
        extract_inputs_ME_no_ref = run_phylo_ME_full_aln_emp_no_ref.out.me_bestree.join(run_full_alns_no_ref.out.msa_code_name).join(remove_ref_species.out.modified_ortho)
        extract_species_submsa_ME_emp_no_ref(extract_inputs_ME_no_ref)
        tobe_merged_ME_no_ref = extract_species_submsa_ME_emp_no_ref.out.species_subtrees.map{
                it -> [(it[0] + " - ME"), it[1]]
        }
        ready_for_superfine_no_ref = tobe_merged_ML_no_ref.concat(tobe_merged_ME_no_ref)
        superfine_emp_no_ref(ready_for_superfine_no_ref)

        // Step in common to ST and SM for compuattion of running times. 
        // gets the MSA computed on all species but ref. Then from this MSA it extracts the species subMSAs.
        // this allows to compare the number of columnns and seq in the MSA matrix between SM, ST and BigTree no ref since the starting MSA is the same. 
        input_extract_submsa = run_full_alns_no_ref.out.fasta_full_aln.join(remove_ref_species.out.modified_ortho)
        extract_fasta_aln_per_species_emp(input_extract_submsa)

        // SuperMatrix with all species but not reference and no shuffle. Done for computation of running time.
        // gets the list of MSAs form each family extracted above (one family at the time) and concatenates all first sequences of all MSA in the list into the same line, same for the second and so on. Order of concatenation is not the order of mSA filenames in the list. The output is a single concatenated MSA in phylip format (columns are not shuffled yet). 
	only_concatenate_aln_emp(extract_fasta_aln_per_species_emp.out.species_aln)
	run_phylo_ML_supermatrix_aln_emp(only_concatenate_aln_emp.out.phylip_only_concatenate_aln_sim)
        run_phylo_ME_supermatrix_aln_emp(only_concatenate_aln_emp.out.phylip_only_concatenate_aln_sim)

        // SuperTree with all species  but not reference  and no shuffle. Done for computation of running time.
        // gets the list of MSAs form each family extracted above (one family at the time) and for each MSA it computes the (paralog) tree. Then thanks to superfine it merges all  this (6) trees into a singel (paralog) tree.
        // the creation of the species paralogs happens in parallel for each of them (not a for loop). So first the above list is divided in family identifier + single species MSA.
        msas_for_supertree = extract_fasta_aln_per_species_emp.out.species_aln.transpose()
        // then each msa is made into phylip format and then the species paralog tree is made either using ME or ML.
        from_fasta_to_phylip_emp(msas_for_supertree)
        run_phylo_ML_supertree_aln_emp(from_fasta_to_phylip_emp.out.phylip)
        run_phylo_ME_supertree_aln_emp(from_fasta_to_phylip_emp.out.phylip)
        // now put toghther all ML trees coming from the same family, do the same for ME. Then merge the two channles
        tobe_merged_supertree_ML = run_phylo_ML_supertree_aln_emp.out.ml_bestree.map{
                it -> [(it[0] + " - ML"), it[1]]
        }.groupTuple()
        tobe_merged_supertree_ME = run_phylo_ME_supertree_aln_emp.out.me_bestree.map{
                it -> [(it[0] + " - ME"), it[1]]
        }.groupTuple()
        ready_for_superfine_supertree = tobe_merged_supertree_ML.concat(tobe_merged_supertree_ME)
        // finally merge together using superfine all paralog trees from the same family.
        superfine_supertree_emp(ready_for_superfine_supertree)

        // Step in common to supermatrix and supertree. species alignment preparation.
        // gathering all the correct sequence inputs for each PFAM family.
        extract_fasta_per_species(input_fasta)
        aln_input = extract_fasta_per_species.out.selected_fasta.transpose()
        // for each species in each family it alligns  the paralogs by themselves, (no orthologs relationships).
        run_alns(aln_input, params.mode)
        // remove mouse paralogs fron the supermatrix supertree approach. this is done to check for circularity in the analysis step.
	alns_grouped = run_alns.out.aln_output.filter{ it -> !it.toString().contains("${params.ref}_domain") }.groupTuple()

        // The following process will do the whole SuperMatrix and SuperTree computation. With 10 replication as well. it will generate all the relevant trees -> 500 (60 ST and 60 SM). 60 = 6 species/units * 10 replicates/samples. 
        concatenate_alns(alns_grouped, params.species_num)

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
        // computes the ME tree from previous phylip alignmen
        run_phylo_ME_full_aln_sim(extract_fasta_aln_per_species_sim_for_full_aln.out.phylip_full_aln_sim)
        // extracting the submsas per species from the ML trees connecting the tree with the codename file that is nedded to raname the tree tips labels
        extract_inputs_ML = run_phylo_ML_full_aln_sim.out.ml_bestree.join(extract_fasta_aln_per_species_sim_for_full_aln.out.msa_code_name)
        extract_species_submsa_ML_sim(extract_inputs_ML)
        tobe_merged_ML = extract_species_submsa_ML_sim.out.species_subtrees.map{
                it -> [(it[0] + " - ML"), it[1]]
        }
        // extracting the submsas per species from the ME trees connecting the tree with the codename file that is nedded to ranem the tree tips labels
        extract_inputs_ME = run_phylo_ME_full_aln_sim.out.me_bestree.join(extract_fasta_aln_per_species_sim_for_full_aln.out.msa_code_name)
        extract_species_submsa_ME_sim(extract_inputs_ME)
        tobe_merged_ME = extract_species_submsa_ME_sim.out.species_subtrees.map{
                it -> [(it[0] + " - ME"), it[1]]
        }
        // put in the same channel ME and ML extracted trees
        ready_for_superfine = tobe_merged_ML.concat(tobe_merged_ME)
        // use the superfine program to merge all species_submsa into one paralog tree
        superfine_sim(ready_for_superfine)

        
        // from the family MSA extraxct one msa in fasta format per species, each with all sequences in input fasta from the same species. Order of sequences is preserved.
        extract_fasta_aln_per_species_sim(input_aln_sim.combine(orthologs_ids_sim))

        // SuperMatrix with all species and no shuffle. Done for computation of running time.
        // gets the list of MSAs form each family extracted above (one family at the time) and concatenates all first sequences of all MSA in the list into the same line, same for the second and so on. Order of concatenation is not the order of mSA filenames in the list. The output is a single concatenated MSA in phylip format (columns are not shuffled yet). 
        only_concatenate_aln_sim(extract_fasta_aln_per_species_sim.out.all_aln_sim)
        run_phylo_ML_supermatrix_aln_sim(only_concatenate_aln_sim.out.phylip_only_concatenate_aln_sim)
        run_phylo_ME_supermatrix_aln_sim(only_concatenate_aln_sim.out.phylip_only_concatenate_aln_sim)
        
        // SuperTree with all species and no shuffle. Done for computation of running time.
        // gets the list of MSAs form each family extracted above (one family at the time) and for each MSA it computes the (paralog) tree. Then thanks to superfine it merges all  this (25) trees into a singel (paralog) tree.
        // the creation of the species paralogs happens in parallel for each of them (not a for loop). So first the above list is divided in family identifier + single species MSA.
        msas_for_supertree = extract_fasta_aln_per_species_sim.out.all_aln_sim.transpose()
        // then each msa is made into phylip format and then the species paralog tree is made either using ME or ML.
        from_fasta_to_phylip_sim(msas_for_supertree)
        run_phylo_ML_supertree_aln_sim(from_fasta_to_phylip_sim.out.phylip)
        run_phylo_ME_supertree_aln_sim(from_fasta_to_phylip_sim.out.phylip)
        // now put toghther all ML trees coming from the same family, do the same for ME. Then merge the two channles
        tobe_merged_supertree_ML = run_phylo_ML_supertree_aln_sim.out.ml_bestree.map{
                it -> [(it[0] + " - ML"), it[1]]
        }.groupTuple()
        tobe_merged_supertree_ME = run_phylo_ME_supertree_aln_sim.out.me_bestree.map{
                it -> [(it[0] + " - ME"), it[1]]
        }.groupTuple()
        ready_for_superfine_supertree = tobe_merged_supertree_ML.concat(tobe_merged_supertree_ME)
        // finally merge together using superfine all paralog trees from the same family.
        superfine_supertree_sim(ready_for_superfine_supertree)


        // The following process will do the whole SuperMatrix and SuperTree computation. With 10 replication as well. it will generate all the relevant trees -> 500 (250 ST and 250 SM). 250 = 25 species/units * 10 replicates/samples. 
        concatenate_alns_sim(extract_fasta_aln_per_species_sim.out.all_aln_sim, params.species_num_sim)
        
}

workflow swapped_data {

        take:
        input_sim
        orthologs_ids

        main:
        // first prepare the inputs
        input_aln_sim     = Channel.fromPath(input_sim)
                            .map{it -> [it.baseName, it]}
        orthologs_ids_sim = Channel.fromPath(orthologs_ids)
        
        // from the family MSA extraxct one msa in fasta format per species, each with all sequences in input fasta from the same species. Order of sequences is preserved.
        extract_fasta_aln_per_species_sim(input_aln_sim.combine(orthologs_ids_sim))

        // 2 sequences could be mislabeld (swapped) in 1 species, or in 2 species and so on
        // the following lines will take care of creating all possibilities for swap up to 10
        // this pair combination will take care of also being used as a key for a given swap "experiment"
        set_species_to_swap = channel.of(1, 2, 3, 4, 5, 6 ,7, 8, 9, 10)
        set_sequence_swap   = channel.of(2, 3, 4, 5, 6 ,7, 8, 9, 10)
        combination_keys    = set_species_to_swap.combine(set_sequence_swap)

        // create the input for the swapper function, for each family there will be 90 
        // possible swappping combinations. Each of wich has to randomly choose which species MSA to add swaps to
        to_be_swap = combination_keys.combine(
                extract_fasta_aln_per_species_sim.out.all_aln_sim
                ).map { it -> [it[2], it[0], it[1], it[3]]}

        // actually swap how many species and how many sequences as required from combination key
        swap_seq_in_fam(to_be_swap)

        // only_concatenate_aln_sim has only one keyword, so the parameters to the swap have to be added to the family name
        to_be_concatenated = swap_seq_in_fam.out.swapped_fastas.map{
                it -> [ it[0] + "_swap_" + it[1] + "_" + it[2], it[3]]
        }

        // SuperMatrix with all species and no shuffle. This section is done only on unit 25.
        // gets the list of MSAs form each family extracted above (one family at the time) and concatenates all first sequences of all MSA in the list into the same line, same for the second and so on. Order of concatenation is not the order of mSA filenames in the list. The output is a single concatenated MSA in phylip format (columns are not shuffled yet). 
        only_concatenate_aln_sim(to_be_concatenated)
        run_phylo_ML_supermatrix_aln_sim(only_concatenate_aln_sim.out.phylip_only_concatenate_aln_sim)
        run_phylo_ME_supermatrix_aln_sim(only_concatenate_aln_sim.out.phylip_only_concatenate_aln_sim)

        //  SuperTree with all species and no shuffle. Done for computation of running time.
        // gets the list of MSAs form each family extracted above (one family at the time) and for each MSA it computes the (paralog) tree. Then thanks to superfine it merges all  this (25) trees into a singel (paralog) tree.
        // the creation of the species paralogs happens in a for loop (all species in same family together).  
        // each msa is made into phylip format and then the species paralog tree is made either using ME or ML.
        from_fasta_to_phylip_swa(to_be_concatenated) 
        run_phylo_ML_supertree_aln_swa(from_fasta_to_phylip_swa.out.phylip)
        run_phylo_ME_supertree_aln_swa(from_fasta_to_phylip_swa.out.phylip)
        // now add a tag to the keyword to ML trees and ME. Then merge the two channles.
        tobe_merged_supertree_ML = run_phylo_ML_supertree_aln_swa.out.ml_bestree.map{
                it -> [(it[0] + " - ML"), it[1]]
        }
        tobe_merged_supertree_ME = run_phylo_ME_supertree_aln_swa.out.me_bestree.map{
                it -> [(it[0] + " - ME"), it[1]]
        }
        ready_for_superfine_supertree = tobe_merged_supertree_ML.concat(tobe_merged_supertree_ME)
        // finally merge together using superfine all paralog trees from the same family.
        superfine_supertree_sim(ready_for_superfine_supertree)

}

workflow {

        if (params.data == 'emp') {

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
                empirical_data(params.fasta, params.intersecting_genes, params.orthologs_ids)

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
        
        } else if (params.data == 'swa') {

                log.info """\
                Paralogs Analysis - version 0.1
                =====================================
                Input simulated paralogs datasets (FASTA)       : ${params.input_sim}
                Number of species		                : ${params.species_num_sim}
                Input species names		                : ${params.orthologs_ids_sim}
                number of paralogs                              : ${params.paralog_num_sim}
                number of max species to swap sequence          : 10
                number of max sequences swapped per species     : 10
                """
                .stripIndent()
                swapped_data(params.input_sim, params.orthologs_ids_sim)

        } else {

                log.info "the data parameter has to be either 'emp', 'sim' or 'swa', pointing to the type of analysis and data to be used. \ngiven : ${params.data}"
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
