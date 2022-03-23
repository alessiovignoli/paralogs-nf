#!/usr/bin/env bash

for fam in `find PF* -maxdepth 0 -type d`; do
	cd $fam/results/tcoffee
	echo $fam
	
	orgsnum=7
	for num in $(seq 1 $orgsnum); do
		grep "Explained" unit_"$num"_*stat.txt | grep -o '\-*[0-9]\.[0-9]*' > unit_"$num"_explained_variance.dat
		cat `find unit_"$num"_sample_*_concatenated_aln.phylip_fastme_tree.nwk -maxdepth 0 | sort -V | tr "\n" "\t"` > unit_"$num"_all.nwk
		cat `find RAxML_bestTree.unit_"$num"_sample_*_concatenated_aln_raxml.nwk -maxdepth 0 | sort -V | tr "\n" "\t"` > unit_"$num"_all_raxml.nwk
		cat `find unit_"$num"_sample_*_concatenated_aln.phylip_fastme_superfine_tree.nwk -maxdepth 0 | sort -V | tr "\n" "\t"` > unit_"$num"_all_superfine_trees.nwk
		cat `find unit_"$num"_sample_*_concatenated_aln.raxml_superfine_tree.nwk -maxdepth 0 | sort -V | tr "\n" "\t"` > unit_"$num"_all_RAxML_superfine_trees.nwk
	done

	for sample in {1..10}; do 
		cat `find unit*sample_"$sample"_concatenated_aln.phylip_fastme_tree.nwk -maxdepth 0 | sort -V | tr "\n" "\t"` > sample_"$sample"_all.nwk
		cat `find RAxML_bestTree.unit_*_sample_"$sample"_concatenated_aln_raxml.nwk -maxdepth 0 | sort -V | tr "\n" "\t"` > sample_"$sample"_all_raxml.nwk
		cat `find unit*sample_"$sample"_concatenated_aln.phylip_fastme_superfine_tree.nwk -maxdepth 0 | sort -V | tr "\n" "\t"` > sample_"$sample"_all_superfine_trees.nwk
		cat `find unit*sample_"$sample"_concatenated_aln.raxml_superfine_tree.nwk -maxdepth 0 | sort -V | tr "\n" "\t"` > sample_"$sample"_all_RAxML_superfine_trees.nwk
	done

	cat `find unit_*_sample_*_concatenated_aln.phylip_fastme_tree.nwk -maxdepth 0 | sort -V | tr "\n" "\t"` > all_units.nwk
	#cat `find unit_*_sample_*_concatenated_aln.phylip_fastme_MRL_tree.nwk -maxdepth 0 | sort -V | tr "\n" "\t"` > all_units_MRL_trees.nwk
	cat `find unit_*_sample_*_concatenated_aln.phylip_fastme_superfine_tree.nwk -maxdepth 0 | sort -V | tr "\n" "\t"` > all_units_superfine_trees.nwk
	cat `find RAxML_bestTree.unit_*_sample*_concatenated_aln_raxml.nwk -maxdepth 0 | sort -V | tr "\n" "\t"` > all_units_RAxML_trees.nwk
	#cat `find unit_*_sample*_concatenated_aln.raxml_MRL_tree.nwk -maxdepth 0 | sort -V | tr "\n" "\t"` > all_units_RAxML_MRL_trees.nwk
	cat `find unit_*_sample*_concatenated_aln.raxml_superfine_tree.nwk -maxdepth 0 | sort -V | tr "\n" "\t"` > all_units_RAxML_superfine_trees.nwk
	cd ../../../
done

