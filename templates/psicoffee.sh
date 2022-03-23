t_coffee -other_pg seq_reformat -in ${fasta} -output code_name > file.code_name
t_coffee -other_pg seq_reformat -code file.code_name -in ${fasta} > coded.fasta	
t_coffee coded.fasta -mode psicoffee -output fasta -outfile ${output_alns} -thread 2
t_coffee -other_pg seq_reformat -in ${output_alns} -output phylip_aln > ${output_ph}

