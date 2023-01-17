#!/usr/bin/env python2


from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import argparse
import random
import subprocess
from subprocess import check_output
import os
import time
import numpy as np
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

parser = argparse.ArgumentParser()

parser.add_argument("-c","--concat", action='store_true', dest='conc',help='Carry out the main concatenation-based analysis')
parser.add_argument("-m","--method", action='store',choices=['ME', 'ML','B'], dest='method',help='Select tree building method: either ME (FastME) or ML (RaxML) or Both')
parser.add_argument("-s","--shuffle", action='store_true', dest='shuf',help='Randomize the order of input MSAs')
parser.add_argument("-p","--percentage", action='store_true', dest='perc',help='Extract percentages of the columns of each initial alignment and concatenate them')
#parser.add_argument("-b","--bootstrap", action='store', dest='boot', type=int,help='Define the number of BS replicates')
parser.add_argument("-f","--family",action='store', dest='fam',help='Define the family reference name')
parser.add_argument("-o","--output", action='store', dest='output', help='Define an output filename')
parser.add_argument("aln", nargs='+', help='Input MSAs in fasta format')

args = parser.parse_args()
#
# If --shuffle randomize the order of input MSAs
#

if args.shuf==True:
    random.shuffle(args.aln)

#
# Read, concatenate the input MSAs and calculate the number of BS replicates
#
each_len = []
cum_len = []
nor_len = []
cum_avg_len = []
nor_avg_len = []
tot_len=0
avg_tot_len=0
alns_list = []
times = 0
nor_len.append(times)
nor_avg_len.append(times)
for f in args.aln:
	align = AlignIO.read(f,"fasta")
	alns_list.append(align)
	if times == 0:
		supermatrix = align
	else:
		supermatrix = supermatrix + align
	times += 1
	cur_len = align.get_alignment_length()
	tot_len+=cur_len
	cum_len.append(tot_len)
	nor_len.append(tot_len)
	each_len.append(cur_len)
boot = int(tot_len/times)
file = open("BS.dat","w")
file.write("%s\n" % boot)
file.close()

file = open("SD.dat","w")
file.write("%s\n" % np.std(each_len))
file.close()

for f in args.aln:
	align = AlignIO.read(f,"fasta")
	ph_outfile = os.path.splitext(os.path.basename(f))[0] + '.phylip'
	ph_outmat = os.path.splitext(os.path.basename(f))[0] + '.mat'
	ph_reps = os.path.splitext(os.path.basename(f))[0] + '.replicates'
	ph_out = open(ph_outfile,"w")
	AlignIO.write(align,ph_out,"phylip")
	ph_out.close()
	if args.method == 'ME' or args.method == 'B':
		command = "fastme -i " + ph_outfile + " -p -b " + str(boot) + " -g 1.0 -s -n -B " + ph_reps +  " -O " + ph_outmat
		subprocess.call(command, shell=True)
	if args.method == 'ML' or args.method == 'B':
		command = "raxml -D -m PROTGAMMALG -s " + ph_outfile +  " -n " + os.path.splitext(os.path.basename(f))[0] + "_raxml.nwk -p 2233"
		subprocess.call(command, shell=True)
	if args.method == 'ME' or args.method == 'B':
		command = "cat " + ph_outfile + "_fastme_tree.nwk >> " + args.fam + ".all_supertrees.nwk"
		subprocess.call(command, shell=True) 
	if args.method == 'ML' or args.method == 'B':
		command = "cat RAxML_bestTree." + os.path.splitext(os.path.basename(f))[0] + "_raxml.nwk >> " + args.fam + ".all_raxml_supertrees.nwk"
		subprocess.call(command, shell=True)
for i in range(0,times):
	if i == times-1:
		avg_tot_len = tot_len
		cum_avg_len.append(avg_tot_len)
		nor_avg_len.append(avg_tot_len)
	else:
		avg_tot_len += boot
		cum_avg_len.append(avg_tot_len)
		nor_avg_len.append(avg_tot_len)



#
# Extract the sub-MSAs corresponding to each unit and create BS replicate MSAs containing the same number of columns as the first unit
#
if args.conc==True:
	all_len = supermatrix.get_alignment_length()
	align_array = np.array([list(rec) for rec in supermatrix], np.character)
	align_array = align_array.transpose()

	for k in range(1, 11):
		col_index = range(all_len) 
		random.shuffle(col_index)
		np_alns_list = []
		for i in range(1,len(cum_avg_len)+1):
			for j in range(0,cum_avg_len[i-1]):
				if j == 0:
					conc_aln = align_array[col_index[j],:] 
				else:
					conc_aln = np.column_stack((conc_aln,align_array[col_index[j],:]))
			for z in range(nor_avg_len[i-1],nor_avg_len[i]):
				if z == nor_avg_len[i-1]:
					conc_aln_unit = align_array[col_index[z],:]
				else:
					conc_aln_unit = np.column_stack((conc_aln_unit,align_array[col_index[z],:]))
			np_alns_list.append(conc_aln_unit)
			tmp_align = MultipleSeqAlignment([])
			tmp_align_unit = MultipleSeqAlignment([])
			for seq in range(0,len(supermatrix)):
				myseq = conc_aln[seq,:].tolist()
				myseqrecord = SeqRecord(Seq("".join(myseq),IUPAC.protein),id=supermatrix[seq].id)
				tmp_align.append(myseqrecord)
				#tmp_align.add_sequence(supermatrix[seq].id,"".join(myseq))

				myseq_unit = conc_aln_unit[seq,:].tolist()
				myseqrecord_unit = SeqRecord(Seq("".join(myseq_unit),IUPAC.protein),id=supermatrix[seq].id)
				tmp_align_unit.append(myseqrecord_unit)

			outfile = "unit_" + str(i) + "_sample_" + str(k) + "_concatenated_aln.phylip"
			outmat = "unit_" + str(i) + "_sample_" + str(k) + "_concatenated_aln.mat"
			out = open(outfile,"w")
			AlignIO.write(tmp_align,out,"phylip")
			out.close()
			if args.method == 'ME' or args.method == 'B':
				command = "fastme -i " + outfile + " -p -g 1.0 -s -n -O " + outmat
				subprocess.call(command, shell=True)
			if args.method == 'ML' or args.method == 'B':
				command = "raxml -D -m PROTGAMMALG -s " + outfile + " -n " + os.path.splitext(os.path.basename(outfile))[0] + "_raxml.nwk -p 2233"
				subprocess.call(command, shell=True)
			outfile = "single_unit_" + str(i) + "_sample_" + str(k) + "_concatenated_aln.phylip"
			outmat = "single_unit_" + str(i) + "_sample_" + str(k) + "_concatenated_aln.mat"
			out = open(outfile,"w")
			AlignIO.write(tmp_align_unit,out,"phylip")
			out.close()
			if args.method == 'ME' or args.method == 'B':
				command = "fastme -i " + outfile + " -p -g 1.0 -s -n -O " + outmat
				subprocess.call(command, shell=True)		
			if args.method == 'ML' or args.method == 'B':
				command = "raxml -D -m PROTGAMMALG -s " + outfile + " -n " + os.path.splitext(os.path.basename(outfile))[0] + "_raxml.nwk -p 2233"
				subprocess.call(command, shell=True)		
			for o in range(1,boot+1):
				boot_col_index = col_index[0:cum_avg_len[i-1]]
				boot_col_index_with_repl = [random.choice(boot_col_index) for _ in range(cum_avg_len[0]+1)]

				for m in range(0,len(boot_col_index_with_repl)-1):
					if m == 0:
						bs_conc_aln = align_array[boot_col_index_with_repl[m],:]
					else:
						bs_conc_aln = np.column_stack((bs_conc_aln,align_array[boot_col_index_with_repl[m],:]))
			
				tmp_bs_align = MultipleSeqAlignment([])
				for seq in range(0,len(supermatrix)):
					myseq = bs_conc_aln[seq,:].tolist()
					myseqrecord = SeqRecord(Seq("".join(myseq),IUPAC.protein),id=supermatrix[seq].id)
					tmp_bs_align.append(myseqrecord)
					#tmp_bs_align.add_sequence(supermatrix[seq].id,"".join(myseq))

				bs_outfile = "unit_" + str(i) + "_sample_" + str(k) + "_rep_" + str(o) + ".bs"
				bs_out = open(bs_outfile,"w")
				AlignIO.write(tmp_bs_align,bs_out,"phylip")
				bs_out.close()
				if args.method == 'ML' or args.method == 'B':
					command = "raxml -D -m PROTGAMMALG -s " + bs_outfile + " -n " + os.path.splitext(os.path.basename(bs_outfile))[0] + "_raxml.nwk -p 2233"
					subprocess.call(command, shell=True)
					command = "cat RAxML_bestTree." + os.path.splitext(os.path.basename(bs_outfile))[0] + "_raxml.nwk >> unit_" + str(i) + "_sample_" + str(k) + "_raxml_rep.trees"
					subprocess.call(command, shell=True)
					command = "rm RAxML_*." + os.path.splitext(os.path.basename(bs_outfile))[0] + "_raxml.nwk"
					subprocess.call(command, shell=True)
				if args.method == 'ME' or args.method == 'B':
					command = "fastme -i " + bs_outfile + " -p -g 1.0 -s -n"
					subprocess.call(command, shell=True)
					command = "cat " + bs_outfile + "_fastme_tree.nwk >> unit_" + str(i) + "_sample_" + str(k) + "_rep.trees"
					subprocess.call(command, shell=True)
					command = "rm\t" + bs_outfile + "\t" + bs_outfile + "_fastme_tree.nwk\t" + bs_outfile + "_fastme_stat.txt"
					subprocess.call(command, shell=True)
if args.perc==True:
	for percent in range(10,100,10):
		counter = 0
		for aln in alns_list:
			counter+=1
			aln_len = aln.get_alignment_length()
			aln_array = np.array([list(rec) for rec in aln], np.character)
			#aln_array = aln_array.transpose()
			col_num=int(round(aln_len*percent/100))
			if counter == 1:
               	       		short_conc_aln_unit = aln_array[:,:col_num]
			else:
				short_conc_aln_unit = np.concatenate((short_conc_aln_unit,aln_array[:,:col_num]),axis=1)
				
			unit_init_aln_concat_align = MultipleSeqAlignment([])
			for u in range(0,len(aln)):
				single_unit = short_conc_aln_unit[u,:].tolist()
				single_unit_record = SeqRecord(Seq("".join(single_unit),IUPAC.protein),id=aln[u].id)
				unit_init_aln_concat_align.append(single_unit_record)
			outfile = str(percent) + "_percent_unit_" + str(counter) + "_init_alns_concatenated_aln.phylip"
			out = open(outfile,"w")
			AlignIO.write(unit_init_aln_concat_align,out,"phylip")
			out.close()

		"""	
		init_aln_concat_align = MultipleSeqAlignment([])
		for u in range(0,len(aln)):
			single_unit = short_conc_aln_unit[u,:].tolist()
			single_unit_record = SeqRecord(Seq("".join(single_unit),IUPAC.protein),id=aln[u].id)
			init_aln_concat_align.append(single_unit_record)
		outfile = str(percent) + "_percent_init_alns_concatenated_aln.phylip"
		out = open(outfile,"w")
               	AlignIO.write(init_aln_concat_align,out,"phylip")
               	out.close()
		"""

