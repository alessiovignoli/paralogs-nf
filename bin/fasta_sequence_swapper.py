#!/usr/bin/env python3

import argparse

from random import sample

def get_args():

    """get the arguments when using from the commandline"""

    parser = argparse.ArgumentParser(description="Launch fasta sequence swapper. It keeps the header as they are preserving they're order to, but it swaps a given number of sequences around the file randomly.")
    parser.add_argument("-i", "--input", type=str, required=True, metavar="FILE", help="Path to input fasta file. Sequences have to be on one line.")
    parser.add_argument("-o", "--output", type=str, required=True, metavar="FILE", help="Path to output file to write the new swapped fasta.")
    parser.add_argument("-n", "--nseq", type=int, required=True, nargs='?', const=None, default=None, metavar="NUM_SEQ_IN_FILE", help="The number of sequences present in the file. This has to be known a priori so that sequences can be swapped randomly by drawing numbers from 1 to NUM_SEQ_IN_FILE.")
    parser.add_argument("-s", "--swaps", type=int, required=True, nargs='?', const=None, default=None, metavar="NUM_OF_SWAPS", help="The number of sequences to swap around. Swaps are made so that they do not reverse each other. If set to 5, then 5 sequences are selected and swapped between themselves.")
    
    args = parser.parse_args()

    return args


def SeqSwapper(input_fasta : str,
               output_fasta : str,
               total_seq : int,
               num_swaps : int
               ) -> None:

    # check all inputs are what they should be
    if total_seq <= 1:
        raise ValueError("there must be at least two sequences in the file")
    if num_swaps <= 1:
        raise ValueError("number of swaps should be greater or equal than 2. 2 sequences is the minimum so that they can be swapped between themselves.")
    if num_swaps > total_seq:
        raise ValueError("number of swaps can not be greater than the number of sequence in the file")
    
    # randomly draw the sequences indices to be swapped
    index_swap = sample(range(total_seq), num_swaps)

    # randomly decide the new indices for the sequences to be swapped and create a mapping
    swap_index_mapping = None
    while True:
        shuffled_index = sample(index_swap, num_swaps)
        if all(a != b for a, b in zip(index_swap, shuffled_index)):
            swap_index_mapping = dict(zip(index_swap, shuffled_index))
            break

    # open file in input and read it the first time and extract 
    # the sequences that have to be swapped and save them to buffer dict with new index as key
    swap_seq = {}
    with open(input_fasta, 'r') as infasta:
        i = 0
        for first_iter_line in infasta:
            if first_iter_line[0] != ">":
                if i in index_swap:
                    swap_seq[i] = first_iter_line.rstrip()
                i += 1

    # write to output all header and non swapped sequences as they are and
    # the swapped seq in the new positions
    with open(input_fasta, 'r') as infasta, open(output_fasta, 'w') as outfasta:
        j = 0
        for second_iter_line in infasta:
            if second_iter_line[0] == ">":
                outfasta.write(second_iter_line)
            else:
                if j in index_swap:
                    swapped_line = swap_seq[swap_index_mapping[j]] + "\n"
                    outfasta.write(swapped_line)
                else:
                    outfasta.write(second_iter_line)
                j += 1


    print(swap_index_mapping)
    print(swap_seq)



def main(input : str,
         output : str,
         total : int,
         swaps : int
         ) -> None:
    
    SeqSwapper(input, output, total, swaps)


if __name__ == "__main__":
    args = get_args()
    main(args.input,
         args.output,
         args.nseq,
         args.swaps)