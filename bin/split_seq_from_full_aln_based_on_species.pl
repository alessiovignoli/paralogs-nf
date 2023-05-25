#!/usr/bin/env perl


###########################################################################
#
# Split sequences from full alignment based on each species
# 
# 
# Usage: perl split_seq_from_full_aln_based_on_species.pl <INTERSECTING_GENE_IDS> <ORG_IDS> <PFAM_ID>
# 
# ORG_IDS: PFAM abbreviations
#
###########################################################################


use strict;

my $line;
my @org_ids;
my @gene_ids;
my $i;
my $filename;
my $m;
my $output;

open(ORG,"<","$ARGV[1]");
@org_ids = <ORG>;
chomp (@org_ids);
close(ORG);

open(GENES,"<","$ARGV[0]");
@gene_ids = <GENES>;
chomp (@gene_ids);
close(GENES);

$/ = '>';

foreach $i (@org_ids)
{
	$filename = $ARGV[2].".domain_sequences_prior_after_intersection_full.fasta_aln";
	open(IN,"<","$filename") or die "Cannot find the file: $filename\n";
	
	$output = $ARGV[2].".$i"."_domain_sequences_prior_after_intersection.fasta_aln";
	open(OUT,">","$output");
	while ($line = <IN>)
	{
		chomp ($line);
		if ( $line =~ /\_$i\_/ )
		{
			if ($i ne 'MOUSE' ) {$line =~ s/\_$i\_/_MOUSE_/;}
			print OUT ">$line";
		}
	}
	seek (IN,0,0);					
	close(IN);
	close(OUT);
}