#!/usr/bin/env perl


###########################################################################
#
# 
# 
# 
# Usage: perl extract_seq_per_species_sim.pl <ORG_IDS> <MSA file in fasta> <-after / -prior>
# 
# 
#
###########################################################################


use strict;

my $line;
my @org_ids;
my $i;
my $filename = $ARGV[1];
my $output;


open(ORG,"<","$ARGV[0]");
@org_ids = <ORG>;
chomp (@org_ids);
close(ORG);

my $aln = $ARGV[2] or die "Cannot find $ARGV[2]\n";

$/ = '>';

open(IN,"<","$filename") or die "Cannot find the file: $filename\n";

foreach $i (@org_ids)
{
	if ( $aln eq '-prior' ) {
		$output = $filename."_".$i."_for_full_aln.fasta_aln";		
		open(OUT,">","$output") or die "Cannot create the file: $output\n";
		while ($line = <IN>)
		{
			chomp ($line);
			if ( $line =~ /^$i\_\d+/ )
			{
				print OUT ">$line";
			}
		}
		seek (IN,0,0);	
		close(OUT);
	}
	else {
		$output = $filename."_".$i.".fasta_aln";		
		open(OUT,">","$output") or die "Cannot create the file: $output\n";
		while ($line = <IN>)
		{
			chomp ($line);
			if ( $line =~ /^$i\_\d+/ )
			{
				$line =~ s/$i/Seq/;
				print OUT ">$line";
			}
		}
		seek (IN,0,0);	
		close(OUT);
	}
}
close(IN);
