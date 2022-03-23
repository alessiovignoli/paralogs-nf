#!/usr/bin/env perl


###########################################################################
#
# Extract paralogs for different organisms from a PFAM file
# containing all the fasta sequences of a PFAM entry extracted from Pfam-A.fasta
# 
# Usage: perl extract_fasta_for_org_from_PFAM_fasta.pl <PFAM_FASTA> <ORG_IDS>
# 
# ORG_IDS: PFAM abbreviations
#
###########################################################################


use strict;
use File::Copy;

my $file;
my $line;
my @org_ids;
my $i;
my $filename;

open(ORG,"<","$ARGV[1]");
@org_ids = <ORG>;
chomp (@org_ids);
close(ORG);

$/ = '>';
open(FILE,"<","$ARGV[0]");
foreach $i (@org_ids)
{
	$ARGV[0] =~ /(\w+)\.fa.*/;
	$filename = $1.".".$i."_domain_sequences.fasta";
	open(OUT,">","$filename") or die("Cannot create the file $filename\n");
        while ($line = <FILE>)
        {
 	       	chomp $line;
 	       	if ( $line =~ /\_$i/ )
 	       	{
	       		print OUT ">$line";
	       	}					
        }
	close(OUT);
	seek (FILE,0,0);
}
close(FILE);
