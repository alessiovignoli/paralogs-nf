#!/usr/bin/env perl


###########################################################################
#
# Extract sequences from each species after intersection
# 
# 
# Usage: perl extract_intersecting_seq.pl <INTERSECTING_GENE_IDS> <ORG_IDS> <BASE_DIR> < -prior | -after > <PFAM_ID>
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
my $base_dir = $ARGV[2] or die "Cannot find $ARGV[2]\n";

my $aln = $ARGV[3] or die "Cannot find $ARGV[3]\n";

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
	$filename = $ARGV[4].".$i"."_domain_sequences_intersect_with_ref_after_OMA.fasta";
	$output = $ARGV[4].".$i"."_domain_sequences_after_intersection.fasta";
			
	open(IN,"<","$filename") or die "Cannot find the file: $filename\n";
	open(OUT,">","$output");
	
	if ( $aln eq '-prior' )
	{
		foreach $m (@gene_ids)
		{
			while ($line = <IN>)
			{
				chomp ($line);
				if ( $line =~ /^$m\n/ )
				{
					if ($i ne 'HUMAN' ) {$line =~ s/$m\n/$m\_$i\n/;}
					print OUT ">$line";
					last;
				}
			}
			seek (IN,0,0);
		}
	}
	else
	{	
		foreach $m (@gene_ids)
		{
			while ($line = <IN>)
			{
				chomp $line;
				if ( $line =~ /^$m\n/ )
				{
					print OUT ">$line";
					last;
				}
			}
			seek (IN,0,0);						
		}
	}
	close(IN);
	close(OUT);
}

