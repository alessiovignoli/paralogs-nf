#!/usr/bin/env perl


###########################################################################
#
# 
# 
# 
# Usage: perl extract_1_to_1_orthologs_for_each_species_OMA.pl <UNIPROT_OMA-GROUPS_REF> <ORG_IDS> <UNIPROT_ID-GENENAMES> <FAM_ID> <REF_SPECIES>
# 
# ORG_IDS: PFAM abbreviations
#
###########################################################################


use strict;
use Cwd qw(getcwd);


my $line;
my $line2;
my $ensline;

my @org_ids;
my @uni_oma_ids;
my @uni_oma_groups;
my @genenames;

my $i;
my $m;
my $k;
my $found;

my $ref_uni;
my $org_oma;
my $org_uni;
my $ref_genename;
my $org_ense;
my $oma_group;

my $ens_out;
my @tmp;

my $output;
my $base_dir = getcwd;
my $path_to_uni_oma_ids = "/users/cn/abaltzis/projects/1_paralogs/paraconc-nf/OMA/oma-uniprot.txt";
my $path_to_ensembl_oma_ids = "/users/cn/abaltzis/projects/1_paralogs/paraconc-nf/OMA/oma-ensembl.txt";
#my $path_to_pairs = "/users/cn/abaltzis/projects/1_paralogs/paraconc-nf/OMA/HUMAN_oma-pairs-1-to-1.txt";
my $path_to_pairs = "/users/cn/abaltzis/projects/1_paralogs/paraconc-nf/OMA/" . $ARGV[4] . "_oma-pairs-1-to-1.txt";


open(ORG,"<","$ARGV[1]") or die "Cannot open $ARGV[1]\n";
@org_ids = <ORG>;
chomp (@org_ids);
close(ORG);

open(OMA,"<","$ARGV[0]") or die "Cannot open $ARGV[0]\n";
@uni_oma_groups = <OMA>;
chomp (@uni_oma_groups);
close(OMA);

open (GENENAMES,"<","$ARGV[2]") or die "Cannot open $ARGV[2]\n";
@genenames = <GENENAMES>;
chomp (@genenames);
close(GENENAMES);


foreach $i (@org_ids)
{
	$output = $ARGV[3].".".$i."_domain_sequences_intersect_with_ref.tab";
        open(OUT,'>',"$output") or die "Cannot open $output\n";
	foreach $m (@uni_oma_groups)
	{
		$found = 0;
		if ( $i eq 'RAT') {$i = 'RATNO';}
		if ( $m =~ /^(\w+)\s+(\w+)\s*.*\s($i\w+)\s+/ )
		{
			$ref_uni = $1;
			$oma_group = $2;
			$org_oma = $3;
			
			foreach $k (@genenames)
			{
				if ( $k =~ /^$ref_uni\s+(\w+)/ )
				{
					$ref_genename = $1;
					last;	
				}
			}
			$line = `grep "$org_oma" $path_to_uni_oma_ids`;
			$line2 = `grep "$org_oma.*$oma_group" $path_to_pairs`; 
			if ( $line =~ /$org_oma\s+(\w+)/ && $line2 )
			{
				$org_uni = $1;
				$found = 1;
				print OUT "$org_uni\t$ref_genename\t$ref_uni\t$org_oma\n";
			}
			if ( $found != 1 )
			{
				$ensline = `grep "$org_oma" $path_to_ensembl_oma_ids`;
				if ( $ensline =~ /$org_oma\s+(\w+)/ && $line2 )
				{
					$org_ense = $1;
					$ens_out = ensembl_to_uniprot($org_ense);
					@tmp = split ("\t", $ens_out);
					$org_uni = $tmp[1];
					if ($org_uni) {print OUT "$org_uni\t$ref_genename\t$ref_uni\t$org_oma\n";}
				}					
				
			}
		}		
	}
	close(OUT);
}


sub ensembl_to_uniprot {

use strict;
use warnings;
use LWP::UserAgent;

my $base = 'https://www.uniprot.org';
my $tool = 'uploadlists';
my $ids = $_[0];

my $out;
my $send;

my $params = {
  from => 'ENSEMBL_ID',
  to => 'ID',
  format => 'tab',
  query => "$ids"
};

my $contact = 'athanasios.baltzis@crg.eu'; # Please set a contact email address here to help us debug in case of problems (see https://www.uniprot.org/help/privacy).
my $agent = LWP::UserAgent->new(agent => "libwww-perl $contact");
push @{$agent->requests_redirectable}, 'POST';

my $response = $agent->post("$base/$tool/", $params);

while (my $wait = $response->header('Retry-After')) {
  print STDERR "Waiting ($wait)...\n";
  sleep $wait;
  $response = $agent->get($response->base);
}

if ($response->is_success)
{
	$out = $response->content;
	$out =~ /\n(.*\t.*)$/;
	$send = $1;
 	return $send;
}
}
