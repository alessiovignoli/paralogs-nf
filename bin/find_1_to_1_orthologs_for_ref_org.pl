#!/usr/bin/env perl

##############################################
#
#Find and print oma-groups for uniprot ids 
#
#Usage: perl find_1_to_1_orthologs.pl <FASTA HEADERS FILE>
#
##############################################


use strict;


my $ref_fasta_headers = $ARGV[0];
my $oma_groups = "/users/cn/abaltzis/projects/1_paralogs/paraconc-nf/OMA/oma-groups.txt";
my $oma_uniprot = "/users/cn/abaltzis/projects/1_paralogs/paraconc-nf/OMA/oma-uniprot.txt";

my $output;
my $oma;
my $i;
my @fasta_headers;

my $query;

open (REF,"<","$ref_fasta_headers") or die "Cannot find fasta file\n";

@fasta_headers = <REF>;
chomp (@fasta_headers);
close (REF);

open (OMA,"<","$oma_groups") or die "Cannot find oma-groups file\n";


foreach $i (@fasta_headers)
{
    $output = `grep "\t$i" $oma_uniprot | cut -f1`;
    chomp ($output);
    while (<OMA>)
    {
        if ( ($output) && $_ =~ /$output/ )
        {
            print "$i\t$_";
            last;
        }
    
    } 
    
    seek (OMA,0,0);

}

close (OMA);
