#!/usr/bin/env perl

#
#Map Gene names to fasta format files (fasta + expresso templates)
#
#Usage: perl map_genenames_to_template.pl <UNIPROT-GENENAMES> <FASTA>
#
#



use strict;


my $line;
my $i;
my $k;
my $val;
my @genenames;
my @tmp;

my @store;
my @sorted_store;

open (GEN,"<","$ARGV[0]") or die "Cannot find $ARGV[0]\n";

@genenames = <GEN>;
chomp (@genenames);
close (GEN);

open (TEMP,"<","$ARGV[1]") or die "Cannot find $ARGV[1]\n";
$/ = '>';

foreach $i (@genenames)
{
    @store = ();
    @sorted_store = ();
    @tmp = split ("\t",$i);
    while ($line = <TEMP>)
    {
        chomp ($line);
        if ( $line =~ /$tmp[0].*\/([0-9]+)\-[0-9]+/ )
        {
            push (@store,$1);
        }
    }
    if ( scalar (@store) == 1 )
    {
        @sorted_store = sort {$a <=> $b} @store;
        for ($k=0;$k<=$#sorted_store;$k++)
        {
            seek (TEMP,0,0);
            while ($line = <TEMP>)
            {
                chomp ($line);
                if ( $line =~ /$tmp[0].*\/$sorted_store[$k]\-[0-9]+/ )
                {
                    $val = $k+1;
                    $line =~ s/$tmp[0].*\/$sorted_store[$k]\-[0-9]+/$tmp[1]\_$val/;
                    print ">$line";
                }
            }

            
        
        }
    }
    #}
    #else
    #{
    #    seek (TEMP,0,0);
    #    while ($line = <TEMP>)
    #    {
    #        chomp ($line);
    #        if ( $line =~ /$tmp[0].*\/([0-9]+)\-[0-9]+/ )
    #        {
    #            $line =~ s/$tmp[0].*\/([0-9]+)\-[0-9]+/$tmp[1]/;
    #            print ">$line";
    #        }
    #    }
    #}
    seek (TEMP,0,0);
}

close(TEMP);
