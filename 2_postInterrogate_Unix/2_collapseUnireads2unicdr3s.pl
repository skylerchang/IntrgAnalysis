#!/bin/perl
#collapseUnireads2unicdr3s.pl

use strict; use warnings;

#my $targetFolder = "/Users/SKeller/Documents/Projects/feTR_HTS/Lymphomas/all/Clntab_test";
my $targetFolder = "/Users/SKeller/Documents/Sequencing/Runs/k9MultiLoci3-73454401/Clntab_2018-04-25/";


my @loci = qw(IGH TRA TRB TRD TRG noju);
my $cdr3Col = 14;   #sequence.JUNCTION.aa seq MINUS 1 ***** array position
my %stats;

opendir (DIR, "$targetFolder/Data");
my @dirs = grep {! /^\./} readdir(DIR);

foreach my $dir (@dirs){
    
    print "\nParsing directory $dir ...\n";
    foreach my $locus (@loci){
        print "Parsing locus $locus ...\n";
        my %data;
        my $infile =  "$targetFolder/Data/$dir/$dir" . "_$locus.clntab\n";
        open IN, $infile or die "Error opening infile at $infile";
        while (<IN>){
            chomp;
            my @fields = split /\t/, $_;
            my $cdr3Seq = $fields[$cdr3Col];
            push @{$data{$cdr3Seq}}, \@fields;
        }
        close IN;
        
        print "Printing compressed clntab file\n";
        my $outfile =  "$targetFolder/Data/$dir/$dir" . "_$locus" . "_uniqueCdr3s.clntab\n";
        open OUT, ">$outfile" or die "Outfile could not be opened at '$outfile'";
        foreach my $cdr3Seq (keys %data){
            my $unireadsPerCdr3 = scalar @{$data{$cdr3Seq}};
            my $counter = 0;
            foreach my $uniread (@{$data{$cdr3Seq}}){
#                print "$uniread->[26]\n";
                $counter += $uniread->[26];   #field with 'size'
            }
            $stats{$dir}->{$locus}->{$unireadsPerCdr3}++;
            my @s = @{$data{$cdr3Seq}->[0]};
            my $line = join "\t", @s[0 .. 25], $counter, @s[26 .. 32];
            print OUT "$line\n";
        }
        close OUT;
    }
}

print "Printing result files ...\n";
mkdir "$targetFolder/Results/UnireadsPerUniqueCDR3";
foreach my $dir (sort keys %stats){
    mkdir "$targetFolder/Results/UnireadsPerUniqueCDR3/$dir";
    foreach my $locus (sort keys %{$stats{$dir}}){
        my $outfile =  "$targetFolder/Results/UnireadsPerUniqueCDR3/$dir/$dir" . "_$locus" . "_uniqueCdr3s.txt\n";
        open OUT, ">$outfile" or die "Outfile could not be opened at '$outfile'";
        
        print OUT "unireadsPerUniqueCdr3\tfrequency\n";
        foreach my $count (sort {$a <=> $b} keys %{$stats{$dir}->{$locus}}){
            print OUT "$count\t$stats{$dir}->{$locus}->{$count}\n";
        }
        close OUT;
    }
}

print "\nDone";



