#!/usr/bin/env perl -w
# This program calculates the Jaccard index between all the genes, based on the predicted regulators
# As inputfile: first column is gene name, next columns ordered according to rank of inferred regulators 
#perl geneclust_jaccard.pl geneclust_out overlap_geneclust

use strict;
open (FILE, "<gene_reg_colAD.txt") or die $!;
die "won't overwrite existing file" if -e $ARGV[0];
open (OUT, ">$ARGV[0]") or die $!;
#jaccard index
die "won't overwrite existing file" if -e $ARGV[1];
open (OVERLAP, ">$ARGV[1]") or die $!;

my %GR;
my %gene;
my @rows = <FILE>;
foreach my $r (@rows) {
	chomp $r;
	my @R = split /\t/, $r;
	my $gene = shift @R;
	my @regs = @R;
	$GR{$gene}=[@R];
	$gene{$gene}=1;
	}
print "read in OK\t";
my %overlap;
my %jaccard;
foreach my $g (sort keys %gene){
	my @R = @{$GR{$g}};
	my $l1= scalar @R;
	my %R = map {$_=>1} @R;
	foreach my $k (sort keys %GR){
		unless (($g eq $k)|(exists $overlap{$k}{$g})|(exists $overlap{$g}{$k})){
		my @K = @{$GR{$k}};
		my $l2=scalar @K;
		my $count=0;
		foreach my $kk (@K){
			if (exists $R{$kk}){
			$count++;
			push (@{$overlap{$g}{$k}},$kk);
			}
			}
		$jaccard{$g}{$k}=$count/($l1+$l2-$count);
		}
		}
		}	
		
print "jaccard calculation OK\t";
#calculating the Jaccard index: the size of the intersection divided by the size of the union
	
#print OUT "Jaccard index\n";
foreach my $a (sort keys %jaccard){
	foreach my $b (sort keys %{$jaccard{$a}}){
	my $j = $jaccard{$a}{$b};
	if ($j){
	my @o = @{$overlap{$a}{$b}};
	print OUT "$a\t$b\t$j\n";
	my $overlap = join "_",@o;
	print OVERLAP "$a\t$b\t$overlap\n";
	}
	}
	}

