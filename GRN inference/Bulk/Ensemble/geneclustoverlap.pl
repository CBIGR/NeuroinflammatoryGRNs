#!/usr/bin/env perl -w
# This program calculates the overlap in regulators between genes in the clusters
# As inputfile: kmedclusters.txt first column is gene name, second column is cluster nr

#perl geneclustoverlap.pl kmedclusters.txt kmedregoverlap


use strict;
open (R,"<gene_reg_colAD.txt") or die $!;
open (IN, "<$ARGV[0]") or die $!;
die "won't overwrite existing file" if -e $ARGV[1];
open (OUT, ">$ARGV[1]") or die $!;


my %clust;
while (<IN>){
	chomp;
	my ($g, $c) = (split /\t/) [0,1];
	$clust{$c}{$g}=1;
	}
	
	
my %GR;
#hash of arrays: for each gene its set of ranked regulators
my @rows = <R>;
foreach my $r (@rows) {
	chomp $r;
	my @R = split /\t/, $r;
	my $gene = shift @R;
	my @regs = @R;
	$GR{$gene}=[@R];
	}

my %overlap;
my %gperc;
my %rankR;
my %countR;
my %tot;
foreach my $c (sort keys %clust){
	foreach my $g (sort keys %{$clust{$c}}){
	$gperc{$c}++;
	#number of genes per cluster
		foreach my $r (sort @{$GR{$g}}){
			my( $index )= grep { ${$GR{$g}}[$_] eq $r } 0..$#{$GR{$g}};
			#order/rank of the regulator in the list of predicted regulators of a gene
			$rankR{$c}{$r} = $rankR{$c}{$r}+$index+1;
			#sum of all the ranks of the regulator in a cluster
			$countR{$c}{$r}++;
			#number of genes a regulator targets in the cluster
			if (exists $overlap{$c}{$r}){
			$overlap{$c}{$r}++;
			}else{
			$overlap{$c}{$r}=1;
			#number of genes the regulator regulates in the cluster
			}
			}
			}
			}

foreach my $c (sort keys %overlap){
	foreach my $r (sort {$overlap{$c}{$b} <=> $overlap{$c}{$a}} keys %{$overlap{$c}}){
		$tot{$c}{$r}=$rankR{$c}{$r}/$countR{$c}{$r};
		}
		}
#average rank of the regulator per gene in the cluster = > calculated as sum of all ranks for a regulator among all regulators for a cluster, in its interactions with the genes of a cluster divided by the number of genes a regulator targets in a cluster

print OUT "#cluster\t#gperc\t#reg(count_rank)\n";	
foreach my $c (sort keys %overlap){
	my $g = $gperc{$c};
	print OUT "$c\t$g\t";
	my $count=1;
	foreach my $r (sort {$overlap{$c}{$b} <=> $overlap{$c}{$a} || $tot{$c}{$a} <=> $tot{$c}{$b}} keys %{$overlap{$c}}){
		my $C = $overlap{$c}{$r};
		if ($C>2){
		#regs should be shared by more than 2 genes in the cluster
		my $R=int($tot{$c}{$r}+0.5);
		#average rank of the regulator per gene in the cluster
		my $p = $C/$g*100;
		my $P = int($p+0.5);
		my $s = $conv{$r};
		if (($P >= 50)&& $count <=10){
		print OUT "$s";
		print OUT "_";
		#print OUT "$s($P%_$R)_";
		#print OUT "$r($P%_$R)\t";
		$count++;
		}
		}
		}
		print OUT "\n";
		}
		


