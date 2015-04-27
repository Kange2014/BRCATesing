#!/usr/bin/perl -w

use strict;

### This script is used for selecting primers for each exon; that is,
### each primer pair should cover the whole exon region, and if possible,
### for each exon, selecting as much as possible primer pairs that don't
### overlap with each other.

my $input = shift @ARGV;
my $output = shift @ARGV;

my %hash = ();
my @primer = ();

open FILE,"<",$input;
open OUT,">",$output;
while(<FILE>){
	if($. == 1){
		print OUT $_;
		next;
	}
	my @columns1 = split /\t/;
	my $line1 = $_;
	
	my $line2 = <FILE>;
	my @columns2 = split /\t/,$line2;
	
	if(exists $hash{$columns1[0]}){
		my $flag = 0;
		my $len = keys %{$hash{$columns1[0]}};
		foreach my $id (keys %{$hash{$columns1[0]}}){
			if($columns1[2] > $hash{$columns1[0]}{$id}[0] + $hash{$columns1[0]}{$id}[1] - 1 || $columns1[2] + $columns1[3] - 1 < $hash{$columns1[0]}{$id}[0] ||
			   $columns2[2] < $hash{$columns1[0]}{$id}[2] - $hash{$columns1[0]}{$id}[3] + 1 || $columns2[2] - $columns2[3] + 1 > $hash{$columns1[0]}{$id}[2]
			)
			{ next;}
			else
			{ $flag = 1;}
		}
		if($flag == 0){
			push @primer,($line1,$line2);
			$len += 1;
			$hash{$columns1[0]}{$len} = [$columns1[2],$columns1[3],$columns2[2],$columns2[3]];
		}
	}
	else{
		$hash{$columns1[0]}{1} = [$columns1[2],$columns1[3],$columns2[2],$columns2[3]];	
		push @primer,($line1,$line2);
	}
}

foreach my $line(@primer){
	print OUT $line;
}
close OUT;
