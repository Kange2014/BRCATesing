#!/usr/bin/perl -w

use strict;

use Bio::SeqIO;

my %brca1 =();
my $str = Bio::SeqIO->new(-file=> "../BRCA1_exons.fa",-format => 'fasta');
while(my $seq = $str->next_seq()){
	my @ids = split /\:/,$seq->display_id();
	my @locs = split /\,/,$ids[1];
	if($ids[0] =~ /([0-9]+)\-/){ 
		if(exists $brca1{$1}) { ${$brca1{$1}}[1] = $locs[1]; }
		else{$brca1{$1} = [$locs[0],$locs[1]]; }
	}
	else { $brca1{$ids[0]} = [$locs[0],$locs[1]]; }
}

my %brca2 = ();
$str = Bio::SeqIO->new(-file=> "../BRCA2_exons.fa",-format => 'fasta');
while(my $seq = $str->next_seq()){
	my @ids = split /\:/,$seq->display_id();
	my @locs = split /\,/,$ids[1];
	if($ids[0] =~ /([0-9]+)\-/){ 
		if(exists $brca2{$1}) { ${$brca2{$1}}[1] = $locs[1]; }
		else{$brca2{$1} = [$locs[0],$locs[1]]; }
	}
	else { 
		if($ids[0] == 6){
			${$brca2{5}}[1] = $locs[1];
		}else{
			$brca2{$ids[0]} = [$locs[0],$locs[1]]; 
		}
	}

}

my %map1 = ();
for(my $i = 1; $i <= 27; $i++){
	if($i <=2) { $map1{$i} = $i+1; }
	if($i <= 9 && $i > 2) { $map1{$i} = $i+2; }
	if($i <=14 && $i > 9) {$map1{$i} = 11; }
	if($i > 14){ $map1{$i} = $i - 3;}
}

my %map2 = ();
for(my $i = 1; $i <= 35; $i++){
	if($i <=4 || $i == 9 || $i == 10) { $map2{$i} = $i+1; }
	if($i <= 8 && $i > 4) { $map2{$i} = $i+2; }
	if($i <=18 && $i >= 11) {$map2{$i} = 11; }
	if($i > 18 && $i <=34){ $map2{$i} = $i - 7;}
	if($i == 35) { $map2{$i} = 27;}
}

open FILE,"all.seq.blast.out";
open OUT,">all.seq.analysis.txt";
while(<FILE>){
	chomp;
	if($. == 1){
		print OUT $_. "\t"."Exon.ID"."\tExon.start"."\tExon.end\tUpstream\tDownstream\n";
		next;
	}
	print OUT $_;
	my @columns = split /\t/;
	my @ids = split/\-/,$columns[0];
	if($ids[0] eq "B1"){
		print OUT "\t"."$map1{$ids[1]}\t".join("\t", @{$brca1{$map1{$ids[1]}}})."\t";
		my ($up,$down);
		if($columns[8] < $columns[9]) {
			$up = ${$brca1{$map1{$ids[1]}}}[0] - $columns[8];
			$down = $columns[9] - ${$brca1{$map1{$ids[1]}}}[1];
		}
		else{
			$up = ${$brca1{$map1{$ids[1]}}}[0] - $columns[9];
			$down = $columns[8] - ${$brca1{$map1{$ids[1]}}}[1];
		}
		print OUT "$up\t$down\n";
	}
	if($ids[0] eq "B2"){
		print OUT "\t"."$map2{$ids[1]}\t".join("\t", @{$brca2{$map2{$ids[1]}}})."\t";
		my ($up,$down);
		if($columns[8] < $columns[9]) {
			$up = ${$brca2{$map2{$ids[1]}}}[0] - $columns[8];
			$down = $columns[9] - ${$brca2{$map2{$ids[1]}}}[1];
		}
		else{
			$up = ${$brca2{$map2{$ids[1]}}}[0] - $columns[9];
			$down = $columns[8] - ${$brca2{$map2{$ids[1]}}}[1];
		}
		print OUT "$up\t$down\n";
	}
}
close FILE;
close OUT;
