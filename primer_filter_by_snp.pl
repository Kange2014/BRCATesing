#!/usr/bin/perl -w

use strict;

use Text::CSV;

my %brca1_snp = ();
my %brca2_snp = ();

#### dbSNP 142 for a single gene
#open FILE1,"BRCA1_672_gene_report.txt";
#while(<FILE1>){
#	next if $. <= 3;
#	next unless /\t/;
#	my @columns = split /\t/;
#	my @HGVS = split /\;/,$columns[2];
#	foreach my $id(@HGVS){
#		if($id =~ /NG_005905/){ #### BRCA1 genomic sequence accession number.
#			my @array = split /\:/,$id;
#			if($array[1] =~ /g\.([0-9]+)[A-Z]/){ $brca1_snp{$1} = $columns[7] if $columns[7] ne "";}
#			last;
#		}
#	}
#}
#close FILE1;

#open FILE2,"BRCA2_675_gene_report.txt";
#while(<FILE2>){
#	next if $. <= 3;
#	next unless /\t/;
#	my @columns = split /\t/;
#	my @HGVS = split /\;/,$columns[2];
#	foreach my $id(@HGVS){
#		if($id =~ /NG_012772/){ #### BRCA2 genomic sequence accession number.
#			my @array = split /\:/,$id;
#			if($array[1] =~ /g\.([0-9]+)[A-Z]/){ $brca2_snp{$1} = $columns[7] if $columns[7] ne "";}
#			last;
#		}
#	}
#}
#close FILE2;

### GRCh38
### BRCA1 genomic locs are Chromosome 17: 43,024,295-43,217,981 reverse strand
### BRCA2 genomic locs are Chromosome 13: 32,310,480-32,401,672 forward strand

### Due to a somewhat catastrophic hardware failure during our production cycle 
### for Ensembl release 79, we have only been able to release human dbSNP 142 
### incorporating 1000 Genomes phase 3 data on our GRCh37-based services. Pls. 
### refer to http://www.ensembl.info/blog/2015/04/09/dbsnp-142-and-1000-genomes-phase-3/

### GRCh37 
### BRCA1 genomic locs are Chromosome 17: 41,176,312-41,370,000 reverse strand
### BRCA2 genomic locs are Chromosome 13: 32,884,617-32,975,809 forward strand


#open FILE,"BRCA_Ensembl_Variation.txt";
open FILE,"BRCA_Ensembl_Variation79_GRCh37.txt";
while(<FILE>){
	next if $. == 1;
	chomp;
	my @columns = split /\t/;
	my $tmp;
	if($columns[3] > $columns[4]){
		$tmp = $columns[4];
		$columns[4] = $columns[3];
		$columns[3] = $tmp;
	}
	if($columns[2] == 17){		### BRCA1 gene
		for(my $i = $columns[3]; $i <= $columns[4]; $i++){
			#my $repos = 43217981 - ($columns[3] + $i - 1)  + 1;
			my $repos = 41370000 - ($columns[3] + $i - 1)  + 1;
			$brca1_snp{$repos} = $columns[6] if $columns[6];
			$brca1_snp{$repos} = 0 unless $columns[6];
		}
	}
	
	if($columns[2] == 13){		### BRCA2 gene
		for(my $i = $columns[3]; $i <= $columns[4]; $i++){
			#my $repos = $columns[3] + $i - 1 - 32310480 + 1;
			my $repos = $columns[3] + $i - 1 - 32884617 + 1;
			$brca2_snp{$repos} = $columns[6] if $columns[6];
			$brca2_snp{$repos} = 0 unless $columns[6];
		}
	}
}
close FILE;


filter_by_snp("BRCA1_primers.csv","BRCA1_filter_by_SNP.txt",\%brca1_snp);
filter_by_snp("BRCA2_primers.csv","BRCA2_filter_by_SNP.txt",\%brca2_snp);

sub filter_by_snp{
	my ($input,$output,$snp) = @_;
	
	my $junction_size = 200;
	
	my $csv = Text::CSV->new();
	open (CSV, "<", $input) or die $!;
	open OUT,">$output";
	
while(<CSV>){
	if ($. == 1) {
		next unless $csv->parse($_);
		my @columns  = $csv->fields();
		print OUT join("\t",@columns)."\n";
		next;
	}
	my $flag = 0;
	
	next unless $csv->parse($_);
  	my @columns1 = $csv->fields();
	chomp;

	my @ids = split /\:/,$columns1[0];
	my ($exon,$sub) = (0,0);
	if($ids[0] =~ /\-/){ ($exon,$sub) = split /\-/,$ids[0]; }
	my ($start,$end) = split /\,/,$ids[1];
	my $primer_start;
	if ($sub == 0 || $sub == 1){
		$primer_start = $start - $junction_size + $columns1[2] - 1;
	}
	else{
		$primer_start = $start + $columns1[2] - 1;
	}
	$columns1[2] = $primer_start;
	for(my $i = 1; $i <= $columns1[3]; $i++){
		my $position = $primer_start + $i - 1;
		#if(exists $snp->{$position} && $snp->{$position} >= 0.005){ $flag = 1;}
		if(exists $snp->{$position}){ $flag = 1;}
	}
	
	$_ = <CSV>;
	
	next unless $csv->parse($_);
  	my @columns2 = $csv->fields();
	chomp;
	
	@ids = split /\:/,$columns2[0];
	
	($exon,$sub) = (0,0);
	if($ids[0] =~ /\-/){ ($exon,$sub) = split /\-/,$ids[0]; }
	($start,$end) = split /\,/,$ids[1];
	if ($sub == 0 || $sub == 1){
		$primer_start = $start - $junction_size + $columns2[2] - 1;
	}
	else{
		$primer_start = $start + $columns2[2] - 1;
	}
	$columns2[2] = $primer_start;
	for(my $i = $columns2[3]; $i >= 1; $i--){
		my $position = $primer_start + $i - 1;
		#if(exists $snp->{$position} && $snp->{$position} >= 0.005){ $flag = 1;}
		if(exists $snp->{$position} ){ $flag = 1;}
	}
	
	if($flag == 0){
		print OUT join("\t",@columns1)."\n";
		print OUT join("\t",@columns2)."\n";
	}
}
close CSV;
close OUT;

}