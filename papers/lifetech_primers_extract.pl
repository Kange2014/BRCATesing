#!/usr/bin/perl -w

use strict;

use Bio::Location::SplitLocationI;
use Bio::SeqIO;
use Bio::Seq;
use Bio::Perl;

open FILE,"lifetech_brca_primers_locs.txt";
open OUT,">lifetech_brca_primers.txt";
print OUT "Gene\tExon\tStart\tEnd\tFW\tRE\tAmplicon_length\tLocation\n";

### the brca 1 & 2 whole genomic sequence is large
### and the corresponding regions for brca1 is 92501..173689
### for brca2 is 5001..89193
### the following files are regions for these two genes

my @chr13 = read_all_sequences("BRCA2-1.gb",'Genbank');
my @chr17 = read_all_sequences("BRCA1-1.gb",'Genbank');

my $seq13 = $chr13[0]; 
my $seq17 = $chr17[0];

my %brca1 = ();
my %brca2 = ();
my @features = $seq17->all_SeqFeatures();
foreach my $feat (@features){
	if($feat->primary_tag eq 'CDS' ){
		my @value=$feat->each_tag_value('gene');
		if($value[0] eq "BRCA1"){
			if ($feat->location->isa('Bio::Location::SplitLocationI')){
				my $count = 1;
				foreach my $loc ($feat->location->sub_Location){
					$brca1{$count} = $loc->start.",".$loc->end;
					$count++;
				}
			}
			last;
		}
	}
}



@features = $seq13->all_SeqFeatures();
foreach my $feat (@features){
	if($feat->primary_tag eq 'CDS' ){
		my @value=$feat->each_tag_value('gene');
		if($value[0] eq "BRCA2"){
			if ($feat->location->isa('Bio::Location::SplitLocationI')){
				my $count = 1;
				foreach my $loc ($feat->location->sub_Location){
					$brca2{$count} = $loc->start.",".$loc->end;
					$count++;
				}
			}
			last;
		}
	}
}

my $seq_tmp = Bio::Seq->new(-id => '00001',-seq=>'AAAAAA');


my @brca1_locs = (41196312,41277500); ### the total length of chr17 is 81195210 bp.
my @brca2_locs = (32889617,32973809);

my $gene;
while(<FILE>){
	if(/^Gene/){
		$gene = <FILE>;
		chomp($gene);
		$gene =~ s/\r//;
		print OUT $gene."\t";
	}
	if(/^Location/){
		my $start = <FILE>;
		my $end = <FILE>;
		chomp($start);
		chomp($end);
		$start =~ s/\r//;
		$end =~ s/\r//;
		$start =~ s/Chr.1[37]: //;
		$end =~ s/-//;
		
		
		if($gene eq "BRCA1"){
			$start -= 41196312 + 1;
			$end -= 41196312 + 1;
			
			foreach my $id(keys %brca1){
				my ($first,$second) = split /,/,$brca1{$id};
				print $brca1{$id},"\n";
				if($start >= $first && $end <= $second){
					print OUT $id."\t".$first."\t".$second."\t";
					last;
				}
			}
			$seq_tmp->seq($seq17->subseq($end-24,$end));
			print OUT $seq17->subseq($start,$start+24)."\t".$seq_tmp->revcom()->seq()."\t";
		}
		
		if($gene eq "BRCA2"){
			$start -= 32889617 + 1;
			$end -= 32889617 + 1;
		
			foreach my $id(keys %brca2){
				my ($first,$second) = split /,/,$brca2{$id};
				if($start >= $first && $end <= $second){
					print OUT $id."\t".$first."\t".$second."\t";
					last;
				}
			}
			$seq_tmp->seq($seq13->subseq($end-24,$end));
			print OUT $seq13->subseq($start,$start+24)."\t".$seq_tmp->revcom()->seq()."\t";
		}
		
		my $length = <FILE>;
		$length =~ s/Amp. Len.//;
		chomp($length);
		$length =~ s/\r//;
		print OUT $length."\t".$start."\t".$end."\n";
	}
}
close OUT;
close FILE;
