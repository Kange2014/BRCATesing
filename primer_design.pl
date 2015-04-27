#!/usr/bin/perl -w

use strict;

#die "Usage: perl primer_design.pl <sequence_fasta_file> <output_file_name>\n" unless scalar @ARGV >= 2;

# design some primers.
# the output will be put into temp.out
use lib "/home/wangyj/BioPerl-1.6.901/";
use Bio::Perl;
use Bio::Seq;
use Bio::SeqIO;
use Primer::Design;
use Bio::Location::SplitLocationI;

my @chr13 = read_all_sequences("BRCA2.gb",'Genbank');
my @chr17 = read_all_sequences("BRCA1.gb",'Genbank');

my $seq_obj = Bio::Seq->new(-id => '0001', -seq => 'AAAA');
my $brca1_exon = Bio::SeqIO->new(-file=> ">BRCA1_exons.fa",-format=>'fasta');
my $brca2_exon = Bio::SeqIO->new(-file=> ">BRCA2_exons.fa",-format=>'fasta');

my $seq13 = $chr13[0]; 
my $seq17 = $chr17[0];

my $junction_size = 200;
my $overlap = 500;
my $window_size = 600;

my @features = $seq17->all_SeqFeatures();
foreach my $feat (@features){
	if($feat->primary_tag eq 'CDS' ){
		my @value=$feat->each_tag_value('gene');
		if($value[0] eq "BRCA1"){
			if ($feat->location->isa('Bio::Location::SplitLocationI')){
				my $count = 2;
				foreach my $loc ($feat->location->sub_Location){
					if($loc->end - $loc->start >= $window_size){
						my $num = int( ($loc->end - $loc->start)/$window_size ) + 1;
						for(my $i = 1; $i <= $num; $i++){
							$overlap = $junction_size if $i == 1;
							my $start = ($i - 1)*$window_size + $loc->start - $overlap;
							my $end = $loc->start + $i*$window_size - 1;
							$end = $loc->end + $junction_size if $i == $num;
							$seq_obj->seq($seq17->subseq($start, $end));
							$start = $loc->start if $i == 1;
							$end = $loc->end if $i == $num;
							$seq_obj->display_id($count."-$i:".$start.",".$end);
							$brca1_exon->write_seq($seq_obj);
						}
					}
					else{
						$seq_obj->display_id($count.":".$loc->start.",".$loc->end);
						$seq_obj->seq($seq17->subseq($loc->start - $junction_size, $loc->end + $junction_size));
						$brca1_exon->write_seq($seq_obj);
					}
					$count++;
					$count++ if $count == 4;
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
				my $count = 2;
				foreach my $loc ($feat->location->sub_Location){
					if($loc->end - $loc->start >= $window_size){
						my $num = int( ($loc->end - $loc->start)/$window_size ) + 1;
						for(my $i = 1; $i <= $num; $i++){
							$overlap = $junction_size if $i == 1;
							my $start = ($i - 1)*$window_size + $loc->start - $overlap;
							my $end = $loc->start + $i*$window_size - 1;
							$end = $loc->end + $junction_size if $i == $num;
							$seq_obj->seq($seq13->subseq($start, $end));
							$start = $loc->start if $i == 1;
							$end = $loc->end if $i == $num;
							$seq_obj->display_id($count."-$i:".$start.",".$end);
							$brca2_exon->write_seq($seq_obj);
						}
					}
					else{
						$seq_obj->display_id($count.":".$loc->start.",".$loc->end);
						$seq_obj->seq($seq13->subseq($loc->start - $junction_size, $loc->end + $junction_size));
						$brca2_exon->write_seq($seq_obj);
					}
					$count++;
				}
			}
			last;
		}
	}
}

design("BRCA1_exons.fa",$junction_size,$overlap,$window_size,"BRCA1_primers.csv");
design("BRCA2_exons.fa",$junction_size,$overlap,$window_size,"BRCA2_primers.csv");

sub design{
my ($seqfile,$junction_size,$overlap,$window_size,$primer_file) = @_;

my @seq_array_ref = read_all_sequences($seqfile,"fasta");

open OUT,">$primer_file";
print OUT "EXON_ID,OLIGO,start,len,tm,gc%,any,3',seq\n";
foreach my $seq_ref(@seq_array_ref){
	next if($seq_ref->length() <= 100);
	my $primerobj; 
	eval{ $primerobj = Design->new(-seq => $seq_ref); };
	if($@){ print $seq_ref->display_id(),"\n";}
	
###################################################################################
# what are the arguments, and what do they mean?
# you can add following codes into this script to check them (just by deleting the "#"
# in the first of each line).
# 
 #my $args = $primerobj->arguments;
 #print "ARGUMENT\tMEANING\n";
 #foreach my $key (keys %{$args}) {print "$key\t", $$args{$key}, "\n"}
#
###################################################################################

# if you hope to customize primer parameters, you can add any legal value to the package. 
# Please use $primer3->arguments to find a list of all the values that are allowed,
# or see the primer3 docs.
# 
# for example,
#
# set the number of designed primers
 $primerobj->add_targets(PRIMER_NUM_RETURN=>"50"); 

# set the region that primer product must cover
# If one or more targets is specified then a legal primer pair must flank at least one of them
# TRAGET: (interval list, default empty) Regions that must be included in the product. 
# The value should be a space-separated list of <start>,<length>

my $intronic_size = 80; ### reporting variants to +50 or -50 into the intronic regions

 if($seq_ref->display_id() =~ /\-[2-9]/ || $seq_ref->display_id() =~ /\-[1-9]{2}/){ 
	my $start = $overlap/2 + 1;
	if($seq_ref->length <= $window_size + $overlap){  ### sliding window size + overlap 
		my $effect_size = $seq_ref->length() - $overlap/2 - $junction_size + $intronic_size;	
		$primerobj->add_targets(TARGET => "$start,$effect_size");
	}
	else{
		### effect size: sliding window size - overlap + 30*2(30 upstream and downstream of target product)
		my $effect_size = $window_size - $overlap/2 + $intronic_size*2;
		$primerobj->add_targets(TARGET => "$start,$effect_size"); 
	}
 }
 else{ 
	### 
	my $effect_size;
	$effect_size = $seq_ref->length() - $junction_size*2 + $intronic_size*2 unless $seq_ref->display_id() =~ /\-1/;
	$effect_size = $seq_ref->length() - $junction_size - $overlap + $intronic_size if $seq_ref->display_id() =~ /\-1/;
	### start: 150 + 1 - 5
	my $start = $junction_size + 1 - $intronic_size;
	$primerobj->add_targets(TARGET => "$start,$effect_size"); 
 }
 
# set a sub-region of the given sequence to pick primers
# For example, often the first dozen or so bases of a sequence are vector, and should be excluded 
# from consideration. The value for this parameter has the form <start>,<length> where <start> is 
# the index of the first base to consider, and <length> is the number of subsequent bases in 
# the primer-picking region.
# $primerobj->add_targets(INCLUDED_REGION => "523988,900"); 
 
# set the maximum and minimum Tm of the primer
# $primerobj->add_targets(PRIMER_OPT_TM=>"55");
# $primerobj->add_targets(PRIMER_MIN_TM=>"50");
# $primerobj->add_targets(PRIMER_MAX_TM=>"60");
	
# set the range of primer product size
 $primerobj->add_targets(PRIMER_PRODUCT_SIZE_RANGE => "200-800");
	
# set the maximum and minimum size of the primer
# $primerobj->add_targets(PRIMER_OPT_SIZE=>"21");
# $primerobj->add_targets(PRIMER_MIN_SIZE=>"18");
# $primerobj->add_targets(PRIMER_MAX_SIZE=>"25");
	
#####################################################################################
# you can also change the program_name
# 
# $primerobj->program_name('my_suprefast_primer3');
 unless ($primerobj->executable) {
 	print STDERR "primer3 can not be found. Is it installed?\n";
 	exit(-1)
 }
#
# or change the primer3's path (default: /usr/bin/primer3_core): 
#
# $primerobj = Design->new(-seq => $seq_ref, -path => /home/usrname/primer3/primer3_core);
######################################################################################

	my $primer_info = $primerobj->run();
	print OUT $primer_info;
}
close OUT;

}