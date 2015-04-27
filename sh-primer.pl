#!/usr/bin/perl -w

use strict;

system "perl primer_design.pl";
system "perl primer_filter_by_snp.pl";
system "perl primer_filter_by_exon.pl BRCA1_filter_by_SNP.txt BRCA1_filter_by_exon.txt";
system "perl primer_filter_by_exon.pl BRCA2_filter_by_SNP.txt BRCA2_filter_by_exon.txt";
