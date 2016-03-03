#!/usr/bin/perl -w

use warnings;
use strict;

#Aaron Reid Leichty, 2014

my %dg_motifs;
my %cg_motifs;
my %sub_dg_motifs;
my @motif_results;

my $usage = "Usage: SWGA.pl [filename of desired genome] [filename of contaminating genome] [size range, e.g. 8-12 or 8] [Tm limit] [number of best desired-motifs] [results filename]\nProgram for finding motifs for SWGA. Genome files should be in fasta format.\n";

my $dg_file = shift or die "$usage\n";
my $cg_file = shift or die "$usage\n";
my $size_range = shift or die "$usage\n";
my $tm_limit = shift or die "$usage\n";
my $tophitnumber = shift or die "$usage\n";
my $outp_file = shift or die "$usage\n";

print "\nSWGA.pl\nDesired Genome: $dg_file\nContaminating Genome: $cg_file\nMotif Size: $size_range\nTm limit: $tm_limit\nTop motifs: $tophitnumber\nOutput file: $outp_file\n\n";

my $motif_low;
my $motif_high;

if($size_range =~ m/(\d+)-(\d+)/){
	$motif_low = $1;
	$motif_high = $2;
}
elsif($size_range =~ m/(\d+)/){
	$motif_low = $1;
	$motif_high = $1;	
}

open(DESIRED, "< $dg_file") || die "Unable to open: $!";
open (CONTAMINATING, "< $cg_file")  || die "Unable to open: $!";
open(RESULTS,">$outp_file") || die "Unable to create: $!";
print RESULTS "Genomewide analysis for top $tophitnumber most frequent motifs in $dg_file with $cg_file as contaminating genome.\n";
print RESULTS "motif\tnumber of desired sites\tnumber of contaminating sites\tRatio CG:DG\tTm\n";


##Cataloging motifs in desired genome

print "Evaluating desired genome.\n";

local $/ = "\n>";

foreach (<DESIRED>) { #iterate through all chromosomes
	chomp $_;
	my $seq = $_;
	$seq = ">$seq";
	$seq =~ m/^>(.*)\n/; #Chromosome ID
	my $id = $1;
	$seq =~ s/^>.*\n//; #Remove header line
	$seq =~ s/\s//g; #Remove spaces and newlines
	my $length = length $seq;
	print "$id\n";
	my $i;
	for ($i = $motif_low; $i <= $motif_high; $i++ ){ 
		my $j;
		for ( $j = 1; $j < $length-($i-1); $j++){ 
			my $window = substr($seq,$j,$i);
			my $a=($window =~tr/A//);			
			my $b=($window =~tr/C//);
			my $c=($window =~tr/G//);
			my $d=($window =~tr/T//);
			my $motif_temp = $a*2+$d*2+$b*4+$c*4;
			if ($motif_temp > $tm_limit){next;}  
			$dg_motifs{$window}++;
		}		
	}
}

close DESIRED;

my @dg_sorted_motifs = sort SortDescending(keys(%dg_motifs));


#Summing motif and its reverse complement, then creating a new sorted list

foreach (@dg_sorted_motifs) {
	my $motif = $_;
	my $revcom_motif = reverse $motif;
	$revcom_motif =~ tr/ACGTacgt/TGCAtgca/;
	my $motif_count = $dg_motifs{$motif};
	my $revcom_count;
	if (exists ($dg_motifs{$revcom_motif})){
		$revcom_count = $dg_motifs{$revcom_motif};
	}
	else{
		$revcom_count = 0;
	}
	my $total_count = $motif_count + $revcom_count;
	if (exists ($sub_dg_motifs{$revcom_motif})){
		next;
	}
	else{
		$sub_dg_motifs{$motif} = $total_count;
	}
}

my @sub_dg_sorted_motifs = sort SortDescending(keys(%sub_dg_motifs));

print "Finished cataloging desired genome.\n\n";

print "Evaluating contaminating genome.\n";


#Cataloging motifs in contaminating genome

local $/ = "\n>";

foreach (<CONTAMINATING>) { #iterate through all chromosomes
	chomp $_;
	my $seq = $_;
	$seq = ">$seq";
	$seq =~ m/^>(.*)\n/; #Chromosome ID
	my $id = $1;
	$seq =~ s/^>.*\n//; #Remove header line
	$seq =~ s/\s//g; #Remove spaces and newlines
	print "$id\n";
	my $k=0;
	foreach (@sub_dg_sorted_motifs) {
		my $dg_motifs_key = $_;
		my $revcom_key = reverse $dg_motifs_key;
		$revcom_key =~ tr/ACGTacgt/TGCAtgca/;
		while ($seq =~ m/$dg_motifs_key|$revcom_key/g){
			$cg_motifs{$dg_motifs_key}++;	
		}
		$k++;
		last if $k == $tophitnumber;
	}		
}

close CONTAMINATING;

print "Finished with contaminating genome.\n";


#Writing report

my $k = 0;

foreach (@sub_dg_sorted_motifs) {
	my $dg_motifs_key = $_;
	my $ratio;
	my $tm = temp($dg_motifs_key);	#calculating melting temp
	if(exists $cg_motifs{$dg_motifs_key}){
		$ratio = $sub_dg_motifs{$dg_motifs_key}/$cg_motifs{$dg_motifs_key};
		push (@motif_results, "$dg_motifs_key\t" . "$sub_dg_motifs{$dg_motifs_key}\t" . "$cg_motifs{$dg_motifs_key}\t" . "$ratio\t" . "$tm\n");
	}
	else{
		$ratio = $sub_dg_motifs{$dg_motifs_key}/0.1;
		push (@motif_results, "$dg_motifs_key\t" . "$sub_dg_motifs{$dg_motifs_key}\t" . "0\t" . "$ratio\t" . "$tm\n");
	}
	$k++;
	last if $k == $tophitnumber;
}

print RESULTS @motif_results;

print "\nAnalysis complete.\n";

exit;


####SUBROUTINES####

#SORTASCENDING

sub SortAscending {
   $dg_motifs{$a} <=> $dg_motifs{$b};
}

#SORTDESCENDING

sub SortDescending {
   $dg_motifs{$b} <=> $dg_motifs{$a};
}


#TM

sub temp {
	my $motif = $_;
	my $tm="";
	my $a=($motif =~tr/A//);			
	my $b=($motif =~tr/C//);
	my $c=($motif =~tr/G//);
	my $d=($motif =~tr/T//);
	$tm = $a*2+$d*2+$b*4+$c*4;
	return $tm;
}
