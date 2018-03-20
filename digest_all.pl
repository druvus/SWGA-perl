#!/usr/bin/env perl -w

use strict;
use warnings;

# Aaron R. Leichty 2013, adapted from parseREBASE.pl

my $seq = '';
my $regexp = '';
my @list = ();
my %rebase;

my $usage = "Usage: digest_all.pl [genome file] [results filename]\nProgram for cataloging number of restriction sites for every enzyme in REBASE. Enzyme database, \"bionet\", can be downloaded from ftp://ftp.neb.com/pub/rebase/bionet.txt.\n";

my $filein = shift or die "$usage\n";
my $fileout = shift or die "$usage\n";

print "\ndigest_all.pl\nSequence being digested: $filein\nFile with results: $fileout\n\n";


open(FILEIN, "< $filein") || die "Unable to open: $!";

local $/ = "\n>";

while (<FILEIN>) { 
	chomp $_;
	$seq = $_;
	$seq = ">$seq";
	$seq =~ s/^>.*\n//; #Remove header line
	$seq =~ s/\s//g; #Remove spaces and newlines

}

print "Done with sequence input.\n";


# Get the REBASE data into a hash, from file "bionet"

open(DATABASE, "< bionet") || die "Unable to open: $!";
 
local $/ = "\n";

while(<DATABASE>){
	
	next if( 1 .. /Rich Roberts/);
	next unless /\S/;
	my(@names) = split /[\s\)\(]+/;
	my $site = IUB_to_regexp(uc($names[-1]));
	if(not exists $rebase{$site} and not exists $rebase{$site}){
		$rebase{$names[0]}=$site;
	}
}

print "Done parsing REBASE.\n";


##count re-sites

foreach (keys %rebase) {    

  	my $query = $_;
        $regexp = $rebase{$query};
	my $counts = 0;
    	while ( $seq =~ /$regexp/ig ) {
		$counts++;
	}
        push (@list, "$query,$counts\n");
} 

open(OUT, ">$fileout") || die "Unable to open: $!";

print OUT @list;

print "Analysis complete.\n";

exit;


############SUBROUTINES###############
 
sub IUB_to_regexp {
 
my($iub) = @_;
 
my $regular_expression = '';
 
my %iub2character_class = (
A => 'A', C => 'C', G => 'G', T => 'T', R => '[GA]', Y => '[CT]', M => '[AC]',
K => '[GT]', S => '[GC]', W => '[AT]', B => '[CGT]', D => '[AGT]', H => '[ACT]',
V => '[ACG]', N => '[ACGT]',
);
$iub =~ s/\^//g;
 
for ( my $i = 0 ; $i < length($iub) ; ++$i ) {
$regular_expression .= $iub2character_class{substr($iub, $i, 1)};
}
return $regular_expression;
}


