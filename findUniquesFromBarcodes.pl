#!/hsgs/projects/petrov/perlroot/perl/bin/perl

use strict;
use warnings;

my $barcodeFilePrefix = shift;

my %uniqueBarcodePCRDupTracker = ();
my %uniqueBarcodeCounter = ();
open(IN,$barcodeFilePrefix.".barcodes");
while(<IN>){
	chomp;
	my ($fwd,$rev,$fwdpcr,$revpcr) = split("\t",$_);
	my $key = $fwd."_".$rev;
	my $pcrKey = $fwdpcr."_".$revpcr."\t";
	if(! defined $uniqueBarcodePCRDupTracker{$key}){
		$uniqueBarcodeCounter{$key} = 1;
		$uniqueBarcodePCRDupTracker{$key} = "\t".$pcrKey;
	}else{
		next if $uniqueBarcodePCRDupTracker{$key} =~ "\t$pcrKey";
		$uniqueBarcodeCounter{$key} ++;
		$uniqueBarcodePCRDupTracker{$key} .= $pcrKey;
	}
}
close(IN);
open(OUT,">".$barcodeFilePrefix.".uniques");
foreach my $key(keys %uniqueBarcodeCounter){
	print OUT $uniqueBarcodeCounter{$key}."\t".$key."\n";
}
close(OUT);
