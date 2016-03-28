#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(sum max min);

my $infile = shift;

open(FILE, $infile);
my $header = <FILE>;
my @lines = <FILE>;
close(FILE);
my $minReadCutoff = 5;

my %linesTrimmed = ();
my %uniqueBarcodesToKeys = ();
foreach my $line(@lines){
	my @split = split("\t",trim($line));
	my $barcode = shift @split;
	my $foundHigher = 0;
	foreach my $numReads (@split){
		if($numReads > $minReadCutoff){
			$foundHigher = 1;
			last;
		}
	}
	if($foundHigher){
		$linesTrimmed{$barcode} = join("\t",@split);
		my @bcSplit = split("/",$barcode);
		foreach my $bc(@bcSplit){
			if (! defined $uniqueBarcodesToKeys{$bc}){
				$uniqueBarcodesToKeys{$bc} = "$barcode";
			}else{
				$uniqueBarcodesToKeys{$bc}.="\t".$barcode;
			}
		}
	}
}

my %barcodeGroups = ();
my %usedBarcodes = ();
print $header;
foreach my $key(keys %uniqueBarcodesToKeys){
	next if defined $usedBarcodes{$key};
	
	my @barcodeKeysToUse = ($key);
	$usedBarcodes{$key} = 1;
	my @barcodesInGroup = ($key);
	my @barcodeGroupCounts = ();
	while(scalar(@barcodeKeysToUse)>0){
		my $currentBC = shift @barcodeKeysToUse;
		my @groupedBarcodes = split("\t",trim($uniqueBarcodesToKeys{$currentBC}));
		foreach my $origBC(@groupedBarcodes){
			my @bcsp = split("/",$origBC);
			foreach my $bc2(@bcsp){
				next if defined $usedBarcodes{$bc2};
				push @barcodeKeysToUse, $bc2;
				push @barcodesInGroup, $bc2;
				$usedBarcodes{$bc2} = 1;
			}
			
			#add to array
			if(scalar(@barcodeGroupCounts) ==0){
				@barcodeGroupCounts = split("\t",$linesTrimmed{$origBC});
			}else{
				my @barcodeGroupCounts2 = split("\t",$linesTrimmed{$origBC});
				foreach my $i(0..scalar(@barcodeGroupCounts)-1){
					$barcodeGroupCounts[$i]+=$barcodeGroupCounts2[$i];
				}
			}
		}
	}
	print join("/",sort {$a<=>$b} @barcodesInGroup)."\t".join("\t",@barcodeGroupCounts)."\n";
}


sub trim
{
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}