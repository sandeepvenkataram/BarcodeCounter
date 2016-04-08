# Sandeep Venkataram
# Petrov Lab | Stanford University
# March 18th, 2016
# Pipeline to convert the uniques file generated by countBarcodesByMapping.pl to a tab-delimited barcode count file
# This pipeline is set up to run on the Stanford Proclus Cluster, you will need to modify software directories to run it on your own system
# The consensus_tags_paper.fa file must have a blastn index set up before this script will work properly

#!/hsgs/projects/petrov/perlroot/perl/bin/perl

use strict;
use warnings;


################################### Hardcoded software and input directories #####################################

my $blastn = "/hsgs/projects/petrov/software/ncbi-blast-2.2.28+/bin/blastn";
my $barcodeDb = "/hsgs/projects/petrov/sandeep/EvolutionData/consensus_tags_paper.fa";
my $knownDuplicates = "/hsgs/projects/petrov/sandeep/FitnessAssayData/BarcodeCounter/knownDuplicateBarcodes.txt";
###################################################################################################################


my $dir = shift;


#determine total number of barcodes under consideration from fasta file
my $barcodeDbSize = `wc -l $barcodeDb`;
my @sp = split(" ",$barcodeDbSize);
$barcodeDbSize = $sp[0]/2;


# determine mapping of barcode tags to sample names using sampleIndexFile
my %outputHash = ();
my @outputSamples = ();
my %sampleNameHash = ();
open(SAMP, $dir."sampleIndexFile.txt");
while(<SAMP>){
	chomp;
	my @split = split("\t",$_);
	$sampleNameHash{$split[1]."_".$split[2]} = $split[0];
}


#read in the list of .uniques files generated by countBarcodesByMapping.pl
opendir(Dir, $dir) or die "cannot open directory $dir";
my @docs = grep(/\.uniques$/,readdir(Dir));
close(Dir);
@docs = sort(@docs);
my @sampleNames = ();

my $currentIndent = "";


open(KD, $knownDuplicates);
my @kdLines = <KD>;
close(KD);

my %duplicateMapper = ();

foreach my $line(@kdLines){
	$line = substr($line,1,-2);
	my @split = split("\t",$line);
	#print $line."\n";
	#print join("\t",@split)."\n";
	#exit;
	foreach my $j(@split){
		$duplicateMapper{$j} = $split[0];
	}
}

my $numBarcodesRenamed = 0;
#for every uniques file
foreach my $uniquesFile (@docs) {

	#determine the sample name associated with the file
	my $sn = $uniquesFile;
	$sn=~s/\.uniques//g;
	$sn=~s/sample_//g;
	push @sampleNames, $sampleNameHash{$sn};
	
	$uniquesFile = $dir.$uniquesFile;
	
	#if we have not mapped the uniques sequences in this file to the barcode tag database yet, do so. this step takes a while (many hours)
	if(! -e "$uniquesFile.blastOutput"){
		open(IN,$uniquesFile);
		open(OUT,">".$uniquesFile.".fa");
		my $lineCounter = 1;
		while(<IN>){
			chomp;
			my ($count,$barcode) = split("\t",$_);
			print OUT ">$lineCounter"."_$count\n$barcode\n";
			$lineCounter++;
		}
		close(IN);
		close(OUT);

		system "$blastn -query $uniquesFile.fa -db $barcodeDb -outfmt 6 -word_size 12 -evalue 0.0001 >$uniquesFile.blastOutput";
	}
	
	#read through the blast output generated in the previous step and count how many times each barcode in the master database was present in each sample
	open(FILE, "$uniquesFile.blastOutput");
	my $currentQuery = "";
	my $currentBestBarcode = 0;
	my @currentBestBarcodes = ();
	my $currentBestBarcodeEvalue;
	my %localBarcodeCounts = ();
	foreach my $i(1..$barcodeDbSize){
		$localBarcodeCounts{$i}=0;
	}
	while(<FILE>){
		chomp;
		my @split = split("\t",$_);
		my $query = $split[0];
		my $barcode = $split[1];
		if(defined $duplicateMapper{$barcode}){
			$barcode = $duplicateMapper{$barcode};
			$numBarcodesRenamed++;
		}
		my $evalue = $split[10];
		if($currentQuery eq ""){ #if we are in the initialization step, initialize
			$currentQuery = $query;
			$currentBestBarcode = $barcode;
			@currentBestBarcodes = ($barcode);
			$currentBestBarcodeEvalue = $evalue;
		}
		if($currentQuery ne "" && $currentQuery ne $query){ #if we are looking at the blast matches for a new uniques sequence, store the current count and read data for the new uniques line
			my @split2 = split("_",$currentQuery);
			#$localBarcodeCounts{$currentBestBarcode} += $split2[1];
			my $bc = join("/",sort {$a<=>$b} @currentBestBarcodes);
			$localBarcodeCounts{$bc} = 0 if ! defined $localBarcodeCounts{$bc};
			$localBarcodeCounts{$bc} += $split2[1];
			
			$currentQuery = $query;
			$currentBestBarcode = $barcode;
			@currentBestBarcodes = ($barcode);
			$currentBestBarcodeEvalue = $evalue;
		}else{ #if this is the second or later blast hit for the same uniques sequence, see if it is as significant of a hit as the one we have already. If it is more significant or equally significant, deal with it appropriately
			if($evalue < $currentBestBarcodeEvalue){
				$currentBestBarcode = $barcode;
				$currentBestBarcodeEvalue = $evalue;
				@currentBestBarcodes = ($barcode);
			}
			if($evalue == $currentBestBarcodeEvalue && $currentBestBarcode > $barcode){
				$currentBestBarcode = $barcode;
				push @currentBestBarcodes, $barcode;
			}
		}
	}
	#cleanup to deal with the last line in the file
	my @split2 = split("_",$currentQuery);
	my $bc = join("/", sort {$a<=>$b} @currentBestBarcodes);
	$localBarcodeCounts{$bc} = 0 if ! defined $localBarcodeCounts{$bc};
	$localBarcodeCounts{$bc} += $split2[1];

	#append data to the master hash table
	foreach my $i(keys %localBarcodeCounts){
		if(defined $outputHash{$i}){
			$outputHash{$i}.="\t".$localBarcodeCounts{$i};
		}
		else{
			$outputHash{$i}.="\t".$currentIndent.$localBarcodeCounts{$i}
		}
	}
	foreach my $i(keys %outputHash){
		if(!defined $localBarcodeCounts{$i}){
			$outputHash{$i}.="\t0";
		}
	}
	$currentIndent .= "0\t";
}
print "\nNum Renamed Barcodes $numBarcodesRenamed\n\n";
#write data to file
open(OUT, ">".$dir."barcodeCounts.tab");
print OUT "barcode\t".join("\t",@sampleNames)."\n";
foreach my $i(sort keys %outputHash){
	my $outLine = $outputHash{$i};
	if($outLine =~ "[1-9]"){
		print OUT $i.$outLine."\n";
	}
}
close(OUT);

system "/hsgs/projects/petrov/perlroot/perl/bin/perl condenseCounts.pl $dir"."barcodeCounts.tab >$dir"."condensedBarcodeCounts.tab";
