# Sandeep Venkataram
# Petrov Lab | Stanford University
# March 18th, 2016
# Pipeline to get barcode counts from Illumina fastq data.
# This pipeline is set up to run on the Stanford Proclus Cluster, you will need to modify software directories to run it on your own system
# The pipeline expects each lane to be processed in its own directory. 
# Each directory should contain a sampleIndexFile.txt file, which is a tab-delimited file with 3 columns and no header
#  The first column is the sample name (e.g. 7-2-1 for batch 7, timepoint 2, replicate 1), the second and third are the primer names used to index the sample (e.g. P104 and P109)


#!/hsgs/projects/petrov/perlroot/perl/bin/perl

#below is the run commands that I used for our data
=pog
#run for expts 1 to 3
qsub -b y -pe orte 12 /hsgs/projects/petrov/perlroot/perl/bin/perl /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/countBarcodesByMapping.pl -si /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts1to3/sampleIndexFile.txt -ft /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/forwardTags.fasta -rt /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/reverseTags.fasta -forwardReference /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/referenceFwd.fasta -reverseReference /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/referenceRev.fasta -r1 /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts1to3/Yeast-01_1.fastq.gz -r2 /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts1to3/Yeast-01_2.fastq.gz -p 11 -odir /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts1to3/

#run for expts 4 to 5
qsub -b y -pe orte 12 /hsgs/projects/petrov/perlroot/perl/bin/perl /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/countBarcodesByMapping.pl -si /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts4to5/sampleIndexFile.txt -ft /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/forwardTags.fasta -rt /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/reverseTags.fasta -forwardReference /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/referenceFwd.fasta -reverseReference /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/referenceRev.fasta -r1 /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts4to5/150416_PINKERTON_0356_BC6CJNACXX_L4_1_pf.fastq.gz -r2 /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts4to5/150416_PINKERTON_0356_BC6CJNACXX_L4_2_pf.fastq.gz -p 11 -odir /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts4to5/

#run for expts 6 to 8
qsub -b y -pe orte 12 /hsgs/projects/petrov/perlroot/perl/bin/perl /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/countBarcodesByMapping.pl -si /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts6to8/sampleIndexFile.txt -ft /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/forwardTags.fasta -rt /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/reverseTags.fasta -forwardReference /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/referenceFwd.fasta -reverseReference /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/referenceRev.fasta -r1 /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts6to8/150423_MARPLE_0302_BC6B58ACXX_L6_1_pf.fastq.gz -r2 /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts6to8/150423_MARPLE_0302_BC6B58ACXX_L6_2_pf.fastq.gz -p 11 -odir /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts6to8/

#run for expts 6 to 8 no spike in
qsub -b y -pe orte 12 /hsgs/projects/petrov/perlroot/perl/bin/perl /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/countBarcodesByMapping.pl -si /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts6to8noSpikeIn/sampleIndexFile.txt -ft /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/forwardTags.fasta -rt /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/reverseTags.fasta -forwardReference /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/referenceFwd.fasta -reverseReference /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/referenceRev.fasta -r1 /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts6to8noSpikeIn/150423_MARPLE_0302_BC6B58ACXX_L7_1_pf.fastq.gz -r2 /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts6to8noSpikeIn/150423_MARPLE_0302_BC6B58ACXX_L7_2_pf.fastq.gz -p 11 -odir /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts6to8noSpikeIn/

#run for expts 13 to 15 with 16 and 17 undigested
qsub -b y -pe orte 12 /hsgs/projects/petrov/perlroot/perl/bin/perl /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/countBarcodesByMapping.pl -si /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts13to15with16and17undigested/sampleIndexFile.txt -ft /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/forwardTags.fasta -rt /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/reverseTags.fasta -forwardReference /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/referenceFwd.fasta -reverseReference /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/referenceRev.fasta -r1 /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts13to15with16and17undigested/150814_BRISCOE_0251_BC7LJDACXX_L2_1_pf.fastq.gz -r2 /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts13to15with16and17undigested/150814_BRISCOE_0251_BC7LJDACXX_L2_2_pf.fastq.gz -p 11 -odir /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts13to15with16and17undigested/

#run for expts 16 to 17 with 13 to 15 undigested
qsub -b y -pe orte 12 /hsgs/projects/petrov/perlroot/perl/bin/perl /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/countBarcodesByMapping.pl -si /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts16to17with13to15undigested/sampleIndexFile.txt -ft /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/forwardTags.fasta -rt /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/reverseTags.fasta -forwardReference /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/referenceFwd.fasta -reverseReference /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/referenceRev.fasta -r1 /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts16to17with13to15undigested/150814_BRISCOE_0251_BC7LJDACXX_L3_1_pf.fastq.gz -r2 /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts16to17with13to15undigested/150814_BRISCOE_0251_BC7LJDACXX_L3_2_pf.fastq.gz -p 11 -odir /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts16to17with13to15undigested/

#run for expts 16 to 19 pool (really just expts 18 and 19)
qsub -b y -pe orte 12 /hsgs/projects/petrov/perlroot/perl/bin/perl /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/countBarcodesByMapping.pl -si /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts16to19pool/sampleIndexFile.txt -ft /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/forwardTags.fasta -rt /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/reverseTags.fasta -forwardReference /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/referenceFwd.fasta -reverseReference /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/referenceRev.fasta -r1 /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts16to19pool/151015_MARPLE_0336_BC7EU8ACXX_L4_1_pf.fastq.gz -r2 /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts16to19pool/151015_MARPLE_0336_BC7EU8ACXX_L4_2_pf.fastq.gz -p 11 -odir /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts16to19pool/

#run for Expts 20-21-23
qsub -b y -pe orte 12 /hsgs/projects/petrov/perlroot/perl/bin/perl /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/countBarcodesByMapping.pl -si /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts20-21-23/sampleIndexFile.txt -ft /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/forwardTags.fasta -rt /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/reverseTags.fasta -forwardReference /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/referenceFwd.fasta -reverseReference /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/referenceRev.fasta -r1 /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts20-21-23/160127_TENNISON_0406_BC8C2DACXX_L7_1_pf.fastq.gz -r2 /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts20-21-23/160127_TENNISON_0406_BC8C2DACXX_L7_2_pf.fastq.gz -p 11 -odir /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts20-21-23/


#run for expts 25 to 27
qsub -b y -pe orte 12 /hsgs/projects/petrov/perlroot/perl/bin/perl /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/countBarcodesByMapping.pl -si /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts25to27/sampleIndexFile.txt -ft /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/forwardTags.fasta -rt /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/reverseTags.fasta -forwardReference /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/referenceFwd.fasta -reverseReference /hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/referenceRev.fasta -r1 /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts25to27/160122_LYNLEY_0532_BC88VCACXX_L4_1_pf.fastq.gz -r2 /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts25to27/160122_LYNLEY_0532_BC88VCACXX_L4_2_pf.fastq.gz -p 11 -odir /hsgs/projects/petrov/sandeep/FitnessAssayData/Expts25to27/

=cut

use strict;
use warnings;

################################### Hardcoded software and output directories #####################################

my $localDir = "/hsgs/projects/petrov/sandeep/FitnessAssayData/bowtiePipeline/";
my $bowtieDir = "/hsgs/projects/petrov/software/bowtie2/";
my $perl = "/hsgs/projects/petrov/perlroot/perl/bin/perl";

###################################################################################################################

use Getopt::Long;

my $dataOutputDirectory;
my $sampleIndexFile;
my $forwardTagsFile;
my $reverseTagsFile;
my $referenceGenomeFwd;
my $referenceGenomeRev;
my $forwardReadsFile;
my $reverseReadsFile;
my $numThreads = 1; #hard coded single threading. Typically takes ~12-24 hours to run on a single lane of hiseq 2000
my %sampleNamingHash = ();
my %outputStringBuffer = ();
my %outputBufferCounter = ();

################################### Select subset of code to run ###############################

my $runBowtieCommands = 1;
my $createBarcodeFiles = 1;
my $generateUniqueCounts = 1;

################################### Preset Variables ###############################

my $bufferLength = 1000;
my $f_bar_start = 64;
my $r_bar_start = 50;
my $barcode_length = 26;

my $fwdTagName = "fwdTags";
my $revTagName = "revTags";
my $referenceNameFwd = "doubleBarcodeReferenceFwd";
my $referenceNameRev = "doubleBarcodeReferenceRev";

################################### Get command line parameters ###############################

GetOptions("sampleIndex|si=s" => \$sampleIndexFile,
	   "forwardTags|ft=s"  => \$forwardTagsFile,
	   "reverseTags|rt=s"  => \$reverseTagsFile,
	   "forwardReference|fr=s"   => \$referenceGenomeFwd,
	   "reverseReference|rr=s"   => \$referenceGenomeRev,
	   "forwardReads|r1=s" => \$forwardReadsFile,
	   "threads|p=s" => \$numThreads,
	   "outputDir|odir=s" => \$dataOutputDirectory,
	   "reverseReads|r2=s" => \$reverseReadsFile);

if(! (-e $sampleIndexFile && -e $forwardTagsFile && -e $dataOutputDirectory && -e $reverseTagsFile && -e $referenceGenomeFwd && -e $referenceGenomeRev && -e $forwardReadsFile && -e $reverseReadsFile)){
	print STDERR "input arguments: \n\n".join("\n",@ARGV)."\nEnd input arguments\n\n";
	print STDERR "Incorrect arguments\n";
	print STDERR "sample index file does not exist\n" if ! -e $sampleIndexFile;
	print STDERR "forward tags file does not exist\n" if ! -e $forwardTagsFile;
	print STDERR "data output directory does not exist\n" if ! -e $dataOutputDirectory;
	print STDERR "reverse tags file does not exist\n" if ! -e $reverseTagsFile;
	print STDERR "ref genome fwd file does not exist\n" if ! -e $referenceGenomeFwd;
	print STDERR "ref genome rev file does not exist\n" if ! -e $referenceGenomeRev;
	print STDERR "forward reads file does not exist\n" if ! -e $forwardReadsFile;
	print STDERR "reverse reads file does not exist\n" if ! -e $reverseReadsFile;
	exit(1);
}


my $fwdTagsSamFile = $dataOutputDirectory."fwdTagsMapping.sam";
my $revTagsSamFile = $dataOutputDirectory."revTagsMapping.sam";
my $doubleBCFwdSamFile = $dataOutputDirectory."doubleBCMappingFwd.sam";
my $doubleBCRevSamFile = $dataOutputDirectory."doubleBCMappingRev.sam";



################################### Run code ###############################


&read_sample_index_file();

if($runBowtieCommands){
	&run_bowtie_commands();
}
if($createBarcodeFiles){
	&extract_barcodes_from_alignments();
}

if($generateUniqueCounts){
	&compute_uniques_from_barcodes();
}




#################################### FUNCTIONS ##############################

sub compute_uniques_from_barcodes{
	foreach my $sampleKey(keys %outputStringBuffer){
		my $samplePrefix = $dataOutputDirectory."sample_$sampleKey";
		system "qsub -b y $perl $localDir"."findUniquesFromBarcodes.pl $samplePrefix";
	}
}

sub read_sample_index_file{
	open(FILE, $sampleIndexFile);
	while(<FILE>){
		chomp;
		next if $_ eq "";
		my @split = split("\t",$_);
		if(! defined $split[1]){
			print "undefined split 1\n".$_."\n";
		}
		my $key = $split[1]."_".$split[2];
		$sampleNamingHash{$key} = $split[0];
		$outputStringBuffer{$key} = "";
		$outputBufferCounter{$key} = 0;
	}
	close(FILE);
}

sub run_bowtie_commands{
#	&system_call($bowtieDir."bowtie2-build $forwardTagsFile $fwdTagName");
#	&system_call($bowtieDir."bowtie2-build $reverseTagsFile $revTagName");




	&system_call($bowtieDir."bowtie2 -p $numThreads --reorder --sensitive-local -x $fwdTagName -U $forwardReadsFile -S $fwdTagsSamFile");
	&system_call($bowtieDir."bowtie2 -p $numThreads --reorder --sensitive-local -x $revTagName -U $reverseReadsFile -S $revTagsSamFile");

	#then map each raw read to barcode sequence itself

	&system_call($bowtieDir."bowtie2-build $referenceGenomeFwd $referenceNameFwd");
	&system_call($bowtieDir."bowtie2-build $referenceGenomeRev $referenceNameRev");
	&system_call($bowtieDir."bowtie2 -p $numThreads --reorder --n-ceil L,0,0.35 --sensitive-local -x $referenceNameFwd -U $forwardReadsFile -S $doubleBCFwdSamFile");
	&system_call($bowtieDir."bowtie2 -p $numThreads --reorder --n-ceil L,0,0.35 --sensitive-local -x $referenceNameRev -U $reverseReadsFile -S $doubleBCRevSamFile");
}


sub extract_barcodes_from_alignments{
	#extract out barcode region (ensure that it is not mapping to restriction site ancestor), and dump barcodes in files split by sample
	foreach my $key(keys %outputBufferCounter){
		system "rm  $dataOutputDirectory"."sample_$key".".barcodes";
	}
	open(FWDTAG, $fwdTagsSamFile);
	open(REVTAG, $revTagsSamFile);
	open(BCMAPFWD, $doubleBCFwdSamFile);
	open(BCMAPREV, $doubleBCRevSamFile);

	my $fwdTagLine = <FWDTAG>;
	my $revTagLine = <REVTAG>;
	my $bcmapFirstLine = <BCMAPFWD>;
	my $bcmapSecondLine = <BCMAPREV>;
	while($fwdTagLine =~ "^@"){
		$fwdTagLine = <FWDTAG>;
	}
	while($revTagLine =~ "^@"){
		$revTagLine = <REVTAG>;
	}
	while($bcmapFirstLine =~ "^@"){
		$bcmapFirstLine = <BCMAPFWD>;
	}
	while($bcmapSecondLine =~ "^@"){
		$bcmapSecondLine = <BCMAPREV>;
	}



	while(1){
		my @fwdTagSplit = split("\t",$fwdTagLine);
		my @revTagSplit = split("\t",$revTagLine);
		my @firstbcmapSplit = split("\t",$bcmapFirstLine);
		my @secondbcmapSplit = split("\t",$bcmapSecondLine);
		
		my $sampleKey = $fwdTagSplit[2]."_".$revTagSplit[2];
		#make sure sample is defined and read maps correctly to dbc locus
		if(defined $sampleNamingHash{$sampleKey} && $firstbcmapSplit[2] !~ "ApaL1" && $secondbcmapSplit[2] !~ "ApaL1" && $firstbcmapSplit[2] ne "*" && $secondbcmapSplit[2] ne "*"){
			my $sampleName = $sampleNamingHash{$sampleKey};
			my ($fseq_id, $fflag, $fstart, $fcigar, $fsequence, $fquality) = @firstbcmapSplit[0,1,3,5,9,10];
			my ($rseq_id, $rflag, $rstart, $rcigar, $rsequence, $rquality) = @secondbcmapSplit[0,1,3,5,9,10];
			my $fwdPCRDupTag = substr($fsequence,0,8);
			my $revPCRDupTag = substr($rsequence,0,8);
			my $fwdbc = get_barcode($fseq_id, $fflag, $fstart, $fcigar, $fsequence, $fquality,$f_bar_start);
			my $revbc = get_barcode($rseq_id, $rflag, $rstart, $rcigar, $rsequence, $rquality,$r_bar_start);
			
			#if($fwdbc ne "" && $revbc ne ""){
				$outputStringBuffer{$sampleKey}.="$fwdbc\t$revbc\t$fwdPCRDupTag\t$revPCRDupTag\n";
				$outputBufferCounter{$sampleKey}++;
				if($outputBufferCounter{$sampleKey} >= $bufferLength){
					write_buffer_to_file($sampleKey,$outputStringBuffer{$sampleKey});
					$outputStringBuffer{$sampleKey} = "";
					$outputBufferCounter{$sampleKey} = 0;
				}
			#}
			
			
		}
		
		
		$fwdTagLine = <FWDTAG>;
		$revTagLine = <REVTAG>;
		$bcmapFirstLine = <BCMAPFWD>;
		$bcmapSecondLine = <BCMAPREV>;
		if(! defined $fwdTagLine || ! defined $revTagLine || ! defined $bcmapFirstLine || ! defined $bcmapSecondLine){
			last;
		}
	}

	foreach my $key(keys %outputStringBuffer){
		if($outputBufferCounter{$key} > 0){
			write_buffer_to_file($key,$outputStringBuffer{$key});
			$outputStringBuffer{$key} = "";
			$outputBufferCounter{$key} = 0;
		}
	}
}

#for each of these files, find consensus sequences of barcodes, make a table of consensus sequences and their frequencies. map consensus to known bc database to get barcode ID
sub write_buffer_to_file{
	my ($sampleKey, $outputString) = @_;
	open(OUT,">>".$dataOutputDirectory."sample_$sampleKey".".barcodes");
	print OUT $outputString;
	close(OUT);
}

sub get_barcode{
	my ($seq_id, $flag, $start, $cigar, $sequence, $quality, $input_bar_start) = @_;
	my $pos = 0;
	my $bar_start = $input_bar_start;
	$bar_start -= $start;
	my $bar_end = $bar_start + $barcode_length;

	if ($cigar eq "*"){
		return "";
	}elsif ($cigar ne "143M"){

	    # the read doesn't align without gaps, so we need to
	    # potentially modify the barcode start and ends based on
	    # the contents of the cigar string, which will enable us
	    # to determine where the gaps in the alignment are

	    # split the cigar string into its constituent parts,
	    # looking for M, I and D operations (and assuming there
	    # are no others - might be useful to check that)

	    my @operations = split(/(M|I|D|S)/, $cigar);

	    for (my $i = 0; $i < @operations; $i+=2){ 

		# our decision as to what to do, based on the observed
		# operation, will change depending on whether we are
		# before or in the barcode

		if ($pos < $bar_start){

		    # we're before the barcode

		    if ($operations[$i+1] eq 'M'){

			# Counts the number of ' M ' and moves current
			# position positively accordingly; barcode
			# position unchanged

			$pos += $operations[$i];

		    }elsif ($operations[$i+1] eq 'I'){

			# Counts number of ' I ' and skips current
			# position positively accordingly, as well as
			# updating the barcode start and end
			# accordingly

			$bar_start += $operations[$i];

			$bar_end += $operations[$i];

			$pos += $operations[$i];

		    }elsif ($operations[$i+1] eq 'D'){

			# Counts number of 'D' and moves current
			# position negatively accordingly

			$bar_start -= $operations[$i];

			$bar_end -= $operations[$i];

			$pos -= $operations[$i];

		    }elsif($operations[$i+1] eq 'S'){
			$pos-=$operations[$i];
			$bar_start+=$operations[$i];
			$bar_end+=$operations[$i];
		    }

		}elsif ($pos >= $bar_start && $pos <= $bar_end){

		    # we have an operation within the barcode

		    if ($operations[$i+1] eq 'I'){
			
			# we have an insertion operation within the
			# barcode, so extend its end by the number of
			# insertions

			$bar_end += $operations[$i];
			$pos     += $operations[$i];

		    }elsif ($operations[$i+1] eq 'D'){
	    
			# we have a deletion operation within the barcode,
			# so decrease its end by the number of deletions

			$bar_end -= $operations[$i];
			$pos     -= $operations[$i];

		    }elsif ($operations[$i+1] eq 'M'){

			$pos   += $operations[$i];

		    }elsif($operations[$i+1] eq 'S'){
				#$pos-=$operations[$i];
				$bar_end+=$operations[$i];
		    }

		}

		# could possible have a last here, if the pos is after
		# the barcode, though its omission shouldn't make any
		# difference
		
	    }

	}

	my $barcode = substr $sequence, $bar_start, $bar_end - $bar_start;
	return $barcode;
	
}

###########################################################################
sub system_call{
###########################################################################
# This subroutine runs any passed in string as a system call, and dies
# if there is an error, otherwise it returns any output from the
# system call as an array of lines.

    my $command = shift;

    my @lines = `$command`;

    if ($?){

	die "An error occured executing the command:\n\n$command\n\n$? : $!\n\n".@lines."\n";

    }

    return @lines;

}
