#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

sub Info
	{
	print "
#######################################################################
# Script for automatic Quality Control processing and merging for
# 	Illumina paired end sequencing files in FastQ format.
# Dependencies:
# 	-cutadapt: <link to cutadapt>
# 	-TrimGalore: <link to TrimGalore>
# 	-FLASH: <link to FLASH>
#
# Needs paired end Illumina files (one file for each read pair)
# Author: piyuranjan\@gmail.com
# More info / Wiki: Please read <link to wiki page> for a full
# 	documentation with example usage
#######################################################################
	\n";
	}
sub Usage
	{
	print "
Usage:\n$0 [options]
\nOptions:
 -f|forward	[string:required] Comma separated list of forward read (R1) file/s in FastQ.
 -r|reverse	[string:required] Comma separated list of reverse read (R2) file/s in FastQ.
 -i|intrlv	[string:required-conditional] Comma separated list of interleaved file/s in FastQ.
					Default: paired-end input takes precedence, interleaved file/s ignored
						if -f, -r available.
 -l|label	[string:required] Comma separated list of unique label/s for FastQ datasets.
 -s|stats	[string:optional] Stats file to export QC, merging statistics.
					Default: No stats calculated
 -w|workDir	[string:optional] Path location for output files/folders
					Default: ./ (present working directory)
 -fa|fAdapt	[string:optional] Custom adapter sequence expected at 3' of R1.
					Default: Illumina nextera CCGAGCCCACGAGAC
 -ra|rAdapt	[string:optional] Custom adapter sequence expected at 3' of R2.
					Default: Illumina nextera GACGCTGCCGACGA
 -ml|minLen	[integer:optional] Min length of reads to keep.
					Default: 50
 -t|threads	[integer:optional] Max threads dependency programs can use.
					Default: 1
 -ri|rmIlv	Remove paired-end files after finishing. Works only with -i. Saves disc space.
 -v|verbose	print logging information.
 -h|help	show help (information and usage); exit 0.
	\n";
	}
sub MoreInfo
	{print "\nPlease check out <link to wiki page> for detailed information.\n";}

## All input options
my ($forFiles,$revFiles,$intlFiles,$label,$statFile,$wDir,$fAdapt,$rAdapt,$minLen,$threads,$rmIntl,$verbose,$help);
$wDir='./';
$fAdapt='CCGAGCCCACGAGAC';
$rAdapt='GACGCTGCCGACGA';
$minLen=50;
$threads=1;
if(!GetOptions('f|forward=s' => \$forFiles,
				'r|reverse=s' => \$revFiles,
				'i|intrlv=s' =>\$intlFiles,
				'l|label=s' => \$label,
				's|stats=s' => \$statFile,
				'w|workDir=s' => \$wDir,
				'fa|fAdapt=s' => \$fAdapt,
				'ra|rAdapt=s' => \$rAdapt,
				'ml|minLen=i' => \$minLen,
				't|threads=i' => \$threads,
				'ri|rmIlv' => \$rmIntl,
				'v|verbose' => \$verbose,
				'h|help' => \$help))
	{Usage;exit 1;} #quit with error code
if($help) #quit with help
	{Info;Usage;MoreInfo;exit 0;}
# if((!defined $forFiles)||(!defined $revFiles)||(!defined $label)) #these options are necessary
#	{Usage;MoreInfo;exit 1;}
unless ((defined $label)&&((defined $intlFiles)||((defined $forFiles)&&(defined $revFiles)))) #these options are necessary
	{Usage;MoreInfo;exit 1;}

## Check for package dependencies
# specify everything as a separate function
CheckCutadapt(); CheckTrimGalore(); CheckFlash();

## Store forward and reverse fastQ files
my @labels=split(/,/,$label);
my (@forFqs,@revFqs,@intlFqs);
if((defined $forFiles)&&(defined $revFiles)&&(defined $intlFiles))
	{MoreInfo;die "\nError: Can't process paired and interleaved reads simultaneously.\n";}
elsif((defined $forFiles)&&(defined $revFiles)) #in case of paired input
	{
	@forFqs=split(/,/,$forFiles);
	@revFqs=split(/,/,$revFiles);
	}
elsif(defined $intlFiles) #in case of interleaved input
	{
	@intlFqs=split(/,/,$intlFiles);
	foreach my $intl(@intlFqs)
		{
		my $tempLabel=$intl;
		$tempLabel=~s/\.f(ast)?q$//;
		push(@forFqs,$tempLabel."-R1.fq");
		push(@revFqs,$tempLabel."-R2.fq");
		}
	}
else #die otherwise
	{MoreInfo;die "\nFatal Error: Input FastQ is neither paired nor interleaved.\n";}


## Create directory structure
$wDir.='/' if($wDir !~ /\/$/); #place trailing slash
my $qcDir=$wDir.'01.qced';
my $mergeDir=$wDir.'02.merged';
my $qualDir=$wDir.'03.qualified';
if(! -d $wDir)
	{
	print "\n## Working directory $wDir does not exist..." if $verbose;
	mkdir $wDir,0755;
	print "Created successfully\n" if $verbose;
	}
else
	{print "\n## Working directory $wDir already exists... datasets will be added/overwritten\n" if $verbose;}
mkdir $qcDir,0755 if(! -d $qcDir);
mkdir $mergeDir,0755 if(! -d $mergeDir);
mkdir $qualDir,0755 if(! -d $qualDir);
CheckLeftOver($wDir); #check for any files that look similar to result files in work dir

## Open a stats file if option given
my $STAT;
if(defined $statFile)
	{
	open($STAT,">",$wDir.$statFile) or die $!; #open a file to record qc/merge stats
	print $STAT "Dataset\tRawReads\tQCPaired\tQCSingle1\tQCSingle2\tPairsCombined\tPairsNotCombined\tTotalSingles\tTotalPaired";
	}

## Process each PE FastQ pair
foreach my $i(0..$#labels)
	{
	print "\n## Processing $labels[$i] ##";
	
	## Split interleaved input to pairs if provided
	if(defined $intlFiles) #if interleaved files are provided as input
		{
		print "\n# Splitting interleaved to paired-end for $labels[$i]..." if $verbose;
		unless(-r $intlFqs[$i]){print "\n#Cannot find/read $intlFqs[$i]... Skipping";next;}
		open(my $I,$intlFqs[$i]) or die $!;open(my $F,">",$forFqs[$i]) or die $!;open(my $R,">",$revFqs[$i]) or die $!;
		my ($ifqRec,$ffqLine,$rfqLine)=(0)x3; #FQ file sanity checks
		while(<$I>)
			{
			$ifqRec++;
			my ($fline,$rline);
			$fline=$_; $ffqLine++;
			foreach(0..2){$fline.=<$I>; $ffqLine++;}
			foreach(0..3){$rline.=<$I>; $rfqLine++;}
			print $F $fline; print $R $rline;
			}
		die "\nFatal Error: Inconsistency found during interleaved FastQ splitting for $labels[$i]\nPaired-end records in interleaved file: $ifqRec\nR1 lines: $ffqLine\tR2 lines: $rfqLine\n" if(($ifqRec!=$ffqLine/4)||($ifqRec!=$rfqLine/4)); #sanity check for interleave split
		print " successful" if $verbose;
		}
	##delete paired end files if initial input was interleaved, maybe give a flag for that
	
	## Confirm input files are valid/readable files
	unless(-r $forFqs[$i]){print "\n#Cannot find/read $forFqs[$i]... Skipping";next;}
	unless(-r $revFqs[$i]){print "\n#Cannot find/read $revFqs[$i]... Skipping";next;}
	
	## Adapter and qual trim PE FastQ files
	print "\n# Performing QC/trimming for $labels[$i]..." if $verbose;
	my $trimLen=$minLen-10;
	my $trimGaloreLog=$labels[$i]."_trimGalore.log";
	my $trimGalore="trim_galore -q 25 --length $trimLen --paired --retain_unpaired -a $fAdapt -a2 $rAdapt -r1 $minLen -r2 $minLen $forFqs[$i] $revFqs[$i] >$trimGaloreLog 2>&1";
	my $trimGaloreOut=`$trimGalore`;
	die "\nCould not finish trimGalore QC for $labels[$i]\nPlease check $trimGaloreLog" if($?); #trimGalore sanity check
	`mv $trimGaloreLog $forFqs[$i]_trimming_report.txt $revFqs[$i]_trimming_report.txt $qcDir`; #move all trimming reports
	my $trimGaloreVal1=$forFqs[$i]; $trimGaloreVal1=~s/\.f(ast)?q$//;$trimGaloreVal1.="_val_1.fq";
	my $trimGaloreVal2=$revFqs[$i]; $trimGaloreVal2=~s/\.f(ast)?q$//;$trimGaloreVal2.="_val_2.fq";
	my $trimGaloreUnpaired1=$forFqs[$i]; $trimGaloreUnpaired1=~s/\.f(ast)?q$//;$trimGaloreUnpaired1.="_unpaired_1.fq";
	my $trimGaloreUnpaired2=$revFqs[$i]; $trimGaloreUnpaired2=~s/\.f(ast)?q$//;$trimGaloreUnpaired2.="_unpaired_2.fq";
	print " successful" if $verbose;
	
	## Merging PE FastQ files
	print "\n# Performing merging for $labels[$i]..." if $verbose;
	my $flashLog=$labels[$i]."_flash.log";
	my $flash="flash -m 25 -x 0.1 -t $threads -o $labels[$i] $trimGaloreVal1 $trimGaloreVal2 >$flashLog";
	my $flashOut=`$flash`;
	die "\nCould not finish FLASH merging for $labels[$i]\nPlease check $flashLog" if($?); #flash sanity check
	`mv $flashLog $labels[$i].hist $labels[$i].histogram $mergeDir`; #move all merging reports
	my $flashExtended=$labels[$i].'.extendedFrags.fastq';
	my $flashNCombined1=$labels[$i].'.notCombined_1.fastq';
	my $flashNCombined2=$labels[$i].'.notCombined_2.fastq';
	print " successful" if $verbose;
	
	## Compiling qualified reads with sequencing pair information preserved
	print "\n# Compiling reads for $labels[$i]..." if $verbose;
	my $currTime=time(); #use current time as an end label to avoid accidental overwrite of raw reads, in case the filename is same as the one produced by the script.
	#compile single end reads
	my $qualifiedSeFinal=$labels[$i]."-SE.fq";
	my $qualifiedSe=$labels[$i]."-SE.fq".$currTime;
	`touch $qualifiedSe`;
	foreach my $file($flashExtended,$trimGaloreUnpaired1,$trimGaloreUnpaired2)
		{
		#my $se="sed 's/\\s/-/' $file >>$qualifiedSe"; #this was to incorporate strand info in read header, useful for mapping but could be harmful for assembly
		my $se="cat $file >>$qualifiedSe";
		`$se`;
		}
	#compile paired end reads
	my $qualifiedR1Final=$labels[$i]."-R1.fq";
	my $qualifiedR2Final=$labels[$i]."-R2.fq";
	my $qualifiedR1=$labels[$i]."-R1.fq".$currTime;
	my $qualifiedR2=$labels[$i]."-R2.fq".$currTime;
	#my $pairedR1="sed 's/\\s/-/' $flashNCombined1 >$qualifiedR1"; #this was to incorporate strand info in read header, useful for mapping but could be harmful for assembly
	my $pairedR1="cp $flashNCombined1 $qualifiedR1";
	`$pairedR1`;
	#my $pairedR2="sed 's/\\s/-/' $flashNCombined2 >$qualifiedR2"; #this was to incorporate strand info in read header, useful for mapping but could be harmful for assembly
	my $pairedR2="cp $flashNCombined2 $qualifiedR2";
	`$pairedR2`;
	print " successful" if $verbose;
	
	## Calculate procedure statistics
	if(defined $statFile)
		{
		print "\n# Calculating qc/merging stats for $labels[$i]..." if $verbose;
		#Dataset RawReads QCPaired QCSingle1 QCSingle2 PairsCombined PairsNotCombined TotalSingles TotalPaired
		print $STAT "\n$labels[$i]";
		my @readStats;
		push(@readStats,$forFqs[$i]); #RawReads
		push(@readStats,$trimGaloreVal1,$trimGaloreUnpaired1,$trimGaloreUnpaired2); #QCPaired QCSingle1 QCSingle2
		push(@readStats,$flashExtended,$flashNCombined1); #PairsCombined PairsNotCombined
		push(@readStats,$qualifiedSe,$qualifiedR1); #TotalSingles TotalPaired
		foreach my $j(@readStats)
			{
			my $val=`wc -l $j`; #calculate number of lines
			my $value=(split(/\s+/,$val))[0]; #get the number
			chomp($value);
			$value/=4; #since this is a FastQ file
			print $STAT "\t$value";
			}
		print " successful" if $verbose;
		}
	
	## Move all files to respective folders
	`mv $trimGaloreVal1 $trimGaloreVal2 $trimGaloreUnpaired1 $trimGaloreUnpaired2 $qcDir`;
	`mv $flashExtended $flashNCombined1 $flashNCombined2 $mergeDir`;
	`mv $qualifiedSe $qualifiedR1 $qualifiedR2 $qualDir`;
	#rename final files to remove the timestamp inside the qualifiedReads DIR
	`mv $qualDir/$qualifiedSe $qualDir/$qualifiedSeFinal`; `mv $qualDir/$qualifiedR1 $qualDir/$qualifiedR1Final`; `mv $qualDir/$qualifiedR2 $qualDir/$qualifiedR2Final`;
	
	## Remove paired-end files if interleaved input and option provided
	`rm $forFqs[$i] $revFqs[$i]` if((defined $intlFiles)&&(defined $rmIntl));
	
	print "\n# $labels[$i] finished #\n";
	#remove any leftover files for example test-R1_trimmed.fq created when trimGalore didn't work
	}
print "\n";
if(defined $statFile){print $STAT "\n"; close($STAT);}


#########################
### Subroutines below ###
#########################

sub CheckCutadapt #check installation of the package Cutadapt, quit if not found
	{
	print "\n## Checking Cutadapt... " if $verbose;
	`cutadapt -h 1>/dev/null 2>&1`;
	die "Fatal error: cutadapt not installed (or not available in PATH)\n$!".MoreInfo if($?);
	print "cutadapt ok\n" if $verbose;
	}
sub CheckTrimGalore #check installation of the package TrimGalore, quit if not found
	{
	print "\n## Checking TrimGalore... " if $verbose;
	`trim_galore -h 1>/dev/null 2>&1`;
	die "Fatal error: trim_galore not installed (or not available in PATH)\n$!".MoreInfo if($?);
	print "trim_galore ok\n" if $verbose;
	}
sub CheckFlash #check installation of the package FLASH, quit if not found
	{
	print "\n## Checking FLASH... " if $verbose;
	`flash -h 1>/dev/null 2>&1`;
	die "Fatal error: flash not installed (or not available in PATH)\n$!".MoreInfo if($?);
	print "flash ok\n" if $verbose;
	}
sub CheckLeftOver #Check if an file is left over from a separate analysis/execution
	{
	my ($wDir)=@_;
	my $leftoverFlag=0;
	opendir(my $DIR,$wDir) or die $!;
	while(my $file=readdir($DIR))
		{
		die "\nError: Detecting leftover files from a previous analysis in $wDir\nPlease remove $file or provide a fresh location\n" if($file=~/(_trimming_report\.txt|_val_.\.fq|_unpaired_.\.fq|extendedFrags\.fastq|notCombined_.\.fastq|\.hist|\.histogram)/);
		}
	}

