#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

sub Info
	{
	print "
#######################################################################
# Script for calculating abundance of a genome/bin and it's sequencing
#	depth (per NT) statistics in a meta-omic dataset using BAM files.
# Dependencies:
#	-samtools >1.3: http:\//www.htslib.org
#	-bedtools: http:\//bedtools.readthedocs.io/en/latest
#	-R: https:\//www.r-project.org/
#
# Needs a BAM alignment file for sequencing reads mapped to a
#	genome/bin as a reference. Only reads that align are recommended.
#
# Author: piyuranjan\@gmail.com
# More info / Wiki: Please read <link to wiki page> for a full
# 	documentation with example usage
#######################################################################
	\n";
	}
sub Usage
	{
	print "
Usage:\n$0 [options] >[abundDepthStatsFile]
\nOptions:
 -b|bam		[string:required] BAM alignment file with reads that align to
				a reference genome/bin.
 -l|label	[string:optional] A label for the output record.
				Default: BAM fileName (without .bam)
 -s|scale	[integer:optional] Total reads in the omic dataset to normalize
				abundance. If given, abundance will be reported as a
				percentage of total reads.
				Default: off; abundance is reported as # reads mapped.
 -n|noHeader	Silence column headers in output.
 -header	Print only header; exit 0.
 -debug		No need to use this unless you know what it does.
 -h|help	Show help (information and usage); exit 0.

Output order:
#MappedOmicDataset\t#Abundance\t#MinDep\t#1stQu\t#MedDep\t#3rdQu\t#MaxDep\t#MeanDep\t#SDDep\n
	\n";
	}

## All input options
my ($bamFile,$label,$scale,$noHeader,$header,$debug,$help);
#$scale=0;
my $outHeader="#MappedOmicDataset\t#Abundance\t#MinDep\t#1stQu\t#MedDep\t#3rdQu\t#MaxDep\t#MeanDep\t#SDDep\n";
my $output="";
if(!GetOptions(
	'b|bam=s' => \$bamFile,
	'l|label=s' => \$label,
	's|scale=i' => \$scale,
	'n|noHeader' =>\$noHeader,
	'header' =>\$header,
	'debug' =>\$debug,
	'h|help' => \$help))
	{Usage;exit 1;} #quit with error code
if($help) #quit with help
	{Info;Usage;exit 0;}
if($header) #quit with output header if only header requested
	{print "$outHeader\n";exit 0;}
unless (defined $bamFile) #these options are necessary
	{Usage;exit 1;}

## Assign label
if(!defined $label)
	{
	$label=$bamFile;
	$label=~s/.*\///; #remove path to file
	$label=~s/\.bam$//; #remove .bam extension
	}
## Start building output
$outHeader=~s/#Abundance/#Abundance\%/ if(defined $scale); #change header to percent if scale given
$output.=$outHeader if(!defined $noHeader); #add header only if noHeader is not specified
$output.="$label"; #add dataset label

## Convert BAM to SAM & count reads mapped
my $samFile=$bamFile.time().int(rand(100)).".sam";
my $samtools=`samtools view -o $samFile $bamFile`; #convert BAM to SAM
my %mappings;
open(my $SAM,$samFile) or die $!;
while(<$SAM>)
	{
	my @line=split(/\t/);
	next if($line[2] eq "*"); #skip unmapped reads
	$mappings{$line[0]} //= 1; #record unique read ids
	}
close($SAM);
unlink($samFile) unless $debug;
my $mapCount=keys %mappings;
if(defined $scale) #convert reads to percent if scale given and add to output
	{
	my $mapPerc=($mapCount/$scale)*100;
	$mapPerc=sprintf("%.3f",$mapPerc);
	$output.="\t$mapPerc";
	}
else #add read numbers to output
	{$output.="\t$mapCount";}
%mappings=() unless $debug; #clear memory taken by read ids

## Calculate per base depth
my $covFile=$bamFile.time().int(rand(100)).".cov";
my $bedtools=`bedtools genomecov -ibam $bamFile -d >$covFile`; #calculate per base depth
my @rScript=`cut -f 3 $covFile |Rscript -e 'x<-as.numeric(readLines(\"stdin\"));summary(x);sd(x);'`; #calculate stats in R
my @stat=split(/\s+/,$rScript[1]);
shift(@stat) if(! $stat[0]); #sometimes R introduces an empty tab if numbers are long
my $sd=(split(/\s+/,$rScript[2]))[1];
#$output.="\t$stat[0]\t$stat[1]\t$stat[2]\t$stat[4]\t$stat[5]\t$stat[3]\t$sd\n";
$output.=sprintf("\t%d\t%d\t%d\t%d\t%d\t%.3f\t%.3f\n",$stat[0],$stat[1],$stat[2],$stat[4],$stat[5],$stat[3],$sd); #add depth statistics to output
unlink($covFile) unless $debug;

## Print everything at the same time
print "$output";
