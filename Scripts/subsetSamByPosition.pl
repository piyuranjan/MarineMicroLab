#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;
use List::Util qw(min max);
use List::MoreUtils qw(uniq);

sub Info
	{
	print "
#######################################################################
# Script to subset alignments in SAM format by a given set of
# 	coordinates.
# Script requires SAM format data file to subset and a file with contigs
# 	and coordinates which can be used to subset.
# Script will make two SAM format output files (with same headers if
# 	detected in original file), one which contians aligned reads in the
# 	given coordinates (positive set) and other containing the ones which
# 	didn't overlap the coordinates (negative set).
# 
# Author: piyuranjan\@gmail.com
#######################################################################
	\n";
	}
sub Usage
	{
	print "
Usage:\n$0 [options] >[statsFile]
\nOptions:
 -s|sam		[string:required] SAM format dataset file
 -c|coord	[string:required] Coordinates file
				(default expectation: GFF2 format)
 -p|pos		[string:required] File for positive/matching set output
 -n|neg		[string:required] File for negative/non-matching set output
 -f|format	[string:optional] Comma separated list of column numbers
				for finding coordinates. Format order: Contig,Start,End
				Default: 1,4,5 (for GFF2)
 -h|help	show help and exit 0.
	\n";
	}

my ($samFile,$coordFile,$posFile,$negFile,$help);
my $format="1,4,5";
if(!GetOptions(	's|sam=s' => \$samFile,
				'c|coord=s' => \$coordFile,
				'p|pos=s' => \$posFile,
				'n|neg=s' => \$negFile,
				'f|format=s' => \$format,
				'h|help' => \$help))
	{Usage;exit 1;} #quit with error code
if($help) #quit with help
	{Info;Usage;exit 0;}
if((!defined $samFile)||(!defined $coordFile)||(!defined $posFile)||(!defined $negFile)) #these options are necessary
	{Usage;exit 1;}

##Define columns based on the format definition
my ($contigCol,$startCol,$endCol)=split(/,/,$format);
$contigCol--;$startCol--;$endCol--; #Adjusting for 0-based offset in Perl

#die "Usage:\n$0 [SAM-alignFile] [refPositions-defaultGFF2] [+veSubsetFileName] [-veSubsetFileName]" if(! defined $ARGV[3]);
#Code currently only uses GFF2 format specification for contig, start and end positions.

##Scan and store reference positions to subset
my (@refs,%coord);
open(my $REF,$coordFile) or die $!;
while(<$REF>)
	{
	next if(/#|^$/); #skip comments and blanks
	#Modify section below to match generic positions.
	my @line=split(/\t/);
	push(@refs,$line[$contigCol]); #pushing ContigName
	#store begin and end in a hash with an arbitrary key to differentiate between records on the same contig.
	$coord{$line[$contigCol]}{$#refs}{begin}=$line[$startCol];
	$coord{$line[$contigCol]}{$#refs}{end}=$line[$endCol];
	$coord{$line[$contigCol]}{$#refs}{readCount}=0;
	}
close($REF);
my @uniqContigs=uniq(@refs);

# foreach my $contig(@uniqContigs)
	# {
	# foreach my $key(sort(keys(%{$coord{$contig}})))
		# {
		# print "\n!$contig\t$key\t$coord{$contig}{$key}{begin}\t$coord{$contig}{$key}{end}\t:\t$coord{$contig}{$key}{readCount}!";
		# }
	# print "\n";
	# }

##Open SAM file and estimate overlap of SAM records with coord provided
##If a match, place the SAM in the +vs set
##If not, place the SAM in the -ve set.
open(my $SAM,$samFile) or die $!;
open(my $POS,">$posFile") or die $!;
open(my $NEG,">$negFile") or die $!;
while(<$SAM>)
	{
	if(/^@/) #If header, print to both files
		{
		for my $fh ($POS,$NEG){print $fh $_;}
		next;
		}
	my $samLine=$_;
	my @samEntry=split(/\t/,$samLine);
	if (grep {$_ eq $samEntry[2]} @uniqContigs) #If ContigName match found
		{
		my $printFlag=0;
		foreach my $key(sort(keys(%{$coord{$samEntry[2]}}))) #For each set of coord provided on the same contig
			{
			my $begin=$coord{$samEntry[2]}{$key}{begin};
			my $end=$coord{$samEntry[2]}{$key}{end};
			my $start=$samEntry[3];
			my $stop=$start+length($samEntry[9])-1;
			#print "!!$begin\t$end\t$start\t$stop!!\n";
			my $overlap=min($stop,$end)-max($start,$begin); #Calculate overlap of the read to ref
			if($overlap>0) #Match overlap condition
				{
				print $POS $samLine;
				$printFlag=1;
				$coord{$samEntry[2]}{$key}{readCount}++;
				last; #If successful, break after reporting the SAM record printed to +ve file
				}
			}
		if($printFlag==0) #If none of the refs on this contig matched, meaning read mapped to same contig but a different location than the coord given.
			{print $NEG $samLine;}
		}
	else
		{print $NEG $samLine;}
	}
close($SAM);
close($POS);
close($NEG);

##Print stats for each coord
print "Read stats for positive set:\nEntry#\tContig\tStart\tEnd\t:\tReadCounts";
foreach my $contig(@uniqContigs)
	{
	foreach my $key(sort(keys(%{$coord{$contig}})))
		{
		print "\n$key\t$contig\t$coord{$contig}{$key}{begin}\t$coord{$contig}{$key}{end}\t:\t$coord{$contig}{$key}{readCount}";
		}
	print "\n";
	}
