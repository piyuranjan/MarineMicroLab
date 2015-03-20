#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

sub Info
	{
	print "
#######################################################################
# Script for automatic Quality Control processing for paired end
#	sequencing files in FastQ format
# Runs cutadapt for adapter removal from both ends
# Runs TrimGalore for trimming sequences based on Phred quality scores
# Optionally runs PrinSeq to generate processed read statistics report
#
# Currently encoded adapters/indexes for: Illumina Nextera Index Kit
#
# Author: piyuranjan\@gmail.com
# Created: Mar 20, 2015; Updated: Check on 
#######################################################################
	\n";
	}
sub Usage
	{
	print "
Usage:\n$0 [options]
\nOptions:
 -f|forward		[string:required] provide forward read (read1) file in FastQ.
 -r|reverse		[string:required] provide reverse read (read2) file in FastQ.
 -fa|fadapter	[string:required] provide index code for the forward adapter.
 -ra|radapter	[string:required] provide index code for the reverse adapter.
 -v|verbose		print logging information.
 -h|help		print help (information and usage); exit 0.
	\n";
	}

my ($forwardFile,$reverseFile,$fadapter,$radapter,$help);
if(!GetOptions('f|forward=s' => \$forwardFile,
				'r|reverse=s' => \$reverseFile,
				'fa|fadapter=s' => \$fadapter,
				'ra|radapter=s' => \$radapter,
				'h|help' => \$help))
	{Usage;exit 1;} #quit with error code
if($help) #quit with help
	{Info;Usage;exit 0;}
my $prinSeqStatus=CheckPackages;
if((!defined $forwardFile)||(!defined $reverseFile)||(!defined $fadapter)||(!defined $radapter)) #all of these options are necessary
	{Usage;exit 1;}


#####################################
######### Subroutines Below #########
#####################################

sub CheckPackages #checks the availability of packages required
	{
	print "\n### Checking dependencies ###\n" if $verbose;
	
	##Check Cutadapt & TrimGalore
	`cutadapt -h 1>/dev/null 2>&1`;
	die "Fatal error: cutadapt not installed (or not available in PATH)\n$!" if($?);
	`trim_galore -h 1>/dev/null 2>&1`;
	die "Fatal error: trim_galore not installed (or not available in PATH)\n$!" if($?);
	
	##Check PrinSeq(Lite)
	`prinseq-lite.pl 1>/dev/null 2>&1`;
	my $prinSeqStatus=1;
	if($?)
		{
		warn "Warning: prinseq-lite.pl not installed (or not available in PATH)\ngeneration of reports will be skipped";
		$prinSeqStatus=0;
		}
	print "All packages working fine!\n" if $verbose;
	return $prinSeqStatus;
	}