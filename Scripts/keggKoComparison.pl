#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

sub Info
	{
	print "
#######################################################################
# Script to compare KEGG KO identifiers between annotated genomes
#	or metagenomes. The script uses KO# and definitions from the
#	annotated genomes to produce a unique list of KO and their counts
#	in the annotated genomes.
# Genomes/metagenomes should be annotated using the
#	BlastKOALA/GhostKOALA with the KEGG KO assignment (including
#	description) stored in a flat tab separated text file.
#	Order of cols: GeneID, KEGGKO#, KODefinition.
# Filenames for multiple genome annotations (at least two) should be
#	given as a single string delimited by comma (,) and no spaces.
#	Additionally, labels can be provided in the same fashion and order.
# The output of the script is the copy number of each KO category in
#	each genome compared (aligned) between all the genomes. In the end,
#	it will also provide how many genomes had each KO category.
#	Output will look like a table in a tab separated flat file.
#
# Author: piyuranjan\@gmail.com
#######################################################################
	\n";
	}
sub Usage
	{
	print "
Usage:\n$0 [options] >[outputFile]
\nOptions:
 -g|genomes	[string:required] list of all genome annotation files
					delimited by comma (,) and no spaces. Minimum 2
					files required.
 -l|labels	[string:optional] list of labels to use, instead of the
					annotation filenames (Default), in the same order
					delimited by comma (,) and no spaces.
 -h|help	show help and exit 0.
\nExample:\n$0 -g g1.txt,g2.txt,g3.txt -l Genome1,Genome2,Genome3
	\n";
	}

my ($fileList,$labelList,$help);
if(!GetOptions('g|genomes=s' => \$fileList,
				'l|labels=s' => \$labelList,
				'h|help' => \$help))
	{Usage;exit 1;} #quit with error code
if($help) #quit with help
	{Info;Usage;exit 0;}
if(!defined $fileList) #these options are necessary
	{Usage;exit 1;}

#Store annotation files and their labels
my @files=split(/,/,$fileList);
if($#files<1) #check if less than two files provided
	{
	Usage;
	print STDERR "\nError: need more than one annotated genome files\nWorking for a single genome is useless\n";
	exit 1;
	}
if(!defined $labelList)
	{$labelList=$fileList;}
my @labels=split(/,/,$labelList);
die "\nError: # labels don't match # files\n" if($#labels!=$#files);
#Store fileNames and labels as a hash
my %genomes;
foreach my $key(0..$#labels)
	{
	die "\nError: $files[$key]\n$!" if(! -f $files[$key]); #Die if annotation file is not found (/not a plain file)
	$genomes{$labels[$key]}=$files[$key];
	}
#print %labels;

##Scanning all files/genomes for KO categories
my %ko;
foreach my $label(sort(keys(%genomes))) #Do for all genomes
	{
	open(INP,$genomes{$label}) or die $!;
	while(<INP>)
		{
		chomp;
		my @annotation=split(/\t/);
		next if(!defined $annotation[1]); #Skip if no KO annotation available
		if(!defined $ko{$annotation[1]}{definition}) #Add definition to KO if not available
			{$ko{$annotation[1]}{definition}=$annotation[2];}
		if(!defined $ko{$annotation[1]}{$label}) #Initialize KO count for the genome
			{$ko{$annotation[1]}{$label}=0;}
		$ko{$annotation[1]}{$label}++; #Add/Register occuence of KO in the genome
		#$ko{$annotation[1]}{$label}
		}
	close(INP);
	}

##Record consistency of KO for genomes and Print output in a tabular format simultaneously
#prepare headers
print "KEGG-KO#\tKODescription";
foreach my $label(sort(keys(%genomes)))
	{print "\t$label";}
print "\tConsistency";
#record consistency
foreach my $category(sort(keys(%ko)))
	{
	print "\n$category\t$ko{$category}{definition}"; #print KO# and definition
	my $count=0;
	foreach my $label(sort(keys(%genomes)))
		{
		if(!defined $ko{$category}{$label}) #Record absence of a KO for a genome
			{$ko{$category}{$label}=0;}
		else #Record consistency of the KO
			{$count++;}
		print "\t$ko{$category}{$label}"; #print KO occurence for a genome
		}
	$ko{$category}{consistency}=$count;
	print "\t$ko{$category}{consistency}"; #print consistency value of a KO for all genomes
	}




#########################
### Subroutines below ###
#########################


