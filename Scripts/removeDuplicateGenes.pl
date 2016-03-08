#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

sub Info
	{
	print "
#######################################################################
# Script to remove duplicate gene prediction calls in a dataset
# Needs the contig sequences in FastA and gene predictions in GFF2 format
# Produces a new GFF2 file with unique gene calls and a FastA file for
#	predicted gene sequences.
# Author: piyuranjan\@gmail.com
#######################################################################
	\n";
	}
sub Usage
	{
	print "
Usage:\n$0 [options] >[GFF2FormatUniqueGeneCalls]
\nOptions:
 -c|contig	[string:required]	provide name for the FastA format contig file.
 -g|gff		[string:required]	provide name for the GFF2 format gene call file.
 -o|outGff	[string:required]	provide name for the file to output unique gene record in GFF2.
 -f|fasta	[string:optional]	provide name for the file to output unique gene NT sequences in FastA.
	\n";
	}

my ($contigFile,$gffFile,$outGffFile,$help);
my $outFastaFile=0;
my %contigLength;
my %geneReferences;
if(!GetOptions('c|contig=s' => \$contigFile,
				'g|gff=s' => \$gffFile,
				'o|outGff=s' => \$outGffFile,
				'f|fasta=s' => \$outFastaFile,
				'h|help' => \$help))
	{Usage;exit 1;} #quit with error code
if($help) #quit with help
	{Info;Usage;exit 0;}
if((!defined $contigFile)&&(!defined $gffFile)&&(!defined $outGffFile)) #these options are necessary
	{Usage;exit 1;}

$contigFile=SingleLineFasta($contigFile);

####To open GFF file to capture the sequences and gene header in a hash###
open(GFF,$gffFile) or die $!;
open(OGFF,">$outGffFile") or die $!;
while(<GFF>)
	{
	if(/^##/){print OGFF $_; next;} #print all headers to output GFF
	next if(/(^#|^$)/); #skip all comments and blank lines; This line should always come after the header lines are parsed.
	my $gffLine=$_;
	my ($geneSequence,$contigLength)=GetSeqFromGFF($gffLine,$contigFile);
	
	if(defined $geneReferences{$geneSequence}{gff}) #compare contig lengths for duplicate genes and skip smaller source contig
		{next if($contigLength<=$geneReferences{$geneSequence}{contigLength});}
	
	#store gene information in the hash
	$geneReferences{$geneSequence}{gff}=$gffLine; #store gff entry
	$geneReferences{$geneSequence}{contigLength}=$contigLength; #store length of the source contig
	#last if $.==5;
	}
close(GFF);
if($outFastaFile){open(OFA,">$outFastaFile") or die $!;} #open FastA file if option provided
#print all the unique gene information in GFF (and FastA if option given)
foreach my $geneSequence(keys(%geneReferences))
	{
	print OGFF $geneReferences{$geneSequence}{gff};
	if($outFastaFile)
		{
		my $header=">".$geneReferences{$geneSequence}{gff};
		$header=~s/\t/|/g;
		print OFA $header;
		print OFA $geneSequence."\n";
		}
	}
if($outFastaFile){close(OFA);}
close(OGFF);
unlink $contigFile; #remove the temp Contig sequence file

#########################
### Subroutines below ###
#########################

sub SingleLineFasta #convert multi-line fasta file in to a temporary single line fasta file and record sequence length for each contig
	{
	my $contigFile=$_[0]; #fileName of the original contig file
	my $tempContig=time().'.tempContig.fa'; #a temp file with contigs in single line FastA
	open(CTG,$contigFile) or die $!;
	open(OUT,">$tempContig") or die $!;
	while(<CTG>)
		{
		next if(/(^$|#)/);
		if(/^>/)
			{
			if($.!=1){print OUT "\n";}
			print OUT $_;
			}
		else{chomp();print OUT $_;}
		}
	close(OUT);
	close(CTG);
	return($tempContig);
	}

sub GetSeqFromGFF #find and extract the gene sequence using the GFF2 record and the contig FastA file
	{
	#order of input: GFF2 record line, name of the contig FastA (single-line sequence) file.
	my $contigFile=$_[1];
	my @gffRecord=split(/\t/,$_[0]);
	#get the contig sequence from fasta
	my $contigSequence=`grep -A 1 '$gffRecord[0]' $contigFile|tail -1`;
	chomp($contigSequence);
	$contigSequence=uc($contigSequence);
	my $contigLength=length($contigSequence);
	#extract gene sequence from contig
	my $geneSequence=substr($contigSequence,($gffRecord[3])-1,($gffRecord[4]-$gffRecord[3]+1));
	if($gffRecord[6] eq '-') #reverse translate
		{$geneSequence=~tr/ATGC/TACG/;$geneSequence=reverse($geneSequence);}
	return ($geneSequence,$contigLength);
	#order of output: sequence of the gene, length of it's source contig.
	}
