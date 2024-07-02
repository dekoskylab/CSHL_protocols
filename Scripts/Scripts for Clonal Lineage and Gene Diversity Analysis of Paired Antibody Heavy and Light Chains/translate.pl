#! /opt/local/bin/perl
use warnings;
use strict;
use Text::Wrap;
use Tie::IxHash;

## Program Info:
#
# Name: trans
#
# Purpose: 
#
# Author: John Nash
#
# Copyright (c) Government of Canada, 2000-2012,
#   all rights reserved.
#
# Licence: This script may be used freely as long as no fee is charged 
#   for use, and as long as the author/copyright attributions 
#   are not removed.
#
# History:
# 01 June 2001: Initial alpha
# 31 January 2012: Munged for new home.
#
# Usage:
#
##
#Function:  Translates fasta files of nucleotide sequence to peptide sequences.  The fasta file can be a multiple fasta file. The reading frame is specified on #the command line.  If the frame "-f" parameter is omitted, translate defaults to reading frame "+1".
#Syntax:  translate -f frame {1 2 3 -1 -2 -3} < fasta_file
#   or    translate -h for help
#Input is captured from stdin, so the program can be a pipe.
#Output is written to stdout, so it can be captured to a file or piped; e.g. translate -f 1  < sequence.fasta > peptide.fasta
#Warning:  No attempt is made to ensure that the input file is a valid multiple FASTA file.




my ($title) = "translate";
my ($version) = "1.0";
my ($date) = "31 January 2012";
$Text::Wrap::columns = 65;

# Error message:
my $error_msg = "Type \"$title -h\" for help.";

# Get and process the command line params:
# Returns array of $fasta_file and $orf_file;
my (@cmd_line, $frame);

@cmd_line = ();
@cmd_line = process_command_line();
$frame = $cmd_line[0];

# validate $frame:
if (($frame > 3) || ($frame < -3) || ($frame == 0)) {
    die ("Error: Invalid frame value. Use: 1, 2, 3, -1, -2 or -3\n", 
	 $error_msg, "\n");
}

## Read in the sequence from a multiple FASTA file:
my (%seq_name, %seq_str);
my ($count) = 0;

tie %seq_name, "Tie::IxHash";
tie %seq_str, "Tie::IxHash";

# read the header:
while (<>) {
    s/\r\n/\n/g;
    chomp;
    if (/^>/)  {
	$seq_name{$count} = substr($_, 1, length $_);
	$count++;
    }
    else {
	$seq_str{$count - 1} .= $_;
    }
}

foreach my $item (keys %seq_name) {
    
# send it off for translation:
    my $peptide = translate ($seq_str{$item}, $frame);
    
# print it:
    print ">$seq_name{$item}\n";
    print wrap ('', '', "$peptide\n");
}

### end of main:

## Subroutines:
sub translate  {
    
# Does:  Translates incoming string from nucleic acid to peptide
# Expects: nucleic acid string ($sequence_str); 
#          reading frame value: +1, +2, +3, -1, -2, -3
# Returns: string of translated peptide
# Uses:
    
# vars for translation:
    my ($codon, $peptide, $seq_str, $frame);
    
# reading frame:
    my $startpoint;
    
# incoming:
    $seq_str = uc $_[0];
    
# assess the reading frame:
    $frame = $_[1];
    $seq_str =~ tr/natgcbdkrvhmyxwsNATGCBDKRVHMYXWS
    /ntacgvhmybdkrxswNTACGVHMYBDKRXSW/ if ($frame < 0);
    $startpoint = (abs $frame) - 1;
    
# Perform the translation:
    for (my $j=$startpoint; $j < length $seq_str; $j+=3) {
	$codon=(substr($seq_str, $j, 3));
	if (($codon eq "GCA") || 
	    ($codon eq "GCC") ||
	    ($codon eq "GCG") ||
	    ($codon eq "GCT") ||
	    ($codon eq "GCR") ||
	    ($codon eq "GCY") ||
	    ($codon eq "GCN")) {
	    $peptide.="A";
	}   # Ala A
	
	elsif (($codon eq "CGA") || 
	       ($codon eq "CGC") ||
	       ($codon eq "CGG") ||
	       ($codon eq "CGT") ||
	       ($codon eq "CGR") ||
	       ($codon eq "CGY") ||
	       ($codon eq "CGN") ||
	       ($codon eq "AGA") ||
	       ($codon eq "AGG") ||
	       ($codon eq "AGR")) {
	    $peptide.="R";
	}   # Arg R
	
	elsif (($codon eq "AAC") || 
	       ($codon eq "AAT") ||
	       ($codon eq "AAY")) {
	    $peptide.="N";
	}   # Asn N
	
	elsif (($codon eq "GAC") || 
	       ($codon eq "GAT") ||
	       ($codon eq "GAY")) {
	    $peptide.="D";
	}   # Asp D
	
	elsif (($codon eq "TGC") || 
	       ($codon eq "TGT") ||
	       ($codon eq "TGY")) {
	    $peptide.="C";
	}   # Cys C
	
	elsif (($codon eq "CAA") || 
	       ($codon eq "CAG") ||
	       ($codon eq "CAR")) {
	    $peptide.="Q";
	}   # Gln Q
	
	elsif (($codon eq "GAA") || 
	       ($codon eq "GAG") ||
	       ($codon eq "GAR")) {
	    $peptide.="E";
	}   # Glu E
	
	elsif (($codon eq "GGA") || 
	       ($codon eq "GGC") ||
	       ($codon eq "GGG") ||
	       ($codon eq "GGT") ||
	       ($codon eq "GGR") ||
	       ($codon eq "GGY") ||
	       ($codon eq "GGN")) {
	    $peptide.="G";
	}   # Gly G
	
    elsif (($codon eq "CAC") || 
	   ($codon eq "CAT") ||
	   ($codon eq "CAT")) {
	$peptide.="H";
    }   # His H
	
	elsif (($codon eq "ATA") || 
	       ($codon eq "ATC") ||
	       ($codon eq "ATT") ||
	       ($codon eq "ATH")) {
	    $peptide.="I";
	}   # Ile I
	
	elsif (($codon eq "CTA") || 
	       ($codon eq "CTC") ||
	       ($codon eq "CTG") ||
	       ($codon eq "CTT") ||
	       ($codon eq "CTR") ||
	       ($codon eq "CTY") ||
	       ($codon eq "CTN") ||
	       ($codon eq "TTA") ||
	       ($codon eq "TTG") ||
	       ($codon eq "TTR")) {
	    $peptide.="L";
	}   # Leu L
	
	elsif (($codon eq "AAA") || 
	       ($codon eq "AAG") ||
	       ($codon eq "AAR")) {
	    $peptide.="K";
	}   # Lys K
	
	elsif (($codon eq "ATG")) {
	    $peptide.="M";
	}   # Met M
	
	elsif (($codon eq "TTC") || 
	       ($codon eq "TTT") ||
	       ($codon eq "TTY")) {
	    $peptide.="F";
	}   # Phe F
	
	elsif (($codon eq "CCA") || 
	       ($codon eq "CCC") ||
	       ($codon eq "CCG") ||
	       ($codon eq "CCT") ||
	       ($codon eq "CCR") ||
	       ($codon eq "CCY") ||
	       ($codon eq "CCN")) {
	    $peptide.="P";
	}   # Pro P
	
	elsif (($codon eq "TCA") || 
	       ($codon eq "TCC") ||
	       ($codon eq "TCG") ||
	       ($codon eq "TCT") ||
	       ($codon eq "TCR") ||
	       ($codon eq "TCY") ||
	       ($codon eq "TCN") ||
	       ($codon eq "AGC") ||
	       ($codon eq "AGT") ||
	       ($codon eq "AGY")) {
	    $peptide.="S";
	}   # Ser S
	
	elsif (($codon eq "ACA") || 
	       ($codon eq "ACC") ||
	       ($codon eq "ACG") ||
	       ($codon eq "ACT") ||
	       ($codon eq "ACR") ||
	       ($codon eq "ACY") ||
	       ($codon eq "ACN")) {
	    $peptide.="T";
	}   # Thr T
	
    elsif (($codon eq "TGG")) {
	$peptide.="W";
    }   # Trp W
	
	elsif (($codon eq "TAC") || 
	       ($codon eq "TAT") ||
	       ($codon eq "TAY")) {
	    $peptide.="Y";
	}   # Tyr Y
	
	elsif (($codon eq "GTA") || 
	       ($codon eq "GTC") ||
	       ($codon eq "GTG") ||
	       ($codon eq "GTT") ||
	       ($codon eq "GTR") ||
	       ($codon eq "GTY") ||
	       ($codon eq "GTN")) {
	    $peptide.="V";
	}   # Val V
	  
	elsif (($codon eq "TAA") || 
	       ($codon eq "TAG") ||
	       ($codon eq "TAR") ||
	       ($codon eq "TGA")) {
	    $peptide.="*";
	}   # Stop *
	
	else {
	    ;
	} # do nothing for now...
    } # end of  for ($j=0; etc...
    
# return the value:
    return $peptide;
} # end of sub translate

sub process_command_line {
#
# Expects: 
# Returns: @my_cmd_line = ($frame)
# Uses:
    
# Variables:
    my %opts = ();    # command line params, as entered by user
    my @my_cmd_line;  # returned value
    my @list;         # %opts as an array for handling
    my $cmd_args;	    # return value for getopts()
    
# Holders for command line's files:
    my $frame = 1;
    
# Scratch:
    my $item;
    
# Get the command=line parameters:
    use vars qw($opt_f $opt_h);
    use Getopt::Std;
    $cmd_args = getopts('f:h', \%opts);
    
# Die on illegal argument list:
    if ($cmd_args == 0) {
	die ("Error: Missing or incorrect command line parameter(s)!\n",
	     $error_msg, "\n");
    }
    
# Check and enact each command-line argument:
    if (!%opts)  {
	die ($error_msg, "\n");
    }
    
# Make the hashes into an array:
    @list = keys %opts;
    
# Do a quick check for "help" and the compulsory parameters:
#   If the compulsory files are not there, squawk and die:
    foreach $item (@list)  {
# Help:
	if ($item eq "h")  {
	    help();
	}
# frame:
	elsif ($item eq "f") {
	    $frame = $opts{$item};
	}
    }
    
# Put it in an array:
    @my_cmd_line = ($frame);
    return @my_cmd_line;
    
} #end of sub process_command_line()

sub help {
    
print <<EOHelp;
$title, version $version, $date

Function:  Translates fasta files of nucleotide sequence to peptide 
   sequences.  The fasta file can be a multiple fasta file.
   The reading frame is specified on the command line.  If the frame
   \"-f\" parameter is omitted, $title defaults to reading frame \"+1\".

Syntax:  $title -f frame {1 2 3 -1 -2 -3} < fasta_file
   or    $title -h for help

Input is captured from stdin, so the program can be a pipe.

    Output is written to stdout, so it can be captured to a file or piped;
   e.g. $title -f 1  < sequence.fasta > peptide.fasta

Warning:  No attempt is made to ensure that the input file is a valid 
   multiple FASTA file.

EOHelp
die ("\n");
} # end of sub help
