#! /usr/bin/perl
# This script is to identify the CDR3 for human and to report the associated information to a tab-delimited file
# This script is useful for NGS data for pairing experiments as this will try to link them based on the header label
# Also, this will use the motif search so limited sequence length/info will not affect results too much
# Note: this version adds the IgBlast V gene usage assignments
# Note: this version will also adds the isotype identified info using fastx_barcode_splitter.pl and also it requires the barcodes.txt
# Note: the CDR3 reported will contain the Cystein C and Tryptophan W; a -1 and a +1 AA surrounding the CDR3
# Note: you still need internal_data from the IgBlast folder to locate to your working directory
# The output file contains the following layout:
# "Header\tCDRH3\tCDRL3\tH3_nt\tL3_nt\tCDRH3_junc_len\tCDRL3_junc_len\tVDJ_nt\tVDJ_aa\tVJ_nt\tVJ_aa\tHeavyRead\tLightRead\tR1vgene\tR2vgene\tR1dgene\tR2dgene\tR1jgene\tR2jgene\tHIsotype\tLIsotype\n"
# Usage: perl vpairhumanalysis.pl FASTA_R1.fna FASTA_R2.fna barcodes.txt

use Tie::File;
use Cwd qw(abs_path);
use File::Basename qw(dirname);


$infile1=$ARGV[0]; #R1 reads
$infile2=$ARGV[1]; #R2 reads
$infile3=$ARGV[2]; #The barcode file
$infile1 =~ m/^(.*)\.(.*)/ ;
$filename=$1;

# Read in isotype labels from barcodes.txt
@barcodes=();
open(InBarcodes,$infile3) or die "Error opening barcode file $infile3 !\n";
while ($line=<InBarcodes>) {
	chomp($line);
	@tmp=split(/\t/,$line);
	push(@barcodes,$tmp[0]);
}
push(@barcodes,'unmatched');
close InBarcodes;

# Run fastx_barcode_splitter.pl to split files to barcodes and then merge with incorporation of the isotype info to header
print "Running barcode splitter with R1 reads\n";
system("cat $infile1 | fastx_barcode_splitter.pl --bcfile $infile3 --prefix $infile1"."R1_ --bol --mismatches 12 --partial 8 --quiet");
open(OutTemp,">$infile1-R1reads.fna") or die $!;
foreach $isotype (@barcodes) {
	open(InTemp,"$infile1"."R1_$isotype") or die "Error opening $infile1_$isotype !\n";
	while($line=<InTemp>){
		chomp($line);
		if ($line =~ />/) {
			print OutTemp "$line<$isotype\n";
		}
		else {
			print OutTemp "$line\n";
		}
	}
}
system("rm $infile1"."R1_*");
close OutTemp;

print "Running barcode splitter with R2 reads\n";
system("cat $infile2 | fastx_barcode_splitter.pl --bcfile $infile3 --prefix $infile2"."R2_ --bol --mismatches 12 --partial 8 --quiet");
open(OutTemp,">$infile2-R2reads.fna") or die $!;
foreach $isotype (@barcodes) {
	open(InTemp,"$infile2"."R2_$isotype") or die "Error opening $infile2_$isotype !\n";
	while($line=<InTemp>){
		chomp($line);
		if ($line =~ />/) {
			print OutTemp "$line<$isotype\n";
		}
		else {
			print OutTemp "$line\n";
		}
	}
}
system("rm $infile2"."R2_*");
close OutTemp;

# Generate the IgBlast results
print "Generating IgBlast results\n";
$igblastpre=dirname(abs_path($0)).'/./igblastn -germline_db_V '.dirname(abs_path($0)).'/database/human_gl_V -germline_db_D '.dirname(abs_path($0)).'/database/human_gl_D -germline_db_J '.dirname(abs_path($0)).'/database/human_gl_J -auxiliary_data '.dirname(abs_path($0)).'/optional_file/human_gl.aux -domain_system imgt -num_alignments_V 4 -num_alignments_D 4 -num_alignments_J 4 -outfmt 3 -query ';
$igblastoutput1=$infile1.'-R1reads-igblast.txt';
$igblastoutput2=$infile2.'-R2reads-igblast.txt';
$igblastcmd1=$igblastpre.$infile1.'-R1reads.fna > '.$igblastoutput1;
$igblastcmd2=$igblastpre.$infile2.'-R2reads.fna > '.$igblastoutput2;
system($igblastcmd1);
system($igblastcmd2);

# Generate dictionary storing best hit V gene usage
tie @OutIgBlast1,'Tie::File',$igblastoutput1;
%R1vgene=[];
%R1dgene=[];
%R1jgene=[];
for ($i=0;$i<@OutIgBlast1;$i++) {
	$line=$OutIgBlast1[$i];
	chomp($line);
	if ($line =~ m/^Query= (.*)\s.{1}:N:/) {
		$fnaheader=$1;
	}
	if ($line =~ m/V\(D\)J rearrangement summary for query sequence/) {
		$line=$OutIgBlast1[$i+1];
		@tmp=split(/\t/,$line);
		if (@tmp==5) {
			$R1vgene{$fnaheader}=$tmp[0];
			$R1dgene{$fnaheader}='-';
			$R1jgene{$fnaheader}=$tmp[1];
		}
		else {
			$R1vgene{$fnaheader}=$tmp[0];
			$R1dgene{$fnaheader}=$tmp[1];
			$R1jgene{$fnaheader}=$tmp[2];
		}
	}
}
untie @OutIgBlast1;

tie @OutIgBlast2,'Tie::File',$igblastoutput2;
%R2vgene=[];
%R2dgene=[];
%R2jgene=[];
for ($i=0;$i<@OutIgBlast2;$i++) {
	$line=$OutIgBlast2[$i];
	chomp($line);
	if ($line =~ m/^Query= (.*)\s.{1}:N:/) {
		$fnaheader=$1;
	}
	if ($line =~ m/V\(D\)J rearrangement summary for query sequence/) {
		$line=$OutIgBlast2[$i+1];
		@tmp=split(/\t/,$line);
		if (@tmp==5) {
			$R2vgene{$fnaheader}=$tmp[0];
			$R2dgene{$fnaheader}='-';
			$R2jgene{$fnaheader}=$tmp[1];
		}
		else {
			$R2vgene{$fnaheader}=$tmp[0];
			$R2dgene{$fnaheader}=$tmp[1];
			$R2jgene{$fnaheader}=$tmp[2];
		}
	}
}
untie @OutIgBlast2;

# Amino acid translation hash #Note: X is Ambiguous N and Z is stop codon
%base2aa =  ("AAA" => "K", "AAC" => "N", "AAG" => "K", "AAT" => "N", "ACA" => "T", "ACC" => "T", "ACG" => "T", "ACT" => "T", "AGA" => "R", "AGC" => "S", "AGG" => "R", "AGT" => "S", "ATA" => "I", "ATC" => "I", "ATG" => "M",  "ATT" => "I", "CAA" => "Q", "CAC" => "H", "CAG" => "Q", "CAT" => "H", "CCA" => "P", "CCC" => "P", "CCG" => "P", "CCT" => "P", "CGA" => "R", "CGC" => "R", "CGG" => "R", "CGT" => "R", "CTA" => "L", "CTC" => "L", "CTG" => "L", "CTT" => "L", "GAA" => "E", "GAC" => "D", "GAG" => "E", "GAT" => "D", "GCA" => "A", "GCC" => "A", "GCG" => "A", "GCT" => "A", "GGA" => "G", "GGC" => "G", "GGG" => "G", "GGT" => "G", "GTA" => "V", "GTC" => "V", "GTG" => "V", "GTT" => "V", "TAA" => "Z", "TAC" => "Y", "TAG" => "Z", "TAT" => "Y", "TCA" => "S", "TCC" => "S", "TCG" => "S",  "TCT" => "S", "TGA" => "Z", "TGC" => "C", "TGG" => "W", "TGT" => "C", "TTA" => "L", "TTC" => "F", "TTG" => "L", "TTT" => "F");

print "++++++++Please wait while run is preparing to start++++++++\n";

tie @InFile1,'Tie::File',$infile1.'-R1reads.fna';
open (OutTemp1,">$infile1".'-temp1.txt') or die $!;
# Process the 1st file
$n=0;
#%seqs=();
for ($i=0;$i<@InFile1;$i++) {
	$line=$InFile1[$i];
	if ($line =~ m/^\>(.*)\s(.{1}):N:.*\<(.*)/){
		$header=$1;
		$read=$2;
		$iso=$3;
		$n++;
#		print "Processing R$read ...$n\n";
		next;
	}
	$current=$line;
	$current_r=reverse($current);
	$current_r=~tr/ATGC/TACG/;
	$len=length($current);
	@AA=();@cdna=();
	$cdna[1] = substr ($current, 0); for ($j=0; $j<$len; $j=$j+3) {$triplet = substr ($current, $j, 3); if ($triplet =~ "N") {$aa="X";} else {$aa=$base2aa{$triplet};} $AA[1]=$AA[1].$aa;}
    $cdna[2] = substr ($current, 1); for ($j=1; $j<$len; $j=$j+3) {$triplet = substr ($current, $j, 3); if ($triplet =~ "N") {$aa="X";} else {$aa=$base2aa{$triplet};} $AA[2]=$AA[2].$aa;}
    $cdna[3] = substr ($current, 2); for ($j=2; $j<$len; $j=$j+3) {$triplet = substr ($current, $j, 3); if ($triplet =~ "N") {$aa="X";} else {$aa=$base2aa{$triplet};} $AA[3]=$AA[3].$aa;}
    $cdna[4] = substr ($current_r, 0); for ($j=0; $j<$len; $j=$j+3) {$triplet = substr ($current_r, $j, 3); if ($triplet =~ "N") {$aa="X";} else {$aa=$base2aa{$triplet};} $AA[4]=$AA[4].$aa;}
    $cdna[5] = substr ($current_r, 1); for ($j=1; $j<$len; $j=$j+3) {$triplet = substr ($current_r, $j, 3); if ($triplet =~ "N") {$aa="X";} else {$aa=$base2aa{$triplet};} $AA[5]=$AA[5].$aa;}
    $cdna[6] = substr ($current_r, 2); for ($j=2; $j<$len; $j=$j+3) {$triplet = substr ($current_r, $j, 3); if ($triplet =~ "N") {$aa="X";} else {$aa=$base2aa{$triplet};} $AA[6]=$AA[6].$aa;}
	$max_P_H1=1e-99;     $max_P_H2=1e-99;
	$found=0;
	
	# Motif search #1
	for ($j=1; $j<=6; $j++) { 
		$current_AA=$AA[$j]; 
		$len_AA=length($current_AA);
		if ($len_AA<10) {$best_segment_H1="BLANK";$best_segment_H2="BLANK";}
		else {
			# Left Flank
			for ($k=0; $k<($len_AA-8); $k++) {
				$segment = substr ($current_AA, $k, 8); $p_H=1e-99;
	
				for ($kk=0; $kk<8; $kk++) {$m[$kk+1] = substr ($segment,$kk,1);}
 		
				if ($m[1] eq "D") {$p_H=0.99;} else {$p_H=0.01/19;} 
	
				if ($m[2] eq "T") {$p_H=$p_H*0.96;} elsif ($m[2] eq "S") {$p_H=$p_H*0.015;} elsif (($m[2] eq "A") or ($m[2] eq "M")) {$p_H=$p_H*0.0075;} else {$p_H=$p_H*0.01/16;} 

				if ($m[3] eq "A") {$p_H=$p_H*0.97;} elsif ($m[3] eq "G") {$p_H=$p_H*0.02;} else {$p_H=$p_H*0.01/18;} 

				if ($m[5] eq "Y") {$p_H=$p_H*0.99;} elsif ($m[5] eq "F") {$p_H=$p_H*0.005;} else {$p_H=$p_H*0.005/18;} 

				if ($m[6] eq "Y") {$p_H=$p_H*0.91;} elsif ($m[6] eq "F") {$p_H=$p_H*0.06;} elsif ($m[6] eq "S") {$p_H=$p_H*0.01;} elsif ($m[6] eq "H") {$p_H=$p_H*0.005;} else {$p_H=$p_H*0.015/16;} 
		
				if ($m[7] eq "C") {$p_H=$p_H*0.99;} else {$p_H=$p_H*0.01/19;} 

				if ($m[8] eq "A") {$p_H=$p_H*0.86;} elsif ($m[8] eq "T") {$p_H=$p_H*0.05;} elsif ($m[8] eq "V") {$p_H=$p_H*0.04;} elsif (($m[8] eq "G") or ($m[8] eq "I")) {$p_H=$p_H*0.01;} elsif (($m[8] eq "S") or ($m[8] eq "L")) {$p_H=$p_H*0.005;} else {$p_H=$p_H*0.02/13;} 
 
				if ($p_H>$max_P_H1) {$max_P_H1=$p_H; $loc1_H1=$j; $loc2_H1=$k;} 
			}
			# Right Flank
			for ($k=0; $k<($len_AA-6); $k++) {
				$segment = substr ($current_AA, $k, 6); $p_H=1e-99; 
	
				for ($kk=0; $kk<6; $kk++) {$m[$kk+1] = substr ($segment,$kk,1);}
 		
				if ($m[1] eq "Y") {$p_H=0.40;} elsif ($m[1] eq "V") {$p_H=0.27;} elsif (($m[1] eq "I") or ($m[1] eq "L") or ($m[1] eq "P") or ($m[1] eq "S")) {$p_H=0.06;} elsif (($m[1] eq "F") or ($m[1] eq "H")) {$p_H=0.02;} else {$p_H=0.05/12;}
 
				if ($m[2] eq "W") {$p_H=$p_H*0.99;} else {$p_H=$p_H*0.01/19;} 

				if ($m[3] eq "G") {$p_H=$p_H*0.98;} else {$p_H=$p_H*0.02/19;} 

				if ($m[4] eq "Q") {$p_H=$p_H*0.83;} elsif ($m[4] eq "K") {$p_H=$p_H*0.07;} elsif ($m[4] eq "R") {$p_H=$p_H*0.05;} elsif ($m[4] eq "P") {$p_H=$p_H*0.02;} else {$p_H=$p_H*0.03/16;} 
	
				if ($m[5] eq "G") {$p_H=$p_H*0.98;} else {$p_H=$p_H*0.02/19;} 

				if ($m[6] eq "T") {$p_H=$p_H*0.93} elsif ($m[6] eq "S") {$p_H=$p_H*0.016;} elsif ($m[6] eq "I") {$p_H=$p_H*0.01;} elsif ($m[6] eq "G") {$p_H=$p_H*0.008;} elsif (($m[6] eq "P") or ($m[6] eq "N")) {$p_H=$p_H*0.005;} else {$p_H=$p_H*0.026/14} 

				if ($p_H>$max_P_H2) {$max_P_H2=$p_H; $loc1_H2=$j; $loc2_H2=$k;} 
			}
		}
	}
	
	if (($max_P_H1*$max_P_H2>1e-10) && ($loc1_H1==$loc1_H2) && ($loc2_H2>$loc2_H1) && ($loc2_H2-$loc2_H1-6 > 0)) {
			$found=1;
			$cdr3temp=substr($AA[$loc1_H1],$loc2_H1+6,$loc2_H2-$loc2_H1-4);
			$aatemp=$AA[$loc1_H1];
			$cdnatemp=$cdna[$loc1_H1];
			$pos=$read-1;
			$cdr3len=length($cdr3temp);
			$cdr3cdna=substr($cdnatemp,($loc2_H1+6)*3,($loc2_H2-$loc2_H1-4)*3);
			print OutTemp1 "H\t$header,$cdr3temp\t$cdr3cdna\t$cdr3len\t$cdnatemp\t$aatemp\tR1H\t$iso\n";
			#$seqs{$header}[0]=$cdr3temp;
			#$seqs{$header}[$pos+2]=$cdr3temp;
			#$seqs{$header}[$pos+3]=$cdnatemp;
			#$seqs{$header}[$pos+4]=$aatemp;
	}
	
	# Motif search #2
	if ($found==0) {
		$max_P_H1=1e-99;     $max_P_H2=1e-99; 
		for ($j=1; $j<=6; $j++) { 
			$current_AA=$AA[$j]; 
			$len_AA=length($current_AA);
			if ($len_AA<10) {$best_segment_H1="BLANK";$best_segment_H2="BLANK";}
			else {
				# Left Flank
				for ($k=0; $k<($len_AA-10); $k++) {
					$segment = substr ($current_AA, $k, 10); $p_H=1e-99;
	
					for ($kk=0; $kk<10; $kk++) {$m[$kk+1] = substr ($segment,$kk,1);}
 		
					if ($m[1] eq "E") {$p_H=0.91;} elsif ($m[1] eq "D") {$p_H=0.07;} elsif ($m[1] eq "A" || $m[1] eq "G") {$p_H=0.008/2;} else {$p_H=0.012/16;}
	
					if ($m[2] eq "D") {$p_H=$p_H*0.99;} else {$p_H=$p_H*0.01/19;} 

					if ($m[6] eq "Y") {$p_H=$p_H*0.99;} else {$p_H=$p_H*0.01/19;} 

					if ($m[7] eq "Y") {$p_H=$p_H*0.93;} elsif ($m[7] eq "F") {$p_H=$p_H*0.05;} else {$p_H=$p_H*0.02/18;} 

					if ($m[8] eq "C") {$p_H=$p_H*0.995;} else {$p_H=$p_H*0.005/19;} 
		
					if ($m[9] eq "Q") {$p_H=$p_H*0.83;} elsif ($m[9] eq "M") {$p_H=$p_H*0.10;} elsif ($m[9] eq "L" || $m[9] eq "H") {$p_H=$p_H*0.028/2;} else {$p_H=$p_H*0.042/16;}

					if ($m[10] eq "Q") {$p_H=$p_H*0.87;} elsif ($m[10] eq "H") {$p_H=$p_H*0.066;} elsif ($m[10] eq "K") {$p_H=$p_H*0.015;} else {$p_H=$p_H*0.049/17;} 
			 
					if ($p_H>$max_P_H1) {$max_P_H1=$p_H; $loc1_H1=$j; $loc2_H1=$k;} 
				}
				# Right Flank
				for ($k=0; $k<($len_AA-6); $k++) {
					$segment = substr ($current_AA, $k, 6); $p_H=1e-99; 
	

					for ($kk=0; $kk<6; $kk++) {$m[$kk+1] = substr ($segment,$kk,1);}
 		
					if ($m[1] eq "T") {$p_H=0.88;} elsif ($m[1] eq "A") {$p_H=0.05;} elsif ($m[1] eq "S") {$p_H=0.034;} else {$p_H=0.036/17;}
 
					if ($m[2] eq "F") {$p_H=$p_H*0.99;} else {$p_H=$p_H*0.01/19;} 

					if ($m[3] eq "G") {$p_H=$p_H*0.99;} else {$p_H=$p_H*0.01/19;} 

					if ($m[4] eq "Q") {$p_H=$p_H*0.54;} elsif ($m[4] eq "G") {$p_H=$p_H*0.26;} elsif ($m[4] eq "P") {$p_H=$p_H*0.18;} else {$p_H=$p_H*0.02/17;} 
	
					if ($m[5] eq "G") {$p_H=$p_H*0.997;} else {$p_H=$p_H*0.003/19;} 

					if ($m[6] eq "T") {$p_H=$p_H*0.99} else {$p_H=$p_H*0.01/19;} 

					if ($p_H>$max_P_H2) {$max_P_H2=$p_H; $best_segment_H2=$segment; $loc1_H2=$j; $loc2_H2=$k; }
				}
			}
		}
	
		if (($max_P_H1*$max_P_H2>1e-16) && ($loc1_H1==$loc1_H2) && ($loc2_H2>$loc2_H1) && ($loc2_H2-$loc2_H1-6 > 0)) {
				$found=1;
				$cdr3temp=substr($AA[$loc1_H1], $loc2_H1+7, $loc2_H2-$loc2_H1-5);
				$aatemp=$AA[$loc1_H1];
				$cdnatemp=$cdna[$loc1_H1];
				$pos=$read-1;
				$cdr3len=length($cdr3temp);
				$cdr3cdna=substr($cdnatemp,($loc2_H1+7)*3,($loc2_H2-$loc2_H1-5)*3);
				print OutTemp1 "L\t$header,$cdr3temp\t$cdr3cdna\t$cdr3len\t$cdnatemp\t$aatemp\tR1L\t$iso\n";
				# $seqs{$header}[1]=$cdr3temp;
				# $seqs{$header}[$pos+2]=$cdr3temp;
				# $seqs{$header}[$pos+3]=$cdnatemp;
				# $seqs{$header}[$pos+4]=$aatemp;
		}
	}	
	
	# Motif search #3
	if ($found==0) {
		$max_P_H1=1e-99;     $max_P_H2=1e-99; 
		for ($j=1; $j<=6; $j++) { 
			$current_AA=$AA[$j]; 
			$len_AA=length($current_AA);
			if ($len_AA<10) {$best_segment_H1="BLANK";$best_segment_H2="BLANK";}
			else {
				# Left Flank
				for ($k=0; $k<($len_AA-8); $k++) {
					$segment = substr ($current_AA, $k, 8); $p_H=1e-99;
	
					for ($kk=0; $kk<8; $kk++) {$m[$kk+1] = substr ($segment,$kk,1);}
 		
					if ($m[1] eq "E") {$p_H=0.67;} elsif ($m[1] eq "G") {$p_H=0.17;} elsif ($m[1] eq "D") {$p_H=0.077;} elsif ($m[1] eq "M") {$p_H=0.059;} else {$p_H=0.024/16;}
	
					if ($m[2] eq "D") {$p_H=$p_H*0.99;} else {$p_H=$p_H*0.01/19;} 

					if ($m[3] eq "E") {$p_H=$p_H*0.98;} else {$p_H=$p_H*0.02/19;} 

					if ($m[4] eq "A") {$p_H=$p_H*0.93;} elsif ($m[4] eq "S") {$p_H=$p_H*0.04;} elsif ($m[4] eq "G") {$p_H=$p_H*0.02;} else {$p_H=$p_H*0.01/17;} 

					if ($m[5] eq "D") {$p_H=$p_H*0.92;} elsif ($m[5] eq "E") {$p_H=$p_H*0.038;} elsif ($m[5] eq "V") {$p_H=$p_H*0.012;} else {$p_H=$p_H*0.03/17;} 
		
					if ($m[6] eq "Y") {$p_H=$p_H*0.99;} elsif ($m[6] eq "F") {$p_H=$p_H*0.0095;} else {$p_H=$p_H*5e-4/18;}

					if ($m[7] eq "Y") {$p_H=$p_H*0.92;} elsif ($m[7] eq "F") {$p_H=$p_H*0.051;} elsif ($m[7] eq "H") {$p_H=$p_H*0.022;} else {$p_H=$p_H*7e-3/17;} 
			 
					if ($m[8] eq "C") {$p_H=$p_H*0.985;} elsif ($m[8] eq "R") {$p_H=$p_H*0.0079;} else {$p_H=$p_H*7.1e-3/18;}
			 
					if ($p_H>$max_P_H1) {$max_P_H1=$p_H; $best_segment_H1=$segment; $loc1_H1=$j; $loc2_H1=$k; }  
				}
				# Right Flank
				for ($k=0; $k<($len_AA-6); $k++) {
					$segment = substr ($current_AA, $k, 6); $p_H=1e-99; 
	
					for ($kk=0; $kk<6; $kk++) {$m[$kk+1] = substr ($segment,$kk,1);}
 		
					if ($m[1] eq "V") {$p_H=0.79;} elsif ($m[1] eq "I") {$p_H=0.076;} elsif ($m[1] eq "L") {$p_H=0.056;} elsif ($m[1] eq "A" || $m[1] eq "G" || $m[1] eq "M") {$p_H=0.014/3;} else {$p_H=0.064/14;}
 
					if ($m[2] eq "F") {$p_H=$p_H*0.97;} else {$p_H=$p_H*0.03/19;} 

					if ($m[3] eq "G") {$p_H=$p_H*0.99;} else {$p_H=$p_H*0.01/19;} 

					if ($m[4] eq "G") {$p_H=$p_H*0.78;} elsif ($m[4] eq "T") {$p_H=$p_H*0.17;} elsif ($m[4] eq "S") {$p_H=$p_H*0.024;} else {$p_H=$p_H*0.026/17;} 
	
					if ($m[5] eq "G") {$p_H=$p_H*0.99;} else {$p_H=$p_H*0.01/19;} 

					if ($m[6] eq "T") {$p_H=$p_H*0.99} else {$p_H=$p_H*0.01/19;} 

					if ($p_H>$max_P_H2) {$max_P_H2=$p_H; $best_segment_H2=$segment; $loc1_H2=$j; $loc2_H2=$k; } 
				}
			}
		}
	
		if (($max_P_H1*$max_P_H2>1e-16) && ($loc1_H1==$loc1_H2) && ($loc2_H2>$loc2_H1) && ($loc2_H2-$loc2_H1-6 > 0)) {
				$found=1;
				$cdr3temp=substr($AA[$loc1_H1], $loc2_H1+7, $loc2_H2-$loc2_H1-5);
				$aatemp=$AA[$loc1_H1];
				$cdnatemp=$cdna[$loc1_H1];
				$pos=$read-1;
				$cdr3len=length($cdr3temp);
				$cdr3cdna=substr($cdnatemp,($loc2_H1+7)*3,($loc2_H2-$loc2_H1-5)*3);
				print OutTemp1 "L\t$header,$cdr3temp\t$cdr3cdna\t$cdr3len\t$cdnatemp\t$aatemp\tR1L\t$iso\n";
				# $seqs{$header}[1]=$cdr3temp;
				# $seqs{$header}[$pos+2]=$cdr3temp;
				# $seqs{$header}[$pos+3]=$cdnatemp;
				# $seqs{$header}[$pos+4]=$aatemp;
		}
	}	
}
untie @InFile1;
close OutTemp1;

tie @InFile2,'Tie::File',$infile2.'-R2reads.fna';
open (OutTemp2,">$infile2".'-temp2.txt') or die $!;
# Process the 2nd file
$n=0;
for ($i=0;$i<@InFile2;$i++) {
	$line=$InFile2[$i];
	if ($line =~ m/^\>(.*)\s(.{1}):N:.*\<(.*)/){
		$header=$1;
		$read=$2;
		$iso=$3;
		$n++;
#		print "Processing R$read ...$n\n";
		next;
	}
	$current=$line;
	$current_r=reverse($current);
	$current_r=~tr/ATGC/TACG/;
	$len=length($current);
	@AA=();@cdna=();
	$cdna[1] = substr ($current, 0); for ($j=0; $j<$len; $j=$j+3) {$triplet = substr ($current, $j, 3); if ($triplet =~ "N") {$aa="X";} else {$aa=$base2aa{$triplet};} $AA[1]=$AA[1].$aa;}
    $cdna[2] = substr ($current, 1); for ($j=1; $j<$len; $j=$j+3) {$triplet = substr ($current, $j, 3); if ($triplet =~ "N") {$aa="X";} else {$aa=$base2aa{$triplet};} $AA[2]=$AA[2].$aa;}
    $cdna[3] = substr ($current, 2); for ($j=2; $j<$len; $j=$j+3) {$triplet = substr ($current, $j, 3); if ($triplet =~ "N") {$aa="X";} else {$aa=$base2aa{$triplet};} $AA[3]=$AA[3].$aa;}
    $cdna[4] = substr ($current_r, 0); for ($j=0; $j<$len; $j=$j+3) {$triplet = substr ($current_r, $j, 3); if ($triplet =~ "N") {$aa="X";} else {$aa=$base2aa{$triplet};} $AA[4]=$AA[4].$aa;}
    $cdna[5] = substr ($current_r, 1); for ($j=1; $j<$len; $j=$j+3) {$triplet = substr ($current_r, $j, 3); if ($triplet =~ "N") {$aa="X";} else {$aa=$base2aa{$triplet};} $AA[5]=$AA[5].$aa;}
    $cdna[6] = substr ($current_r, 2); for ($j=2; $j<$len; $j=$j+3) {$triplet = substr ($current_r, $j, 3); if ($triplet =~ "N") {$aa="X";} else {$aa=$base2aa{$triplet};} $AA[6]=$AA[6].$aa;}
	$max_P_H1=1e-99;     $max_P_H2=1e-99;
	$found=0;
	
	# Motif search #1
	for ($j=1; $j<=6; $j++) { 
		$current_AA=$AA[$j]; 
		$len_AA=length($current_AA);
		if ($len_AA<10) {$best_segment_H1="BLANK";$best_segment_H2="BLANK";}
		else {
			# Left Flank
			for ($k=0; $k<($len_AA-8); $k++) {
				$segment = substr ($current_AA, $k, 8); $p_H=1e-99;
	
				for ($kk=0; $kk<8; $kk++) {$m[$kk+1] = substr ($segment,$kk,1);}
 		
				if ($m[1] eq "D") {$p_H=0.99;} else {$p_H=0.01/19;} 
	
				if ($m[2] eq "T") {$p_H=$p_H*0.96;} elsif ($m[2] eq "S") {$p_H=$p_H*0.015;} elsif (($m[2] eq "A") or ($m[2] eq "M")) {$p_H=$p_H*0.0075;} else {$p_H=$p_H*0.01/16;} 

				if ($m[3] eq "A") {$p_H=$p_H*0.97;} elsif ($m[3] eq "G") {$p_H=$p_H*0.02;} else {$p_H=$p_H*0.01/18;} 

				if ($m[5] eq "Y") {$p_H=$p_H*0.99;} elsif ($m[5] eq "F") {$p_H=$p_H*0.005;} else {$p_H=$p_H*0.005/18;} 

				if ($m[6] eq "Y") {$p_H=$p_H*0.91;} elsif ($m[6] eq "F") {$p_H=$p_H*0.06;} elsif ($m[6] eq "S") {$p_H=$p_H*0.01;} elsif ($m[6] eq "H") {$p_H=$p_H*0.005;} else {$p_H=$p_H*0.015/16;} 
		
				if ($m[7] eq "C") {$p_H=$p_H*0.99;} else {$p_H=$p_H*0.01/19;} 

				if ($m[8] eq "A") {$p_H=$p_H*0.86;} elsif ($m[8] eq "T") {$p_H=$p_H*0.05;} elsif ($m[8] eq "V") {$p_H=$p_H*0.04;} elsif (($m[8] eq "G") or ($m[8] eq "I")) {$p_H=$p_H*0.01;} elsif (($m[8] eq "S") or ($m[8] eq "L")) {$p_H=$p_H*0.005;} else {$p_H=$p_H*0.02/13;} 
 
				if ($p_H>$max_P_H1) {$max_P_H1=$p_H; $loc1_H1=$j; $loc2_H1=$k;} 
			}
			# Right flank
			for ($k=0; $k<($len_AA-6); $k++) {
				$segment = substr ($current_AA, $k, 6); $p_H=1e-99; 
	
				for ($kk=0; $kk<6; $kk++) {$m[$kk+1] = substr ($segment,$kk,1);}
 		
				if ($m[1] eq "Y") {$p_H=0.40;} elsif ($m[1] eq "V") {$p_H=0.27;} elsif (($m[1] eq "I") or ($m[1] eq "L") or ($m[1] eq "P") or ($m[1] eq "S")) {$p_H=0.06;} elsif (($m[1] eq "F") or ($m[1] eq "H")) {$p_H=0.02;} else {$p_H=0.05/12;}
 
				if ($m[2] eq "W") {$p_H=$p_H*0.99;} else {$p_H=$p_H*0.01/19;} 

				if ($m[3] eq "G") {$p_H=$p_H*0.98;} else {$p_H=$p_H*0.02/19;} 

				if ($m[4] eq "Q") {$p_H=$p_H*0.83;} elsif ($m[4] eq "K") {$p_H=$p_H*0.07;} elsif ($m[4] eq "R") {$p_H=$p_H*0.05;} elsif ($m[4] eq "P") {$p_H=$p_H*0.02;} else {$p_H=$p_H*0.03/16;} 
	
				if ($m[5] eq "G") {$p_H=$p_H*0.98;} else {$p_H=$p_H*0.02/19;} 

				if ($m[6] eq "T") {$p_H=$p_H*0.93} elsif ($m[6] eq "S") {$p_H=$p_H*0.016;} elsif ($m[6] eq "I") {$p_H=$p_H*0.01;} elsif ($m[6] eq "G") {$p_H=$p_H*0.008;} elsif (($m[6] eq "P") or ($m[6] eq "N")) {$p_H=$p_H*0.005;} else {$p_H=$p_H*0.026/14} 

				if ($p_H>$max_P_H2) {$max_P_H2=$p_H; $loc1_H2=$j; $loc2_H2=$k;} 
			}
		}
	}
	
	if (($max_P_H1*$max_P_H2>1e-10) && ($loc1_H1==$loc1_H2) && ($loc2_H2>$loc2_H1) && ($loc2_H2-$loc2_H1-6 > 0)) {
			$found=1;
			$cdr3temp=substr($AA[$loc1_H1],$loc2_H1+6,$loc2_H2-$loc2_H1-4);
			$aatemp=$AA[$loc1_H1];
			$cdnatemp=$cdna[$loc1_H1];
			$pos=$read-1;
			$cdr3len=length($cdr3temp);
			$cdr3cdna=substr($cdnatemp,($loc2_H1+6)*3,($loc2_H2-$loc2_H1-4)*3);
			print OutTemp2 "H\t$header,$cdr3temp\t$cdr3cdna\t$cdr3len\t$cdnatemp\t$aatemp\tR2H\t$iso\n";
			# $seqs{$header}[0]=$cdr3temp;
			# $seqs{$header}[$pos+4]=$cdr3temp;
			# $seqs{$header}[$pos+5]=$cdnatemp;
			# $seqs{$header}[$pos+6]=$aatemp;
	}
	
	# Motif search #2
	if ($found==0) {
		$max_P_H1=1e-99;     $max_P_H2=1e-99; 
		for ($j=1; $j<=6; $j++) { 
			$current_AA=$AA[$j]; 
			$len_AA=length($current_AA);
			if ($len_AA<10) {$best_segment_H1="BLANK";$best_segment_H2="BLANK";}
			else {
				# Left Flank
				for ($k=0; $k<($len_AA-10); $k++) {
					$segment = substr ($current_AA, $k, 10); $p_H=1e-99;
	
					for ($kk=0; $kk<10; $kk++) {$m[$kk+1] = substr ($segment,$kk,1);}
 		
					if ($m[1] eq "E") {$p_H=0.91;} elsif ($m[1] eq "D") {$p_H=0.07;} elsif ($m[1] eq "A" || $m[1] eq "G") {$p_H=0.008/2;} else {$p_H=0.012/16;}
	
					if ($m[2] eq "D") {$p_H=$p_H*0.99;} else {$p_H=$p_H*0.01/19;} 

					if ($m[6] eq "Y") {$p_H=$p_H*0.99;} else {$p_H=$p_H*0.01/19;} 

					if ($m[7] eq "Y") {$p_H=$p_H*0.93;} elsif ($m[7] eq "F") {$p_H=$p_H*0.05;} else {$p_H=$p_H*0.02/18;} 

					if ($m[8] eq "C") {$p_H=$p_H*0.995;} else {$p_H=$p_H*0.005/19;} 
		
					if ($m[9] eq "Q") {$p_H=$p_H*0.83;} elsif ($m[9] eq "M") {$p_H=$p_H*0.10;} elsif ($m[9] eq "L" || $m[9] eq "H") {$p_H=$p_H*0.028/2;} else {$p_H=$p_H*0.042/16;}

					if ($m[10] eq "Q") {$p_H=$p_H*0.87;} elsif ($m[10] eq "H") {$p_H=$p_H*0.066;} elsif ($m[10] eq "K") {$p_H=$p_H*0.015;} else {$p_H=$p_H*0.049/17;} 
			 
					if ($p_H>$max_P_H1) {$max_P_H1=$p_H; $loc1_H1=$j; $loc2_H1=$k;} 
				}
				# Right Flank	
				for ($k=0; $k<($len_AA-6); $k++) {
					$segment = substr ($current_AA, $k, 6); $p_H=1e-99; 
	
					for ($kk=0; $kk<6; $kk++) {$m[$kk+1] = substr ($segment,$kk,1);}
 		
					if ($m[1] eq "T") {$p_H=0.88;} elsif ($m[1] eq "A") {$p_H=0.05;} elsif ($m[1] eq "S") {$p_H=0.034;} else {$p_H=0.036/17;}
 
					if ($m[2] eq "F") {$p_H=$p_H*0.99;} else {$p_H=$p_H*0.01/19;} 

					if ($m[3] eq "G") {$p_H=$p_H*0.99;} else {$p_H=$p_H*0.01/19;} 

					if ($m[4] eq "Q") {$p_H=$p_H*0.54;} elsif ($m[4] eq "G") {$p_H=$p_H*0.26;} elsif ($m[4] eq "P") {$p_H=$p_H*0.18;} else {$p_H=$p_H*0.02/17;} 
	
					if ($m[5] eq "G") {$p_H=$p_H*0.997;} else {$p_H=$p_H*0.003/19;} 

					if ($m[6] eq "T") {$p_H=$p_H*0.99} else {$p_H=$p_H*0.01/19;} 

					if ($p_H>$max_P_H2) {$max_P_H2=$p_H; $best_segment_H2=$segment; $loc1_H2=$j; $loc2_H2=$k; }
				}
			}
		}
	
		if (($max_P_H1*$max_P_H2>1e-16) && ($loc1_H1==$loc1_H2) && ($loc2_H2>$loc2_H1) && ($loc2_H2-$loc2_H1-6 > 0)) {
				$found=1;
				$cdr3temp=substr($AA[$loc1_H1], $loc2_H1+7, $loc2_H2-$loc2_H1-5);
				$aatemp=$AA[$loc1_H1];
				$cdnatemp=$cdna[$loc1_H1];
				$pos=$read-1;
				$cdr3len=length($cdr3temp);
				$cdr3cdna=substr($cdnatemp,($loc2_H1+7)*3,($loc2_H2-$loc2_H1-5)*3);
				print OutTemp2 "L\t$header,$cdr3temp\t$cdr3cdna\t$cdr3len\t$cdnatemp\t$aatemp\tR2L\t$iso\n";
				# $seqs{$header}[1]=$cdr3temp;
				# $seqs{$header}[$pos+4]=$cdr3temp;
				# $seqs{$header}[$pos+5]=$cdnatemp;
				# $seqs{$header}[$pos+6]=$aatemp;
		}
	}	
	
	# Motif search #3
	if ($found==0) {
		$max_P_H1=1e-99;     $max_P_H2=1e-99; 
		for ($j=1; $j<=6; $j++) { 
			$current_AA=$AA[$j]; 
			$len_AA=length($current_AA);
			if ($len_AA<10) {$best_segment_H1="BLANK";$best_segment_H2="BLANK";}
			else {
				# Left Flank
				for ($k=0; $k<($len_AA-8); $k++) {
					$segment = substr ($current_AA, $k, 8); $p_H=1e-99;
	
					for ($kk=0; $kk<8; $kk++) {$m[$kk+1] = substr ($segment,$kk,1);}
 		
					if ($m[1] eq "E") {$p_H=0.67;} elsif ($m[1] eq "G") {$p_H=0.17;} elsif ($m[1] eq "D") {$p_H=0.077;} elsif ($m[1] eq "M") {$p_H=0.059;} else {$p_H=0.024/16;}
	
					if ($m[2] eq "D") {$p_H=$p_H*0.99;} else {$p_H=$p_H*0.01/19;} 

					if ($m[3] eq "E") {$p_H=$p_H*0.98;} else {$p_H=$p_H*0.02/19;} 

					if ($m[4] eq "A") {$p_H=$p_H*0.93;} elsif ($m[4] eq "S") {$p_H=$p_H*0.04;} elsif ($m[4] eq "G") {$p_H=$p_H*0.02;} else {$p_H=$p_H*0.01/17;} 

					if ($m[5] eq "D") {$p_H=$p_H*0.92;} elsif ($m[5] eq "E") {$p_H=$p_H*0.038;} elsif ($m[5] eq "V") {$p_H=$p_H*0.012;} else {$p_H=$p_H*0.03/17;} 
		
					if ($m[6] eq "Y") {$p_H=$p_H*0.99;} elsif ($m[6] eq "F") {$p_H=$p_H*0.0095;} else {$p_H=$p_H*5e-4/18;}

					if ($m[7] eq "Y") {$p_H=$p_H*0.92;} elsif ($m[7] eq "F") {$p_H=$p_H*0.051;} elsif ($m[7] eq "H") {$p_H=$p_H*0.022;} else {$p_H=$p_H*7e-3/17;} 
			 
					if ($m[8] eq "C") {$p_H=$p_H*0.985;} elsif ($m[8] eq "R") {$p_H=$p_H*0.0079;} else {$p_H=$p_H*7.1e-3/18;}
			 
					if ($p_H>$max_P_H1) {$max_P_H1=$p_H; $best_segment_H1=$segment; $loc1_H1=$j; $loc2_H1=$k; }  
				}
				# Right Flank
				for ($k=0; $k<($len_AA-6); $k++) {
					$segment = substr ($current_AA, $k, 6); $p_H=1e-99; 
	
					for ($kk=0; $kk<6; $kk++) {$m[$kk+1] = substr ($segment,$kk,1);}
 		
					if ($m[1] eq "V") {$p_H=0.79;} elsif ($m[1] eq "I") {$p_H=0.076;} elsif ($m[1] eq "L") {$p_H=0.056;} elsif ($m[1] eq "A" || $m[1] eq "G" || $m[1] eq "M") {$p_H=0.014/3;} else {$p_H=0.064/14;}
 
					if ($m[2] eq "F") {$p_H=$p_H*0.97;} else {$p_H=$p_H*0.03/19;} 

					if ($m[3] eq "G") {$p_H=$p_H*0.99;} else {$p_H=$p_H*0.01/19;} 

					if ($m[4] eq "G") {$p_H=$p_H*0.78;} elsif ($m[4] eq "T") {$p_H=$p_H*0.17;} elsif ($m[4] eq "S") {$p_H=$p_H*0.024;} else {$p_H=$p_H*0.026/17;} 
	
					if ($m[5] eq "G") {$p_H=$p_H*0.99;} else {$p_H=$p_H*0.01/19;} 

					if ($m[6] eq "T") {$p_H=$p_H*0.99} else {$p_H=$p_H*0.01/19;} 

					if ($p_H>$max_P_H2) {$max_P_H2=$p_H; $best_segment_H2=$segment; $loc1_H2=$j; $loc2_H2=$k; } 
				}
			}
		}
	
		if (($max_P_H1*$max_P_H2>1e-16) && ($loc1_H1==$loc1_H2) && ($loc2_H2>$loc2_H1) && ($loc2_H2-$loc2_H1-6 > 0)) {
				$found=1;
				$cdr3temp=substr($AA[$loc1_H1], $loc2_H1+7, $loc2_H2-$loc2_H1-5);
				$aatemp=$AA[$loc1_H1];
				$cdnatemp=$cdna[$loc1_H1];
				$pos=$read-1;
				$cdr3len=length($cdr3temp);
				$cdr3cdna=substr($cdnatemp,($loc2_H1+7)*3,($loc2_H2-$loc2_H1-5)*3);
				print OutTemp2 "L\t$header,$cdr3temp\t$cdr3cdna\t$cdr3len\t$cdnatemp\t$aatemp\tR2L\t$iso\n";
				# $seqs{$header}[1]=$cdr3temp;
				# $seqs{$header}[$pos+4]=$cdr3temp;
				# $seqs{$header}[$pos+5]=$cdnatemp;
				# $seqs{$header}[$pos+6]=$aatemp;
		}
	}	
}
untie @InFile2;
close OutTemp2;

%seqs=[];
$tempname='./'.$infile1.'-temp1.txt';
tie @InFile,'Tie::File',$tempname;
for ($i=0;$i<@InFile;$i++) {
	$line=$InFile[$i];
	chomp($line);
	@tmp=split(/,/,$line);
	@label=split(/\t/,$tmp[0]);
	if ($label[0] eq 'H') {
		$seqs{$label[1]}[0]=$tmp[1];
	} else {
		$seqs{$label[1]}[1]=$tmp[1];
	}
}
untie @InFile;

$tempname='./'.$infile2.'-temp2.txt';
tie @InFile,'Tie::File',$tempname;
for ($i=0;$i<@InFile;$i++) {
	$line=$InFile[$i];
	chomp($line);
	@tmp=split(/,/,$line);
	@label=split(/\t/,$tmp[0]);
	if ($label[0] eq 'H') {
		$seqs{$label[1]}[0]=$tmp[1];
	} else {
		$seqs{$label[1]}[1]=$tmp[1];
	}
}
untie @InFile;

open (OutFile,">$infile1".'_paired_igblast_isotype.txt') or die $!;
print OutFile "Header\tCDRH3\tCDRL3\tH3_nt\tL3_nt\tCDRH3_junc_len\tCDRL3_junc_len\tVDJ_nt\tVDJ_aa\tVJ_nt\tVJ_aa\tHeavyRead\tLightRead\tR1Vgene\tR2Vgene\tR1Dgene\tR2Dgene\tR1Jgene\tR2Jgene\tHIsotype\tLIsotype\n";
map {
	$header=$_;
	unless (defined($R1vgene{$header})) {$R1vgene{$header}='-';}
	unless (defined($R1dgene{$header})) {$R1dgene{$header}='-';}
	unless (defined($R1jgene{$header})) {$R1jgene{$header}='-';}
	unless (defined($R2vgene{$header})) {$R2vgene{$header}='-';}
	unless (defined($R2dgene{$header})) {$R2dgene{$header}='-';}
	unless (defined($R2jgene{$header})) {$R2jgene{$header}='-';}
	if (defined($seqs{$header}[0]) && defined($seqs{$header}[1])) {
		@tmpH=split(/\t/,$seqs{$header}[0]);
		@tmpL=split(/\t/,$seqs{$header}[1]);
	print OutFile "$header\t$tmpH[0]\t$tmpL[0]\t$tmpH[1]\t$tmpL[1]\t$tmpH[2]\t$tmpL[2]\t$tmpH[3]\t$tmpH[4]\t$tmpL[3]\t$tmpL[4]\t$tmpH[5]\t$tmpL[5]\t$R1vgene{$header}\t$R2vgene{$header}\t$R1dgene{$header}\t$R2dgene{$header}\t$R1jgene{$header}\t$R2jgene{$header}\t$tmpH[6]\t$tmpL[6]\n";
	} 
	elsif (defined($seqs{$header}[0])) {
		@tmpH=split(/\t/,$seqs{$header}[0]);
		print OutFile "$header\t$tmpH[0]\t-\t$tmpH[1]\t-\t$tmpH[2]\t-\t$tmpH[3]\t$tmpH[4]\t-\t-\t$tmpH[5]\t-\t$R1vgene{$header}\t$R2vgene{$header}\t$R1dgene{$header}\t$R2dgene{$header}\t$R1jgene{$header}\t$R2jgene{$header}\t$tmpH[6]\t$tmpL[6]\n";
	} elsif (defined($seqs{$header}[1])) {
		@tmpL=split(/\t/,$seqs{$header}[1]);
		print OutFile "$header\t-\t$tmpL[0]\t-\t$tmpL[1]\t-\t$tmpL[2]\t-\t-\t$tmpL[3]\t$tmpL[4]\t-\t$tmpL[5]\t$R1vgene{$header}\t$R2vgene{$header}\t$R1dgene{$header}\t$R2dgene{$header}\t$R1jgene{$header}\t$R2jgene{$header}\t$tmpH[6]\t$tmpL[6]\n";
	}
} keys (%seqs);
close OutFile;

#unlink "$infile1".'-temp1.txt';
#unlink "$infile2".'-temp2.txt';
#unlink "$infile1".'-R1reads.fna';
#unlink "$infile2".'-R2reads.fna';