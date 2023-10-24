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
 $igblastpre=dirname(abs_path($0)).'/./igblastn -germline_db_V '.dirname(abs_path($0)).'/custom/custom_hu_re_gl_V -germline_db_D '.dirname(abs_path($0)).'/custom/custom_hu_re_gl_D -germline_db_J '.dirname(abs_path($0)).'/custom/custom_hu_re_gl_J -auxiliary_data '.dirname(abs_path($0)).'/custom_gl.aux -domain_system imgt -num_alignments_V 1 -num_alignments_D 1 -num_alignments_J 1 -outfmt 19 -query ';
$igblastoutput1=$infile1.'-R1reads-igblast.txt';
$igblastoutput2=$infile2.'-R2reads-igblast.txt';
 $igblastcmd1=$igblastpre.$infile1.'-R1reads.fna >'.$igblastoutput1;
 $igblastcmd2=$igblastpre.$infile2.'-R2reads.fna > '.$igblastoutput2;
#print ($igblastcmd1); 
system($igblastcmd1);
system($igblastcmd2);
$cmd1 = 'cat '. $igblastoutput1.' | awk -F"\t" \'{if($6=="T" && length($48)>0) print $1 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $14 "\t" $48  "\t" $50 "\t" $51 }\''.'|sed \'s/_F</\t/\''.' > '.$igblastoutput1.'-temp.txt';
system($cmd1);
$cmd2 = 'cat '. $igblastoutput2.' | awk -F"\t" \'{if($6=="T" && length($48)>0) print $1 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\t" $14 "\t" $48  "\t" $50 "\t" $51 }\''.'|sed \'s/_R</\t/\''.' > '.$igblastoutput2.'-temp.txt';
system($cmd2);
$cmd3='sort -k 1,1 '.$igblastoutput1.'-temp.txt'.' > '.$igblastoutput1.'-temp.txt_sorted';
system($cmd3);
$cmd4='sort -k 1,1 '.$igblastoutput2.'-temp.txt'.' > '.$igblastoutput2.'-temp.txt_sorted';
system($cmd4);

$cmd5= 'join -j 1 -o 1.1 1.9 2.8 1.8 2.7 1.10 2.9 1.6 1.7 2.5 2.6 1.3 2.3 1.4 1.5 2.4 1.2 2.2 '. $igblastoutput1.'-temp.txt_sorted'. ' '.$igblastoutput2.'-temp.txt_sorted'. '| awk \'{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" "R1H" "\t" "R2L" "\t" $12 "\t" $13 "\t" $14 "\t" "-" "\t" $15 "\t" $16 "\t" $17 "\t" $18}\' '.' > '.$infile1.'_paired_igblast_isotype.txt';
system($cmd5);



