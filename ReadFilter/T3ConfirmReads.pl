#!/usr/bin/perl
use strict;
use warnings;

chdir "/home/das/Run_1691";

# Align sequence with tRNA
my $infile1 = "Variants_All.txt";	
my $outfile = "Variants_Align.txt";
my $outerr = "Variants_discard.txt";

open(IN1, "<", $infile1)
	or die "cannot open < input.txt: $!";
	
unlink ($outfile, $outerr);	
open(OUT, ">", $outfile) or die "cannot open > output.txt: $!";
open(OUTERR, ">", $outerr) or die "cannot open > output.txt: $!";

my %hash;

####################
# output1: alignment sequences
####################

my $seq1 = "ATGCTGGAA"."GTTCCGTTGGCGTAATGGTAACGCGTCTCCCTCCTAAGGAGAAGACTGCGGGTTCGAGTCCCGTACGGAACG"."ATTGTGCA";

while (<IN1>){
	chomp;
	my $line = $_;
	my @seq = split '\t', $line;
	
	my $seq2 = "ATGCTGGAA".$seq[0]."ATTGTGCA";
	my ($ali1,$ali2)= &runClustal($seq1,$seq2);
	
	if ($ali1 =~ /-/ or $ali2 =~ /-/) {
		print OUTERR $line."\n";
		print OUTERR $ali1."\n";
		print OUTERR $ali2."\n";
		
	}else{
		# write hash output for correct sequences		
		print OUT $line."\n";	
	}

}





############## SUBFUNCTIONS ################
# subfunction for read alignment
=head1 use runClustal
	my $pep = "ATTTTTGG";
	my $pro_tmp = "ATGTTTTGG";
	my ($ali1T,$ali2T)= &runClustal($pep,$pro_tmp);
	print $ali1T;
	print "\n";
=cut

sub runClustal{
    my ($p1, $p2) = @_;
    my $rand = int(rand(9999999)); my $tmpf = "rand$rand\.fa";
    my $dndf = "rand$rand\.dnd";

    open(tmpClu, ">$tmpf");print tmpClu ">p1\n$p1\n>p2\n$p2\n";close(tmpClu);
    $rand = int(rand(9999999)); my $outf = "rand$rand\.aln";                                                                                          
#	print "/usr/bin/clustalw2 -infile=$tmpf -outfile=$outf\n";
	`/usr/bin/clustalw -infile=$tmpf -outfile=$outf`;
    if(! -s "$outf"){goto no_outf;}
    open(tmpClu, "<$outf");my @aln = <tmpClu>;close(tmpClu);
    my $ali1 = ''; my $ali2 = '';
    foreach my $aln (@aln){
        if($aln =~ /^(\S+)\s+(\S+)/){
            my $idT = $1; my $aliT = $2;
            if($idT eq 'p1'){$ali1 .= $aliT;}
            if($idT eq 'p2'){$ali2 .= $aliT;}
        }
    }
   `rm $tmpf $outf $dndf`;
 
    return ($ali1,$ali2);
  no_outf:;
}

