#!/usr/bin/perl
use strict;
use warnings;

# A less verbose version of this code is T4CountsOfVariants.pl
# this is the actual one used...
# extract some sequence as example
chdir "/home/chuan/das/Run_1691/Outputs";

my $infile = "Variants_filtered.txt";

my $outfile0 = "Mut0.txt";
my $outfile1 = "Mut1.txt";
my $outfile2 = "Mut2.txt";
my $outfile3 = "Mut3.txt";
my $outfile4 = "Mut4.txt";
my $outfile5 = "Mut5.txt";
my $outfile6 = "Mut6.txt";
my $outfile7 = "Mut7.txt";
my $outfile8 = "Mut8.txt";
my $outfile9 = "Mut9.txt"; # 9 contains 9+

unlink ($outfile0, $outfile1, $outfile2, $outfile3, $outfile4, $outfile5, $outfile6, $outfile7, $outfile8, $outfile9);

open(INPUT, "<".$infile) or die "cannot open < input.txt: $!";

open(OUT0, ">".$outfile0)  or die "cannot open > output.txt: $!";
open(OUT1, ">".$outfile1)  or die "cannot open > output.txt: $!";
open(OUT2, ">".$outfile2)  or die "cannot open > output.txt: $!";
open(OUT3, ">".$outfile3)  or die "cannot open > output.txt: $!";
open(OUT4, ">".$outfile4)  or die "cannot open > output.txt: $!";
open(OUT5, ">".$outfile5)  or die "cannot open > output.txt: $!";
open(OUT6, ">".$outfile6)  or die "cannot open > output.txt: $!";
open(OUT7, ">".$outfile7)  or die "cannot open > output.txt: $!";
open(OUT8, ">".$outfile8)  or die "cannot open > output.txt: $!";
open(OUT9, ">".$outfile9)  or die "cannot open > output.txt: $!";

my $refseq = "GTTCCGTTGGCGTAATGGTAACGCGTCTCCCTCCTAAGGAGAAGACTGCGGGTTCGAGTCCCGTACGGAACG";

while (<INPUT>){
	chomp;
	my @seq = split '\t', $_;
	my $seq = $seq[0];

	my ($posref, $mutref) = compareString2($seq, $refseq);
	my @pos = @{$posref};
	my @mut = @{$mutref};
	shift @seq;

	if ($pos[0] == -1) {
			print OUT0 "@pos\t";
			print OUT0 "@mut\t";
			print OUT0 join "\t", @seq;
			print OUT0 "\n";
	}elsif(scalar @pos == 1){
			print OUT1 "@pos\t";
			print OUT1 "@mut\t";
			print OUT1 join "\t", @seq;
			print OUT1 "\n";
	}elsif(scalar @pos == 2){
			print OUT2 "@pos\t";
			print OUT2 "@mut\t";
			print OUT2 join "\t", @seq;
			print OUT2 "\n";
	}elsif(scalar @pos == 3){
			print OUT3 "@pos\t";
			print OUT3 "@mut\t";
			print OUT3 join "\t", @seq;
			print OUT3 "\n";
	}elsif(scalar @pos == 4){
			print OUT4 "@pos\t";
			print OUT4 "@mut\t";
			print OUT4 join "\t", @seq;
			print OUT4 "\n";
	}elsif(scalar @pos == 5){
			print OUT5 "@pos\t";
			print OUT5 "@mut\t";
			print OUT5 join "\t", @seq;
			print OUT5 "\n";
	}elsif(scalar @pos == 6){
			print OUT6 "@pos\t";
			print OUT6 "@mut\t";
			print OUT6 join "\t", @seq;
			print OUT6 "\n";
	}elsif(scalar @pos == 7){
			print OUT7 "@pos\t";
			print OUT7 "@mut\t";
			print OUT7 join "\t", @seq;
			print OUT7 "\n";
	}elsif(scalar @pos == 8){
			print OUT8 "@pos\t";
			print OUT8 "@mut\t";
			print OUT8 join "\t", @seq;
			print OUT8 "\n";
	}else{
			print OUT9 "@pos\t";
			print OUT9 "@mut\t";
			print OUT9 join "\t", @seq;
			print OUT9 "\n";
	} 

}

exit;

###########SUBFUNCTIONS#############
sub compareString2{
	my ($seq, $refseq) = @_;
	my @x = split '', $seq;
	my @y = split '', $refseq;
	my $char = 0;
	my @posoutput;
	my @mutoutput;

	if ($seq eq $refseq ) {
		return ([-1], ['V']);
	}

	my $string = join '',
				 map { $x[$_] eq $y[$_] ? $y[$_] :$char}
				 0 .. $#y;

	# find every index
	my $offset = 0;
	my $result = index($string, $char, $offset);

	while ($result != -1) {
		push @posoutput, $result;

		my $mutant = substr $seq, $result, 1;
		push @mutoutput, $mutant;
		$offset = $result + 1;
		$result = index($string, $char, $offset);
		# record the ATGC status

	}
	
	return (\@posoutput, \@mutoutput);
}