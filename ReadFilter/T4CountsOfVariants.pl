#!/usr/bin/perl
use strict;
use warnings;

chdir "/home/chuan/das/Run_1691/Outputs";
my $len;

my $str = "942";
for (my $i = 22; $i < 32; $i++) {
	my $infile = $str.$i."_Align.txt";

	my $outfile = $str.$i."_Variants.txt";

	unlink $outfile;

	open(INPUT, "<".$infile) or die "cannot open < input.txt: $!";

	open(OUT, ">".$outfile)  or die "cannot open > output.txt: $!";

	my $refseq = "GTTCCGTTGGCGTAATGGTAACGCGTCTCCCTCCTAAGGAGAAGACTGCGGGTTCGAGTCCCGTACGGAACG";

	while (<INPUT>){
		chomp;
		my $line = $_;
		my @seq = split '\t', $line;
		my $seq = $seq[0];

		my ($posref, $mutref) = compareString2($seq, $refseq);
		my @pos = @{$posref};
		my @mut = @{$mutref};
		shift @seq;

		if ($pos[0] == -1) {
			$len = 0;
		}else{
			$len = scalar @pos;
		}
		print OUT $len."\t";
		print OUT "@pos\t";
		print OUT "@mut\t";
		print OUT $line."\t";
		print OUT "\n";

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