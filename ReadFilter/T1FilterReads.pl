#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw(min);
use lib '/home/chuan/ChuanPerlLib';
use ChuanAnalysis;
use ChuanIO;

for (my $i = 7; $i < 18; $i++) {
	# get each fastq file in each folder
	my $str;
	if ($i < 10){
		$str = "0";
	}else{
		$str = "";
	}
	my $indir = '/home/das/Run_1724/Sample_763'.$str."$i";
	my $outfile1 = "763".$str."$i"."_R1.txt"; 
	my $outfileerr = "763".$str."$i"."_R1ERR.txt"; 
	my $diffseq1 = "763".$str."$i"."_DIFR1.txt"; 
	my $diffseq2 = "763".$str."$i"."_DIFR2.txt";

	my @files;
	opendir(DIR, $indir) or die $!;

	while (my $file = readdir(DIR)) {

		next unless (-f "$indir/$file");
		
		next unless ($file =~ m/\.fastq$/);

		push (@files, $file);
	}

	closedir(DIR);

	my $count = 0;

	unlink ($outfile1, $outfileerr, $diffseq1, $diffseq2);
	
	# output files for correct reads, incorrect reads, reads that differs between R1 and R2
	open(OUT1, ">".$outfile1)  or die "cannot open > output.txt: $!";
	open(OUTERR, ">".$outfileerr)  or die "cannot open > output.txt: $!";
	open(DIFFSeq1, ">".$diffseq1)  or die "cannot open > output.txt: $!";
	open(DIFFSeq2, ">".$diffseq2)  or die "cannot open > output.txt: $!";

	chdir $indir;
	
	# parse through each pairs of file
	foreach my $infile1 (@files) {

		if ($infile1 =~ /_R1_/) {

			open(INPUT1, "<".$infile1) or die "cannot open < input.txt: $!";
			my $infile2 = $infile1;
			$infile2 =~ s/_R1_/_R2_/g;
			open(INPUT2, "<".$infile2) or die "cannot open < input.txt: $!";

			while (<INPUT1>){
				<INPUT2>;
				if ($_ =~ /^\@D004/) { 
					chomp(my $seq1 = <INPUT1>); # read in next line
					chomp(my $seq2 = <INPUT2>); # read in next line
					if ($seq1 =~ /CCAAGTTG/ & $seq1 =~ /TTGATTAT/) {
						if ($seq2 =~ /ATAATCAA/ & $seq2 =~ /CAACTTGG/) {

							# get sequence1 and sequence2 by the flanking sequence
							$seq1 =~ s/.*CCAAGTTG(.+)TTGATTAT.*/$1/;
							$seq2 =~ s/.*ATAATCAA(.+)CAACTTGG.*/$1/;
							$seq2 = revcom($seq2);

							if ($seq1 eq $seq2) {
								# make sure they have the correct size of tRNA gene			
								if (length($seq1) == 72) {
									$count = $count + 1;
									print OUT1 $seq1."\n";
								}else{
									print OUTERR $seq1."\n";
								}
							}else{
								# when R1 and R2 disagree, disgard
								print DIFFSeq1 $seq1."\n";
								print DIFFSeq2 $seq2."\n";
							}
						}
					}
						
				}
			}
			close INPUT1;
			close INPUT2;

		}
	
	}

	close OUT1;
	close OUTERR;
	close DIFFSeq1;
	close DIFFSeq2;
}

exit; 