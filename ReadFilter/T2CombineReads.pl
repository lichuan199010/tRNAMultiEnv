#!/usr/bin/perl
# combine the readout from run 1724 and 1691
use strict;
use warnings;
use lib '/home/chuan/ChuanPerlLib';
use ChuanAnalysis;
use ChuanIO;

## read in the file for T0 (combining two lanes)
chdir "/home/das/Run_1691/Outputs";
my $str = "72585";

# T0_1 read in First sample of T0 in lane 1691
openInLC($str."_R1.txt");
my %hash;
while (<INPUT>) {
	chomp;
	if (exists $hash{$_}) {
		$hash{$_} = $hash{$_}+1;
	}else{
		$hash{$_} = 1;
	}
}
close INPUT;

# T0-2 read in second sample of T0 in lane 1691
$str = "72592";
openInLC($str."_R1.txt");
while (<INPUT>) {
	chomp;
	if (exists $hash{$_}) {
		$hash{$_} = $hash{$_}+1;
	}else{
		$hash{$_} = 1;
	}
}
close INPUT;

# T0_3 read in first sample of T0 in lane 1724
$str = "76307";
openInLC($str."_R1.txt");
while (<INPUT>) {
	chomp;
	if (exists $hash{$_}) {
		$hash{$_} = $hash{$_}+1;
	}else{
		$hash{$_} = 1;
	}
}
close INPUT;

# T0_4 read in first sample of T0 in lane 1724
$str = "76312";
openInLC($str."_R1.txt");
while (<INPUT>) {
	chomp;
	if (exists $hash{$_}) {
		$hash{$_} = $hash{$_}+1;
	}else{
		$hash{$_} = 1;
	}
}
close INPUT;
my %hash_copy = %hash;

# hash_copy stores the number for T0 combined
# hash store information from all samples


## read in from all samples (including re-reading T0 with separate numbers for each sample)
for (my $dircount = 85; $dircount < 93; $dircount ++) {
	my %hash2 = ();
	openInLC("725".$dircount."_R1.txt");
	
	while (<INPUT>) {
		chomp;
		if (exists $hash2{$_}) {
			$hash2{$_} = $hash2{$_}+1;
		}else{
			$hash2{$_} = 1;
		}
	}
	close INPUT;

	foreach my $key (keys %hash) {

		if (exists $hash2{$key}) {
			$hash{$key} = $hash{$key}."\t".$hash2{$key};
		}else{
			$hash{$key} = $hash{$key}."\t0";
		}
	}
}


for (my $dircount = 7; $dircount < 18; $dircount ++) {
	if ($dircount < 10){
		$str = "0";
	}else{
		$str = "";
	}
	my %hash2 = ();
	openInLC("763".$str.$dircount."_R1.txt");
	
	while (<INPUT>) {
		chomp;
		if (exists $hash2{$_}) {
			$hash2{$_} = $hash2{$_}+1;
		}else{
			$hash2{$_} = 1;
		}
	}
	close INPUT;

	foreach my $key (keys %hash) {

		if (exists $hash2{$key}) {
			$hash{$key} = $hash{$key}."\t".$hash2{$key};
		}else{
			$hash{$key} = $hash{$key}."\t0";
		}
	}
}



## print all outputs
my $outfile = "Variants_All.txt";
unlink $outfile;
openOutLC($outfile);

foreach my $key (keys %hash_copy){
	if ($hash_copy{$key} >= 100){
		print OUT $key."\t";
		print OUT $hash{$key}."\n";
	}
}

close OUT;
exit;

