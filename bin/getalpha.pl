my %nonzeroct;
my %logct;
my %data;
my $metafile = $ARGV[0];
my $countdir = $ARGV[1];
my $outfile = $ARGV[2];
open FP1, "$metafile";
while(<FP1>) {
    chomp();
    my @b = split("\t");
    open FP, "$countdir\/count_$b[0]\.out";
    while(<FP>) {
	chomp();
	my @a = split("\t");
	$data{$b[0]}{$a[0]} = $a[1];
    }
    close FP;
}
close FP1;

foreach my $ccc (keys %data) {
    foreach my $gene (keys %{$data{$ccc}}) {
	$logct{$ccc} = $logct{$ccc} + log($data{$ccc}{$gene});
	$nonzeroct{$ccc}++;
    }
}

open OUT, ">$outfile";
open FP, "$metafile";
while(<FP>) {
    chomp();
    my @a = split("\t");
    if($nonzeroct{$a[0]} > 0) {
	my $aaa = $logct{$a[0]} / $nonzeroct{$a[0]};
	print OUT "$a[0]\t$aaa\t0\t0\t0\n";
    } else {
	print OUT "NA\n";
    }
}
close FP;
close OUT;
