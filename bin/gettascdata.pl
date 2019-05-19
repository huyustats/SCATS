my %cdt;
my $ct=0;
my $cdt0 = $ARGV[0];
my $cdt1 = $ARGV[1];
my $tmpdir = $ARGV[2];
my $outdir = $ARGV[3];
open FP, "$tmpdir\/condition_$cdt0\_$cdt1";
while(<FP>) {
    chomp();
    $ct++;
    $cdt{$ct} = $_;
}
close FP;

my %qgene;
open FP, "$tmpdir\/countdata_$cdt0\_$cdt1";
while(<FP>) {
    chomp();
    my @a = split("\t");
    my %tmp;
    foreach my $i (1..$#a) {
	$tmp{$cdt{$i}} = $tmp{$cdt{$i}} + $a[$i];
    }
    if($tmp{1} >= 10 && $tmp{0} >= 10) {
	$qgene{$a[0]} = 1;
    }
}
close FP;


open OUT, ">$outdir\/tascdata_$cdt0\_$cdt1";
open FP, "$tmpdir\/countdata_$cdt0\_$cdt1";
while(<FP>) {
    chomp();
    if(/ERCC/) {
	next;
    }
    my @a = split("\t");
    if($qgene{$a[0]} == 1) {
	print OUT "$a[0]\t";
	foreach my $i (1..($#a-1)) {
	    print OUT "$a[$i],";
	}
	print OUT "$a[$#a]\t";
	print OUT "1\t1\t1\t1\n";
    }
}
close FP;
close OUT;
