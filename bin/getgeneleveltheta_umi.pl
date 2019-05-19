my %theta;
my %qgene;
my %bursting;
my $ct = 0;
my $comparedir = $ARGV[0];
my $genedir = $comparedir."/gene_script";
my $datadir = $genedir."/data";

open FP, "$comparedir\/comparegroup";
while(<FP>) {
    chomp();
    $ct++;
    my @a = split("\t");
    my $gp = "$a[0]\_$a[1]";
    open FP1, "$datadir\/outgene_$gp";
    while(<FP1>) {
	chomp();
	my @b = split("\t");
	if($b[1] eq "True" && $b[5] ne "nan" && $b[6] ne "nan") {
	    #$qgene{$a[0]}{$b[0]} = 1;
	    $theta{$a[0]}{$b[0]} = $b[5];
	    #$qgene{$a[1]}{$b[0]} = 1;
            $theta{$a[1]}{$b[0]} = $b[6];
	    $bursting{$gp}{$b[0]} = $b[8];

	    #print "$b[8]\n";
	}
    }
    close FP1;
}
close FP;

open OUT, ">$genedir\/geneleveltheta_umi";
open FP, "$comparedir\/celltypes";
while(<FP>) {
    chomp();
    foreach my $g (keys %{$theta{$_}}) {
	#print "$bursting{$_}{$g}\n";
	print OUT "$_\t$g\t$theta{$_}{$g}\n";
    }
}
close FP;
close OUT;

open OUT, ">$genedir\/genelevelbursting_umi";
open FP, "$comparedir\/comparegroup";
while(<FP>) {
    chomp();
    my @a = split("\t");
    my $gp = "$a[0]\_$a[1]";
    foreach my $g (keys %{$bursting{$gp}}) {
	print OUT "$gp\t$g\t$bursting{$gp}{$g}\n";
    }
}
close FP;
close OUT;
