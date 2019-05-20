my %results;
my %qgene;
my $ct = 0;
my $comparedir = $ARGV[0];
my $dasdir = $comparedir."/das_script";
my $datadir = $dasdir."/data";
my $infofile = $ARGV[1];
my $outfile = $ARGV[2];

my %event2exon;
open FP, "$infofile";
while(<FP>) {
    chomp();
    my @a = split("\t");
    my $tmp = $a[0].":".$a[1];
    $event2exon{$tmp} = $a[$#a] if $a[$#a] ne "NA";
}
close FP;


open OUT, ">$outfile";
print OUT "gp1\tgp2\tgene_name\tAS_exons\tPSI_gp1\tPSI_gp2\ttest_stat\tp_value\n";
open FP, "$comparedir\/comparegroup";
while(<FP>) {
    chomp();
    $ct++;
    my @a = split("\t");
    my $gp = "$a[0]\_$a[1]";
    open FP1, "$datadir\/out_$gp";
    while(<FP1>) {
	chomp();
	my @b = split("\t");
	if($b[1] eq "True" && $b[5] ne "nan" && $b[6] ne "nan") {
	    my @c = split(":", $b[0]);
	    print OUT "$a[0]\t$a[1]\t$c[0]\t$event2exon{$b[0]}\t$b[$#b-5]\t$b[$#b-4]\t$b[$#b-1]\t$b[$#b]\n";
	    
	    #$results{$gp}{$b[0]}{"pv"} = $b[$#b];
	    #$results{$gp}{$b[0]}{"stat"} = $b[$#b-1];
	    #$results{$gp}{$b[0]}{0} = $b[$#b-5];
	    #$results{$gp}{$b[0]}{1} = $b[$#b-4];
	}
    }
    close FP1;
}
close FP;
close OUT;
