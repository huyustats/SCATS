my $cdt0 = $ARGV[0];
my $cdt1 = $ARGV[1];
my $tmpdir = $ARGV[2];
my $metafile = $ARGV[3];
my $gpinfofile = $ARGV[4];

my $gpp = $cdt0."_".$cdt1;

my @condition;
my %qout;
open FP, "$tmpdir\/abkt/abkt_umi";
while(<FP>) {
    chomp();
    my @a = split("\t");
    $qout{$a[0]} = 1;
}
close FP;

my %count;
my %quality1;
my %quality2;
open FP, "$metafile";
while(<FP>) {
    chomp();
    my @a = split("\t");
    if($qout{$a[0]} == 1) {
	#print "$a[0]\n";
	open FP1, "$tmpdir\/count_script\/count_$a[0]\.out";
	while(<FP1>) {
	    chomp();
	    my @b = split("\t");
	    #$count{$b[0].":".$b[2]}{$a[0]}{$b[3]} = $b[4] if $qgene1{$b[0]} == 1;
	    $count{$b[0].":".$b[2]}{$a[0]}{$b[3]} = $b[4];
	    if($a[$#a] =~ $cdt0 || $a[$#a] =~ $cdt1) {
		$quality1{$b[0].":".$b[2]}{$b[3]} = $quality1{$b[0].":".$b[2]}{$b[3]} + $b[4];
		$quality2{$b[0].":".$b[2]}{$b[3]}++ if $b[4] > 0;
	    }
	}
	close FP1;
    }
}
close FP;

my %abkt;
open FP, "$tmpdir\/abkt/abkt_umi";
while(<FP>) {
    chomp();
    my @a = split("\t");
    $abkt{$a[0]} = "$a[1]\t$a[2]\t$a[3]\t$a[4]";
}
close FP;

my @abktfile;
my @cell;
open FP, "$metafile";
while(<FP>) {
    chomp();
    my @a = split("\t");
    if($qout{$a[0]} == 1) {
	if($a[1] eq $cdt0) {
	    @condition = (@condition, 0);
	    @cell = (@cell, $a[0]);
	    @abktfile = (@abktfile, $abkt{$a[0]});
	}
	if($a[1] eq $cdt1) {
            @condition = (@condition, 1);
            @cell = (@cell, $a[0]);
	    @abktfile = (@abktfile, $abkt{$a[0]});
        }
    }
}
close FP;

open OUT1, ">$tmpdir\/das_script/data/condition_$gpp";
foreach my $i (0..$#condition) {
    print OUT1 "$condition[$i]\n";
}
close OUT1;

open OUT3, ">$tmpdir\/das_script/data/abktfile_$gpp";
foreach my $i (0..$#abktfile) {
    print OUT3 "$abktfile[$i]\n";
}
close OUT3;

my %prop;
my %qgroup;
open FP, "$gpinfofile";
while(<FP>) {
    chomp();
    my @a = split("\t");
    my $gp = $a[0].":".$a[1];
    $prop{$gp}{$a[3]} = log($a[4]);
    $qgroup{$gp} = 1;
}
close FP;

my %qgene;
my %mean;
my %bursting;
open FP, "$tmpdir\/gene_script/geneleveltheta_umi";
while(<FP>) {
    chomp();
    next if /theta_rd/;
    my @a = split("\t");
    next if ($a[0] ne $cdt0 && $a[0] ne $cdt1);
    $qgene{$a[0]}{$a[1]} = 1;
    $mean{$a[0]}{$a[1]} = $a[2];
}
close FP;

open FP, "$tmpdir\/gene_script/genelevelbursting_umi";
while(<FP>) {
    chomp();
    my @a = split("\t");
    $bursting{$a[0]}{$a[1]} = $a[2];
}
close FP;


open OUT2, ">$tmpdir\/das_script/data/countdata_$gpp";
foreach my $gp (keys %count) {
    my @a = split(":",$gp);
    next if !($qgene{$cdt0}{$a[0]} == 1);
    next if !($qgene{$cdt1}{$a[0]} == 1);
    next if !($qgroup{$gp} == 1);
    #print "yes\n";
    my $genect;
    my $plus;
    my $minus;
    my $tmpquality = 1;
    if($quality1{$gp}{"+"} > 50 || $quality1{$gp}{"-"} > 50) {
	$tmpquality = 0 if !($quality2{$gp}{"+"} > 15 || $quality2{$gp}{"-"} > 15);
    } else {
	$tmpquality = 0 if !($quality2{$gp}{"+"} > 10 && $quality2{$gp}{"-"} > 10);
    }
    #next if $tmpquality == 0;
    foreach my $i (0..$#cell) {
	$genect = $genect."100," if $i < $#cell;
	$genect = $genect."100" if $i == $#cell;
	if(!($count{$gp}{$cell[$i]}{"+"}>0)) {
	    $count{$gp}{$cell[$i]}{"+"} = 0;
	}
	if(!($count{$gp}{$cell[$i]}{"-"}>0)) {
            $count{$gp}{$cell[$i]}{"-"} = 0;
        }
	$plus = $plus.$count{$gp}{$cell[$i]}{"+"}."," if $i < $#cell;
	$plus = $plus.$count{$gp}{$cell[$i]}{"+"} if $i == $#cell;
	$minus = $minus.$count{$gp}{$cell[$i]}{"-"}."," if $i < $#cell;
        $minus = $minus.$count{$gp}{$cell[$i]}{"-"} if $i == $#cell;
    }
    #print OUT2 "$gp\t$genect\t$plus\t$minus\t$mean{$cdt0}{$a[0]},1,$mean{$cdt1}{$a[0]},1,$prop{$gp}{\"+\"},$prop{$gp}{\"-\"}\t+\\-\n";
    my $tmpgp = $cdt0."_".$cdt1;    
    print OUT2 "$gp\t$genect\t$plus\t$minus\t$mean{$cdt0}{$a[0]},1,$mean{$cdt1}{$a[0]},1,$bursting{$tmpgp}{$a[0]}\t+\\-\n";
}
close OUT2;
