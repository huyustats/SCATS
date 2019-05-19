my $cdt0 = $ARGV[0];
my $cdt1 = $ARGV[1];
my @condition;
my %qout;
my $abktfile = $ARGV[2];
my $metafile = $ARGV[3];
my $tmpdir = $ARGV[4];
my $tmpdirgene = $ARGV[5];
open FP, "$abktfile";
while(<FP>) {
    chomp();
    my @a = split("\t");
    $qout{$a[0]} = 1;
}
close FP;

open FP, "$metafile";
while(<FP>) {
    chomp();
    my @a = split("\t");
    if($qout{$a[0]} == 1) {
	#print "$a[0]\n";
	open FP1, "$tmpdir\/count_$a[0]\.out";
	while(<FP1>) {
	    chomp();
	    my @b = split("\t");
	    #$count{$b[0]}{$a[0]} = $b[1] if $qgene{$b[0]} == 1;
	    $count{$b[0]}{$a[0]} = $b[1];
	}
	close FP1;
    }
}
close FP;

my %abkt;
open FP, "$abktfile";
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

open OUT1, ">$tmpdirgene\/condition_$cdt0\_$cdt1";
foreach my $i (0..$#condition) {
    print OUT1 "$condition[$i]\n";
}
close OUT1;

open OUT2, ">$tmpdirgene\/countdata_$cdt0\_$cdt1";
foreach my $gene (keys %count) {
    print OUT2 "$gene\t";
    foreach my $i (0..$#cell) {
	if($i < $#cell) {
	    if($count{$gene}{$cell[$i]} > 0) {
		print OUT2 "$count{$gene}{$cell[$i]}\t";
	    } else {
		print OUT2 "0\t";
	    }
	}
	if($i == $#cell) {
	    if($count{$gene}{$cell[$i]} > 0) {
                print OUT2 "$count{$gene}{$cell[$i]}\n";
            } else {
                print OUT2 "0\n";
            }
	}
    }
}
close OUT2;

open OUT3, ">$tmpdirgene\/abktfile_$cdt0\_$cdt1";
foreach my $i (0..$#abktfile) {
    print OUT3 "$abktfile[$i]\n";
}
close OUT3;
