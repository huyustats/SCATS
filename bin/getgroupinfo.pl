my $input = $ARGV[0];
#my $output = $ARGV[1];
my %genecount;
my $totalweight;
my $tmpweight;
my $tmpisoset;
my $tmpisoindex;
my $tmpisoindex1;
my @tmpiso;
my %isosetinf;
my %genelength;
my %grouplength;
my %isoindex;
my %isoindex_complement;
my %exonset;
open FP, "$input";
while(<FP>) {
    chomp();
    my @a = split("\t");
    my @b = split(",",$a[$#a]);
    $genecount{$a[0]}++;
    if($genecount{$a[0]} == 1) {
	@tmpiso = split(",",$a[$#a]);
	$totalweight = 0;
	foreach my $i (0..$#tmpiso) {
	    $totalweight = $totalweight + 2 ** $i;
	}
    } else {
	$genelength{$a[0]} = $genelength{$a[0]} + $a[4] - $a[3] + 1;
	my $tmpexon = $a[1].",".$a[3].",".$a[4];
	if($a[$#a] =~ /0,/) {
	    $tmpweight = 0;
	    $tmpisoset = "";
	    $tmpisoindex = "";
	    $tmpisoindex1 = "";
	    foreach my $i (0..$#b) {
		$tmpweight = $tmpweight + 2 ** $i if $b[$i] == 1;
		$tmpisoset = $tmpisoset.$tmpiso[$i]."," if $b[$i] == 1;
		$tmpisoindex = $tmpisoindex.$i."," if $b[$i] == 0;
		$tmpisoindex1 = $tmpisoindex1.$i."," if $b[$i] == 1;
	    }
	    if($tmpweight < $totalweight/2) {
		$isoindex{$a[0]}{$tmpweight}{1} = $tmpisoindex1;
		$isoindex_complement{$a[0]}{$tmpweight}{1} = $tmpisoindex;
		$isosetinf{$a[0]}{$tmpweight}{1} = $tmpisoset;
		$grouplength{$a[0]}{$tmpweight}{1} = $grouplength{$a[0]}{$tmpweight}{1} + $a[4] - $a[3] + 1;
		$exonset{$a[0]}{$tmpweight}{1} = $exonset{$a[0]}{$tmpweight}{1}.$tmpexon.";";
	    } else {
		$tmpweight = $totalweight - $tmpweight;
		$isoindex{$a[0]}{$tmpweight}{0} = $tmpisoindex1;
		$isoindex_complement{$a[0]}{$tmpweight}{0} = $tmpisoindex;
		$isosetinf{$a[0]}{$tmpweight}{0} = $tmpisoset;
		$grouplength{$a[0]}{$tmpweight}{0} = $grouplength{$a[0]}{$tmpweight}{0} + $a[4] - $a[3] + 1;
		$exonset{$a[0]}{$tmpweight}{0} = $exonset{$a[0]}{$tmpweight}{0}.$tmpexon.";";
	    }
	}
    }
}
close FP;


foreach my $gene (keys %isosetinf) {
    foreach my $weight (keys %{$isosetinf{$gene}}) {
	my $status;
	if($grouplength{$gene}{$weight}{1} > 0 && $grouplength{$gene}{$weight}{0} > 0) {
	    $status = "both";
	} else {
	    $status = "one";
	}
	my $tmph1 = $grouplength{$gene}{$weight}{1} / $genelength{$gene};
	my $tmph0 = $grouplength{$gene}{$weight}{0} / $genelength{$gene};

	if($status eq "both") {
	    print "$gene\t$weight\t$status\t+\t$tmph1\t$isoindex_complement{$gene}{$weight}{1}\t$exonset{$gene}{$weight}{1}\n";
	    print "$gene\t$weight\t$status\t-\t$tmph0\t$isoindex_complement{$gene}{$weight}{0}\t$exonset{$gene}{$weight}{0}\n";
	} 
	
	if($status eq "one") {
	    if($grouplength{$gene}{$weight}{1} > 0) {
		print "$gene\t$weight\tplus\t+\t$tmph1\t$isoindex_complement{$gene}{$weight}{1}\t$exonset{$gene}{$weight}{1}\n";
	    } else {
		print "$gene\t$weight\tminus\t+\t$tmph0\t$isoindex{$gene}{$weight}{0}\tNA\n";
	    }
	    if($grouplength{$gene}{$weight}{0} > 0) {
		print "$gene\t$weight\tminus\t-\t$tmph0\t$isoindex_complement{$gene}{$weight}{0}\t$exonset{$gene}{$weight}{0}\n";
	    } else {
		print "$gene\t$weight\tplus\t-\t$tmph1\t$isoindex{$gene}{$weight}{1}\tNA\n";
	    }
	}
    }
}
