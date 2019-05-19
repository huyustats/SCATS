open FP, "list";
while(<FP>) {
    chomp();
    my $cmd = "samtools view -bT ~/mm10_fasta/mm10_nochr.fa $_\.sam > $_\.bam";
    #system "echo \"$cmd\" > stb_$_\.sh ";
    #system "bsub sh stb_$_\.sh";
    system "samtools index $_\.bam";
}
close FP;
