#!/usr/bin/perl  -w
#use strict;

open (LST,"<$ARGV[0]"); #open contig id list file
open (OUT,">$ARGV[2]"); #open for write into subseq file
while($contig_id=<LST>){
	$contig_id=~s/[\r\n]$//g;
	open (FASTA,"<$ARGV[1]"); #open fasta format seq file
	$/=">";
	<FASTA>;
	while($seq=<FASTA>){
		chomp $seq;
		if($seq=~/^(\S+)/){
			$id=$1;
			if($id=~/$contig_id/){
					print OUT ">$seq";
					last;
				}
			}
			
		}
	close FASTA;
	$/="\n";
}
close LST;
close OUT;









