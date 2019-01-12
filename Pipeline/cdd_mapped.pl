#!/usr/bin/perl  -w
#use strict;

open (LST,"<$ARGV[0]") 
open (OUT,">$ARGV[2]") ;
while($contig_id=<LST>){
	$contig_id=~s/[\r\n]$//g;
	open (FASTA,"<$ARGV[1]") ;
	while($seq=<FASTA>){
		@a=split/\t/,$seq;
		$name=$a[1];
		$name=~s/[\r\n]$//g;
		if($name=~/$contig_id/){
			print OUT "$seq";
		}
	}
	close FASTA;
	
}
close LST;
close OUT;









