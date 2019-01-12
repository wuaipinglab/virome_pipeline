#!/usr/bin/perl  -w

open (IN,"<$ARGV[0]");
open (OUT,">$ARGV[2]") ;

$/=">";
<IN>;
while($seq=<IN>){
	if($seq=~/^(\S+)/){
		$id=$1;
	}
	$seq=~s/.+\n//;
    $seq=~s/\s//g;
	chomp $seq;
	$len=length($seq);
	print OUT ">$ARGV[1]_$id len=$len\n$seq\n"
	
}
$/="\n";
close IN;
close OUT;









