#!/usr/bin/perl  -w

open (IN,"<$ARGV[0]");
open (OUT,">$ARGV[1]");

while($line=<IN>){
	@a=split(/\t/,$line);
	$tax=$a[13];
	if($tax=~/(.+?)\>/)
	{
		$name=$1;
	}
	else{
		$name=$tax;
	}
	print OUT "$a[0]\t$a[1]\t$a[2]\t$a[3]\t$a[4]\t$a[5]\t$a[6]\t$a[7]\t$a[8]\t$a[9]\t$a[10]\t$a[11]\t$a[12]\t$name\t$a[14]";
}
close OUT;
close IN;
