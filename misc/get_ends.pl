#!/usr/bin/perl

open(INFILE, $ARGV[0]) or die();
while(<INFILE>){
	my @line = split();
	if($line[6] eq '+'){
		print $line[4], "\n";
	}elsif($line[6] eq '-'){
		print $line[3], "\n";
	}

}

