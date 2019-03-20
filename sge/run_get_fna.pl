#!/usr/bin/perl

my $flag = 0;
open(INFILE, $ARGV[0]) or die();
while(<INFILE>){	
	chomp();
	if(m/^LOCUS/){
		print ">";
		s/LOCUS *//;
		s/ .*/ /;
		print;
	}
	if(m/^DEFINITION/){
		s/DEFINITION *//;
		print;
		print "\n";
	}
	
	next unless (m/^ORIGIN/ or $flag);	
	$flag = 1;
	s/[^A-z]//g;
	print;
	print "\n";

}
