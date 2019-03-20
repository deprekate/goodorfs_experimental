#!/usr/bin/perl

my $ids = `cat ../list.tsv`;
my @ids = split(/\n/, $ids);

foreach my $id (@ids){
	`~/opt/genemark/gmsn.pl --phage ../fna/$id.fna`;
}
