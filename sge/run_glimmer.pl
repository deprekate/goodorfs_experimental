#!/usr/bin/perl

my $ids = `cat list.tsv`;
my @ids = split(/\n/, $ids);

foreach my $id (@ids){
	`/home3/katelyn/opt/glimmer3.02/scripts/g3-from-scratch.csh fna/$id.fna $id`;
	unlink("$id.icm");
	unlink("$id.longorfs");
	unlink("$id.train");
	unlink("$id.detail");
}
