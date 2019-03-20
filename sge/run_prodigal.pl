#!/usr/bin/perl

#my $ids = `cat list.tsv`;
#my @ids = split(/\n/, $ids);

my $id = $ARGV[0];

#foreach my $id (@ids){
my $out = `/home3/katelyn/opt/Prodigal-edit/prodigal -q -i fna/$id.fna -f sco -o prodigal/$id.prodigal 2>&1 1>/dev/null`;
if($out){
	my $out = `/home3/katelyn/opt/Prodigal/prodigal -q -p meta -i fna/$id.fna -f sco -o prodigal/$id.prodigal 2>&1 1>/dev/null`;
}
