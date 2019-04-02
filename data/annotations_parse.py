import sys
import os
import gzip


def genemark(file_contents):
	for line in file_contents.splitlines():
		split_line = line.split()
		if len(split_line) == 6 and split_line[1] == '+':
			yield split_line[3]
		elif len(split_line) == 6 and split_line[1] == '-':
			yield split_line[2]
	return
	
def glimmer(file_contents):
	for line in file_contents.splitlines():
		split_line = line.split()
		if len(split_line) == 5 and split_line[3][0] in ['+', '-']:
			yield split_line[2]
	return
	
def phanotate(file_contents):
	for line in file_contents.splitlines():
		split_line = line.split()
		if len(split_line) == 5 and split_line[2] in ['+', '-']:
			yield split_line[1]
	return
	
def prodigal(file_contents):
	for line in file_contents.splitlines():
		split_line = line.split('_')
		if len(split_line) == 4 and split_line[3] == '+':
			yield split_line[2]
		elif len(split_line) == 4 and split_line[3] == '-':
			yield split_line[1]
	return

all_orf_ends = dict()
for output_type in ['genemark', 'glimmer', 'phanotate', 'prodigal']:
	with gzip.open('annotations/' + sys.argv[1] + '.' + output_type + '.gz', 'rb') as fp:
		contents = fp.read() # contents now has the uncompressed bytes of foo.gz
		fp.close()
		u_str = contents.decode('utf-8')
		orf_ends = globals()[output_type](u_str)
		for orf_end in orf_ends:
			orf_end = orf_end.replace('<','')
			orf_end = orf_end.replace('>','')
			orf_count = all_orf_ends.get(orf_end, 0)
			orf_count += 1
			all_orf_ends[orf_end] = orf_count
for end in all_orf_ends:
	if all_orf_ends[end] > 1:
		print(end)


