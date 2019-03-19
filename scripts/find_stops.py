import sys
import os

def get_stop(type, line):
	if type == 'genemark':
		cols = line.split()
		if len(cols) == 6 and cols[1] == '+':
			return cols[3].replace('>','')
		elif len(cols) == 6 and cols[1] == '-':
			return cols[2].replace('<','')
	elif type == 'glimmer':
		
	return


my_stops = dict()
for type in ['genemark','glimmer','phanotate','prodigal']:
	filename = sys.argv[1] + "." + type
	file = open('../gbk/' + filename, "r")
	for line in file:
		my_stop = get_stop(type, line)
		if my_stop:
			my_types = my_stops.get(my_stop, [])
			my_types.append(type)
			my_stops[my_stop] = my_types

print(my_stops)

		
