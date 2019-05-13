import os
import sys
import time
import argparse
import numpy as np
import pandas as pd 
import gzip

import warnings
warnings.filterwarnings('ignore')

#################################
parser = argparse.ArgumentParser()
parser.add_argument('-i','--genome_id', action="store", dest="genome_id", required=True)
args = parser.parse_args()

def get_order(row):
	return min(row[['START', 'STOP']]), max(row[['START','STOP']]), (row['STOP']-row['START'])/abs(row['STOP']-row['START'])
def max_idx(a, b, c):
	if(a > b):
		if(a > c):
			return 1;
		else:
			return 3;
	else:
		if(b > c):
			return 2;
		else:
			return 3;
def most_common(lst):
	if not lst:
		return None 
	else:
    		return max(set(lst), key=lst.count)

#-----------------------READ IN THE DATA
path = '..' #os.path.dirname(os.path.realpath(__file__))
dat = pd.read_csv(path + "/data/aa/" + args.genome_id + ".tsv.gz", compression='gzip', header=0,sep='\t' )
end = pd.read_csv(path + "/data/ends/" + args.genome_id + ".tsv", header=None)

with gzip.open(path + '/genomes/fna/' + args.genome_id + '.fna.gz') as f:
	first_line = f.readline().decode("utf-8")

dat['TYPE'] = dat['STOP'].isin(end[0])
dat['OFFSET'] = None

bases = [[] for i in range(3+dat.loc[:, ['START', 'STOP']].max().max())]
for index, row in dat.iterrows():
	if row['TYPE']:
		left, right, direction = get_order(row)
		for i in range(left, right,3):
			bases[i+0].append( direction * (1+(i-1)%3) )
			bases[i+1].append( direction * (1+(i-1)%3) )
			bases[i+2].append( direction * (1+(i-1)%3) )
		dat.at[index,'OFFSET'] = 0

for index, row in dat.iterrows():
	if not row['TYPE']:
		left, right, direction = get_order(row)
		coding_frames = []
		for i in range(left, right+3, 3):
			coding_frames.append(most_common(bases[i]))
		coding_frame = most_common(coding_frames)
		offset_frame = direction * (1+(i-1)%3)
		if coding_frame is not None:
			dat.at[index,'OFFSET'] = direction*(coding_frame/abs(coding_frame))*(((3+abs(offset_frame)) - abs(coding_frame))%3)

#direction * (1+(((1+(i-1)%3)+
#most_common(coding_frames))%3))

print(dat.to_string())
			
			





