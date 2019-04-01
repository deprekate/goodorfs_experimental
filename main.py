import sys
import time
import argparse
import numpy as np
import pandas as pd 

from matplotlib.mlab import PCA
from sklearn.preprocessing import StandardScaler

from sklearn import decomposition
from sklearn import datasets
from sklearn.mixture import GaussianMixture
from sklearn.mixture import BayesianGaussianMixture
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering

import warnings
warnings.filterwarnings('ignore')

def argmed(a):
    if(len(a) == 2):
	    return -1
    if len(a) % 2 == 1:
        return np.where( a == np.median(a) )[0][0]
    else:
        l,r = len(a)/2 -1, len(a)/2
        left = np.partition(a, l)[l]
        right = np.partition(a, r)[r]
        return [np.where(a == left)[0][0], np.where(a==right)[0][0]]
setattr(np, 'argmed', argmed)

def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation 
    """
    arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.median(arr)
    return np.median(np.abs(arr - med))
setattr(np, 'mad', mad)
def aad(data, axis=None):
	    return np.mean(np.absolute(data - np.mean(data, axis)), axis)
setattr(np, 'aad', aad)
def color(x):
	if(x):
		return 'red'
	else:
		return 'blue'
#################################
parser = argparse.ArgumentParser()
parser.add_argument('-i','--genome_id', action="store", dest="genome_id", required=True)
parser.add_argument('-d','--data_type', action="store", dest="data_type", required=True)
parser.add_argument('-t','--clust_type', action="store", dest="clust_type", required=True)
parser.add_argument('-n','--clust_num', action="store", dest="clust_num", required=True, type=int)
parser.add_argument('-a','--annotate', action="store_true", dest="annotate", required=False)
args = parser.parse_args()

#-----------------------READ IN THE DATA
dat = pd.read_csv("/home3/katelyn/goodorfs/data/aa/" + args.genome_id + ".tsv", header=0,sep='\t' )
end = pd.read_csv("/home3/katelyn/goodorfs/data/ends/prodigal/" + args.genome_id + ".prod", header=None)

with open('/home3/katelyn/goodorfs/genomes/fna/' + args.genome_id + '.fna') as f:
	first_line = f.readline()

dat = dat[dat.CODON != 'CTG']
dat['COUNT'] = dat.groupby(['STOP'])['START'].transform('count')

import re
name = re.search('\[.*\]', first_line).group(0)

dat['TYPE'] = dat['STOP'].isin(end[0])

#dat = dat.groupby(['STOP','TYPE'],as_index=False)[list("ARNDCEQGHILKMFPSTWYV")].sum()

#dat = dat[dat['CODON'].isin(['ATG','GTG','TTG'])]

X = dat[list("ARNDCEQGHILKMFPSTWYV")]
X = X.div(X.sum(axis=1), axis=0)
if(args.data_type == 'se'):
#this is the shannon entropy
	X = X * np.log(X)
	X = X.div(X.sum(axis=1), axis=0)
	X.fillna(0, inplace=True)
	title = 'shannon-entropy'
else:
	title = 'aminoacid-percent'
X['n'] = dat['COUNT']
X_std = StandardScaler().fit_transform(X)

if(args.clust_type == 'km'):
	dat['CLUSTER'] = KMeans(n_clusters=args.clust_num, n_init=50).fit(X_std).labels_
elif(args.clust_type == 'gm'):
	dat['CLUSTER'] = GaussianMixture(n_components=2,covariance_type='tied').fit(X_std).predict(X_std)
elif(args.clust_type == 'bgm'):
	dat['CLUSTER'] = BayesianGaussianMixture(n_components=3,covariance_type='tied').fit_predict(X_std)
elif(args.clust_type == 'ag'):
	dat['CLUSTER'] = AgglomerativeClustering(n_clusters=2).fit_predict(X_std)

pca = decomposition.PCA(n_components=2)
pca.fit(X_std)
P = pca.transform(X_std)
dat['x'] = P[:,0]
dat['y'] = P[:,1]


x0 = X_std[dat.CLUSTER==0]
x1 = X_std[dat.CLUSTER==1]
x2 = X_std[dat.CLUSTER==2]
variance = [np.var(x0),np.var(x1),np.var(x2)]
variance = [v for v in variance if str(v) != 'nan']
index_min = np.argmin(variance)
index_mid = np.argmed(variance)
index_max = np.argmax(variance)
x0, x1, x2 = X_std[dat.CLUSTER==index_min], X_std[dat.CLUSTER==index_mid], X_std[dat.CLUSTER==index_max]
p0 = P[dat.CLUSTER==index_min]
p1 = P[dat.CLUSTER==index_mid]
p2 = P[dat.CLUSTER==index_max]

dat['GOOD'] = dat.STOP.map(dat[dat.CLUSTER==index_min].STOP.value_counts())
dat['TOTAL'] = dat.STOP.map(dat.STOP.value_counts())
dat.fillna(0, inplace=True)

#########THIS MAKES THE VIOLIN DATA###########
#df = dat.drop_duplicates('STOP')
#for index, row in df.iterrows():
#	print(args.genome_id, row.STOP, row.TYPE, row.GOOD/row.TOTAL)
#exit()
#########THIS MAKES THE BAR DATA###########
TP = len(dat[(dat.CLUSTER == index_min) & (dat.TYPE) & (dat.GOOD/dat.TOTAL > 0.5)].STOP.unique() )
FP = len(dat[(dat.CLUSTER == index_min) & (~dat.TYPE) & (dat.GOOD/dat.TOTAL > 0.5)].STOP.unique() )
TN = len(dat[~dat.TYPE].STOP.unique()) - FP
FN = len(dat[dat.TYPE].STOP.unique()) - TP
print(args.genome_id, round(np.var(x0),4), round(np.var(x1),4), round(np.var(x2),4), TP, FP, FN, TN)
exit()

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig = plt.figure()

#fig.suptitle(args.genome_id + " (" + title + ")", fontsize=20)
fig.suptitle(name, fontsize=20)
p = P[~dat.TYPE]
plt.scatter(p[:, 0], p[:, 1], c='r', alpha=0.5, edgecolors='none', label='fake')
p = P[dat.TYPE]
plt.scatter(p[:, 0], p[:, 1], c='g', alpha=0.5, edgecolors='none', label='real')
plt.scatter(p0[:, 0], p0[:, 1], facecolor='none', cmap='Spectral', alpha=0.5, marker='d', edgecolor='black', label='Predicted')
if(args.annotate):
	p = dat[dat.CLUSTER==index_min]
	for i, row in p.iterrows():
		plt.annotate(str(row.STOP), (row.x, row.y), size=4)

#print('             ', 'MAD', 'σ', 'AAD', sep='\t')
#print('black_diamond: ', round(np.mad(x0),4), round(np.var(x0),4), round(np.aad(x0),4), sep='\t')
#plt.scatter(p1[:, 0], p1[:, 1], facecolor='none', cmap='Spectral', alpha=0.5, marker='D', edgecolor='b')
#print(' blue_rectang: ', round(np.mad(x1),4), round(np.var(x1),4), round(np.aad(x1),4), sep='\t')
#plt.scatter(p2[:, 0], p2[:, 1], facecolor='none', cmap='Spectral', alpha=0.5, marker='s', edgecolor='g')
#print('green_square:', round(np.mad(x2),4), round(np.var(x2),4), round(np.aad(x2),4), sep='\t')

plt.legend()
#fig.set_size_inches(20, 10)
fig.savefig(args.genome_id + '.png', dpi=100)
#plt.show(block=False)
#time.sleep(4)
#plt.close("all") 
exit()







