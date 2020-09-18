import os
import sys
import time
import argparse
import numpy as np
import pandas as pd 
import gzip
from math import sqrt
from statistics import mean
from math import log10
from math import log2

#from matplotlib.mlab import PCA
from sklearn.decomposition import PCA
from sklearn.preprocessing import scale
from sklearn.preprocessing import StandardScaler

from sklearn import decomposition
from sklearn import datasets
from sklearn.mixture import GaussianMixture
from sklearn.mixture import BayesianGaussianMixture
from sklearn.cluster import KMeans
from sklearn.cluster import AgglomerativeClustering
from sklearn_extra.cluster import KMedoids
from sklearn.cluster import OPTICS
from sklearn.cluster import SpectralClustering
from sklearn.cluster import Birch
from scipy import stats

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
def amd(data, axis=None):
	    return np.median(np.absolute(data - np.median(data, axis)), axis)
setattr(np, 'amd', amd)

def argmid(x):
    for i, item in enumerate(x):
        if item != min(x) and item != max(x):
            return i
setattr(np, 'argmid', argmid)


def color(x):
	if(x):
		return 'red'
	else:
		return 'blue'

def get_best(variances, counts):
	index = None
	v = float("+Inf")
	for i, (var,count) in enumerate(zip(variances, counts)):
		if var<v and var>=0 and (count/sum(counts))>0.05:
			index = i
			v = var
	return index
			




#################################
parser = argparse.ArgumentParser()
parser.add_argument('-i','--genome_id', action="store", dest="genome_id", required=True)
parser.add_argument('-d','--data_type', action="store", dest="data_type", required=True)
parser.add_argument('-t','--clust_type', action="store", dest="clust_type", required=True)
parser.add_argument('-n','--clust_num', action="store", dest="clust_num", required=True, type=int)
parser.add_argument('-a','--annotate', action="store_true", dest="annotate", required=False)
args = parser.parse_args()

#-----------------------READ IN THE DATA
path = os.path.dirname(os.path.realpath(__file__))
dat = pd.read_csv(path + "/data/aa/" + args.genome_id + ".tsv.gz", compression='gzip', header=0,sep='\t' )

with gzip.open(path + '/genomes/fna/' + args.genome_id + '.fna.gz') as f:
	first_line = f.readline().decode("utf-8")

counts = dat.groupby('STOP').START.nunique().to_dict()

count = dict()
mask = []
longest = dict()
for index,row in dat[['START', 'STOP']].iterrows():
	r = list(row.values)
	count[r[1]] = count.get(r[1], 0) + 1
	#if count[r[1]] <= 30:
	#if count[r[1]] <= 4*log2(counts[r[1]]):
	if count[r[1]] <= 10+sqrt(counts[r[1]]):

	#if count[r[1]] <= 10: #(counts[r[0]]/2):
		mask.append(False)
	else:
		mask.append(True)
	length = max(r[1]-r[0], r[0]-r[1])
	longest[r[1]] = max(longest.get(r[1],0), length) 

dat = dat.drop(dat[mask].index)
#dat = dat[dat.CODON != 'TTG']

nums = []
for index,row in dat[['START','STOP']].iterrows():
	r = list(row.values)
	length = max(r[1]-r[0], r[0]-r[1])
	#nums.append(r[0] / counts[r[1]])
	nums.append(length/longest[r[1]])
	#nums.append(length)


X = dat[list("ARNDCEQGHILKMFPSTWYV")]
X = X.div(X.sum(axis=1), axis=0)
#X.loc[:, list("ARNDCEQGHILKMFPSTWYV")]  = StandardScaler().fit_transform(X)

dat.loc[:,'A':'V'] = dat.loc[:,'A':'V'].div(dat.loc[:,'A':'V'].sum(axis=1), axis=0)

if(args.data_type == 'se'):
	#this is the shannon entropy
	X = X * np.log(X)
	X.fillna(0, inplace=True)
	#X = np.nan_to_num(X, nan=0)
	X = X.div(X.sum(axis=1), axis=0)
	#X = X / X.sum(axis=1, keepdims=True)
	title = 'shannon-entropy'
else:
	X = StandardScaler().fit_transform(X)
	title = 'aminoacid-percent'

#X['NUM'] = dat['NUM']
#X['NUM'] = dat['NUM'].div(dat['NUM'].mean()/X.to_numpy().mean(), axis=0
#X['NUM'] = np.array(nums) * -np.log(np.array(nums))
#X['NUM'] = nums
#X = StandardScaler().fit_transform(X)

args.clust_num = int(log10(len(X)))

print('len', len(X), 'clustnum', args.clust_num) 


if(args.clust_type == 'km'):
	dat['CLUSTER'] = KMeans(n_clusters=args.clust_num, n_init=50).fit(X).labels_
elif(args.clust_type == 'gm'):
	#dat['CLUSTER'] = GaussianMixture(n_components=args.clust_num,covariance_type='spherical').fit(X).predict(X)
	model = GaussianMixture(n_components=args.clust_num, n_init=10, covariance_type='spherical', reg_covar=0.01).fit(X)
	print(model.covariances_)
	print(model.converged_)
	print(model.n_iter_)
	print(model.lower_bound_)
	dat['CLUSTER'] = model.predict(X)
elif(args.clust_type == 'bgm'):
	dat['CLUSTER'] = BayesianGaussianMixture(n_components=args.clust_num,covariance_type='spherical').fit_predict(X)
elif(args.clust_type == 'ag'):
	dat['CLUSTER'] = AgglomerativeClustering(n_clusters=args.clust_num).fit_predict(X)
elif(args.clust_type == 'ward'):
	dat['CLUSTER'] = AgglomerativeClustering(n_clusters=args.clust_num, linkage='ward').fit_predict(X)
elif(args.clust_type == 'pam'):
	dat['CLUSTER'] = KMedoids(n_clusters=args.clust_num).fit_predict(X)
elif(args.clust_type == 'op'):
	dat['CLUSTER'] = OPTICS(min_samples=10).fit_predict(X)
elif(args.clust_type == 'sp'):
	dat['CLUSTER'] = SpectralClustering(n_clusters=args.clust_num, assign_labels="discretize").fit_predict(X)
elif(args.clust_type == 'b'):
	dat['CLUSTER'] = Birch(n_clusters=args.clust_num).fit_predict(X)


pca = decomposition.PCA(n_components=2).fit(X)
dat['x'],dat['y'] = pca.fit_transform(X).T

# PICK THE CLUSTER WITH THE LOWEST VARIANCE
variance = [np.sum(np.var(X[dat.CLUSTER==i], axis=0)) for i in range(args.clust_num)]
print("----------")
print('   var', variance)
v = [np.sum(np.aad(X[dat.CLUSTER==i], axis=0)) for i in range(args.clust_num)]
print('sumaad', v)
v = [np.sum(np.amd(X[dat.CLUSTER==i], axis=0)) for i in range(args.clust_num)]
print('   amd', v)
l = [len(X[dat.CLUSTER==i]) for i in range(args.clust_num)]
print('   len', l)
#index_min = np.argmin(variance)
#index_min = np.argmin(v)
index_min = get_best(v, l)
print(index_min)
#index_mid = np.argmed(variance)
#index_max = np.argmax(variance)
title = "[" + ', '.join(str(round(e, 2)) for e in v) + "]"

dat['GOOD'] = dat.STOP.map(dat[dat.CLUSTER==index_min].STOP.value_counts())
dat['TOTAL'] = dat.STOP.map(dat.STOP.value_counts())
dat.fillna(0, inplace=True)

#print(dat.to_string())
#exit()
#########THIS MAKES THE VIOLIN DATA###########
#df = dat.drop_duplicates('STOP')
#for index, row in df.iterrows():
#	print(args.genome_id, row.STOP, row.TYPE, row.GOOD/row.TOTAL)
#exit()
#########THIS MAKES THE BAR DATA###########
#sd
#TP = len(dat[(dat.CLUSTER == index_min) & (dat.TYPE) & (dat.GOOD/dat.TOTAL > 0.5)].STOP.unique() )
#FP = len(dat[(dat.CLUSTER == index_min) & (~dat.TYPE) & (dat.GOOD/dat.TOTAL > 0.5)].STOP.unique() )
#TN = len(dat[~dat.TYPE].STOP.unique()) - FP
#FN = len(dat[dat.TYPE].STOP.unique()) - TP
#print(args.genome_id, round(np.var(x0),4), round(np.var(x1),4), round(np.var(x2),4), TP, FP, FN, TN)
#exit()

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
fig, ax = plt.subplots()

# PLOT
colors = {True:'#3CC9CF', False:'#F2766E'}
markers = {0:'o', 1:'v', 2:'^', 3:'<', 4:'>', 5:'s'}

fig.suptitle(args.genome_id + " (" + title + ")")
ax.scatter(dat['x'], dat['y'], c=dat['TYPE'].apply(lambda x: colors[x]), marker='.', linewidths=0.0, alpha=0.4, zorder=5)


if(args.annotate):
	p = dat[dat.CLUSTER==index_min]
	for i, row in p.iterrows():
		plt.annotate(str(row.STOP), (row.x, row.y), size=4, zorder=10)
		print(row.STOP)
else:
	ax.scatter(dat[dat.CLUSTER==index_min].x, dat[dat.CLUSTER==index_min].y, facecolor='none', cmap='Spectral', s=80, linewidths=0.3, alpha=0.3, marker='d', edgecolor='black', label='Predicted', zorder=10)
	for i in range(args.clust_num):
		if i != index_min:
			ax.scatter(dat[dat.CLUSTER==i].x, dat[dat.CLUSTER==i].y, facecolor='none', cmap='Spectral', s=80, linewidths=0.3, alpha=0.5, marker=markers[i], edgecolor='green', label='Predicted', zorder=10)
	pass

## visualize projections
#xvector = pca.components_[0] # see 'prcomp(my_data)$rotation' in R
#yvector = pca.components_[1]
#xs = pca.transform(X)[:,0] # see 'prcomp(my_data)$x' in R
#ys = pca.transform(X)[:,1]
#for i in range(len(xvector)):
	# arrows project features (ie columns from csv) as vectors onto PC axes
	#plt.arrow(0, 0, xvector[i]*max(xs), yvector[i]*max(ys), color='black', alpha=0.7, width=0.0005, head_width=0.0025, zorder=11)
	#plt.text(xvector[i]*max(xs)*1.2, yvector[i]*max(ys)*1.2, list("ARNDCEQGHILKMFPSTWYV")[i], color='black', alpha=0.7, zorder=11)
#	pass


plt.legend(handles=[mpatches.Patch(color=col, label=str(lab)) for lab,col in colors.items()])

fig.set_size_inches(4, 2.5)
fig.savefig(args.genome_id + '.png', dpi=300)
plt.show() #block=False)
#time.sleep(4)
#plt.close("all") 
exit()







