import os
import sys
import time
import argparse
import numpy as np
import pandas as pd 
import gzip

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
path = os.path.dirname(os.path.realpath(__file__))
dat = pd.read_csv(path + "/data/aa/" + args.genome_id + ".tsv.gz", compression='gzip', header=0,sep='\t' )

with gzip.open(path + '/genomes/fna/' + args.genome_id + '.fna.gz') as f:
	first_line = f.readline().decode("utf-8")

X = dat[list("ARNDCEQGHILKMFPSTWYV")]
X = X.div(X.sum(axis=1), axis=0)

dat.loc[:,'A':'V'] = dat.loc[:,'A':'V'].div(dat.loc[:,'A':'V'].sum(axis=1), axis=0)

if(args.data_type == 'se'):
#this is the shannon entropy
	X = X * np.log(X)
	X = X.div(X.sum(axis=1), axis=0)
	X.fillna(0, inplace=True)
	title = 'shannon-entropy'
else:
	title = 'aminoacid-percent'

X = StandardScaler().fit_transform(X)

if(args.clust_type == 'km'):
	dat['CLUSTER'] = KMeans(n_clusters=args.clust_num, n_init=50).fit(X).labels_
elif(args.clust_type == 'gm'):
	dat['CLUSTER'] = GaussianMixture(n_components=2,covariance_type='tied').fit(X).predict(X)
elif(args.clust_type == 'bgm'):
	dat['CLUSTER'] = BayesianGaussianMixture(n_components=3,covariance_type='tied').fit_predict(X)
elif(args.clust_type == 'ag'):
	dat['CLUSTER'] = AgglomerativeClustering(n_clusters=2).fit_predict(X)

pca = decomposition.PCA(n_components=2).fit(X)
dat['x'],dat['y'] = pca.fit_transform(X).T

# PICK THE CLUSTER WITH THE LOWEST VARIANCE
variance = [np.var(X[dat.CLUSTER==i]) for i in range(3)]
index_min = np.argmin(variance)
#index_mid = np.argmed(variance)
#index_max = np.argmax(variance)

dat['GOOD'] = dat.STOP.map(dat[dat.CLUSTER==index_min].STOP.value_counts())
dat['TOTAL'] = dat.STOP.map(dat.STOP.value_counts())
dat.fillna(0, inplace=True)

#########THIS MAKES THE VIOLIN DATA###########
#df = dat.drop_duplicates('STOP')
#for index, row in df.iterrows():
#	print(args.genome_id, row.STOP, row.TYPE, row.GOOD/row.TOTAL)
#exit()
#########THIS MAKES THE BAR DATA###########
#TP = len(dat[(dat.CLUSTER == index_min) & (dat.TYPE) & (dat.GOOD/dat.TOTAL > 0.5)].STOP.unique() )
#FP = len(dat[(dat.CLUSTER == index_min) & (~dat.TYPE) & (dat.GOOD/dat.TOTAL > 0.5)].STOP.unique() )
#TN = len(dat[~dat.TYPE].STOP.unique()) - FP
#FN = len(dat[dat.TYPE].STOP.unique()) - TP
#print(args.genome_id, round(np.var(x0),4), round(np.var(x1),4), round(np.var(x2),4), TP, FP, FN, TN)
#exit()

#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
fig, ax = plt.subplots()

# PLOT
colors = {True:'#3CC9CF', False:'#F2766E'}

fig.suptitle(args.genome_id + " (" + title + ")")
ax.scatter(dat['x'], dat['y'], c=dat['TYPE'].apply(lambda x: colors[x]), marker='.', linewidths=0.0, alpha=0.8, zorder=5)

ax.scatter(dat[dat.CLUSTER==index_min].x, dat[dat.CLUSTER==index_min].y, facecolor='none', cmap='Spectral', alpha=0.5, marker='d', edgecolor='black', label='Predicted', zorder=10)

if(args.annotate):
	p = dat[dat.CLUSTER==index_min]
	for i, row in p.iterrows():
		plt.annotate(str(row.STOP), (row.x, row.y), size=4)

plt.legend(handles=[mpatches.Patch(color=col, label=str(lab)) for lab,col in colors.items()])

fig.set_size_inches(5, 5)
fig.savefig(args.genome_id + '.png', dpi=200)
#plt.show() #block=False)
#time.sleep(4)
#plt.close("all") 
exit()







