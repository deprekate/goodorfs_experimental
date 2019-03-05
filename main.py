import sys
import time
import argparse
import numpy as np
import pandas as pd 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.patches as mpatches

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
parser.add_argument('-t','--clust_type', action="store", dest="clust_type", required=True)
parser.add_argument('-n','--clust_num', action="store", dest="clust_num", required=True)
args = parser.parse_args()
print(args.genome_id)
exit()

#-----------------------READ IN THE DATA
dat = pd.read_csv("aa/" + sys.argv[1] + ".aa", header=0,sep=' ' )
end = pd.read_csv("prodigal/" + sys.argv[1] + ".ends", header=None)

dat['TYPE'] = dat['STOP'].isin(end[0])

tag = sys.argv[1].split("/")[-1].split(".")[0]
#X = dat.groupby('STOP')[list("ARNDCEQGHILKMFPSTWYV")].sum()

X = []
X = dat[list("ARNDCEQGHILKMFPSTWYV")]
X = X.div(X.sum(axis=1), axis=0)
if(sys.argv[2] == 'se'):
#this is the shannon entropy
	X = X * np.log(X)
	X = X.div(X.sum(axis=1), axis=0)
	X.fillna(0, inplace=True)
	title = 'shannon-entropy'
else:
	title = 'aminoacid-percent'
X_std = StandardScaler().fit_transform(X)


pca = decomposition.PCA(n_components=2)
pca.fit(X_std)
P = pca.transform(X_std)


if(sys.argv[3] == 'km'):
	gm = KMeans(n_clusters=4).fit(X_std).labels_
#DO THE MIXTUREMODEL
#gm = GaussianMixture(n_components=2,covariance_type='tied').fit(X_std).predict(X_std)
#gm = BayesianGaussianMixture(n_components=3,covariance_type='tied').fit_predict(X_std)
#gm = AgglomerativeClustering(n_clusters=2).fit_predict(X_std)

fig = plt.figure()
fig.suptitle(tag + " (" + title + ")", fontsize=20)


if(sys.argv[3] == 'annotate'):
	for i,p in enumerate(P):
		print(dat.TYPE[i])
		plt.scatter(p[0], p[1], c=color(dat.TYPE[i]), alpha=0.5, edgecolors='none', label='fake')
		plt.annotate(str(dat.STOP[i]), (p[0], p[1]), size=8)
		#print(dat.STOP[i])
		#plt.pause(0.2)
	red_patch = mpatches.Patch(color='red', label='coding')
	blue_patch = mpatches.Patch(color='blue', label='noncoding')
	plt.legend(handles=[red_patch, blue_patch])
	plt.show()
	exit()

x0 = X_std[gm==0]
x1 = X_std[gm==1]
x2 = X_std[gm==2]
variance = [np.var(x0),np.var(x1),np.var(x2)]
variance = [v for v in variance if str(v) != 'nan']
index_min = np.argmin(variance)
index_mid = np.argmed(variance)
index_max = np.argmax(variance)
x0, x1, x2 = X_std[gm==index_min], X_std[gm==index_mid], X_std[gm==index_max]
p0 = P[gm==index_min]
p1 = P[gm==index_mid]
p2 = P[gm==index_max]

#X['STOP'] = X.index
#X['CLUSTER'] = gm
#X['TYPE'] = X.index.isin(end[0])
#dat = X

dat['CLUSTER'] = gm
dat['GOOD'] = dat.STOP.map(dat[dat.CLUSTER==index_min].STOP.value_counts())
dat['TOTAL'] = dat.STOP.map(dat.STOP.value_counts())
#dat.fillna(0, inplace=True)
#print(dat[gm==index_min].query('TYPE == "real"').TYPE.count(), end = "\t" )
#print(dat[gm==index_min].TYPE.count(), end = "\t" )
#print(dat[gm==index_mid].query('TYPE == "real"').TYPE.count(), end = "\t" )
#print(dat[gm==index_mid].TYPE.count(), end = "\t" )
#print(dat[gm==index_max].query('TYPE == "real"').TYPE.count(), end = "\t" )
#print(dat[gm==index_max].TYPE.count(), end = "\n" )
#exit()
#df = dat.drop_duplicates('STOP')
#for index, row in df.iterrows():
#	print(sys.argv[1], row.STOP, row.TYPE, row.GOOD/row.TOTAL)

#exit()
#TP = len(dat[(dat.CLUSTER == index_min) & (dat.TYPE) & (dat.GOOD/dat.TOTAL > 0.5)].STOP.unique() )
#FP = len(dat[(dat.CLUSTER == index_min) & (~dat.TYPE) & (dat.GOOD/dat.TOTAL > 0.5)].STOP.unique() )
#TN = len(dat[~dat.TYPE].STOP.unique()) - FP
#FN = len(dat[dat.TYPE].STOP.unique()) - TP
#print(sys.argv[1], TP, FP, FN, TN)
#exit()
p = P[dat.TYPE]
plt.scatter(p[:, 0], p[:, 1], c='g', alpha=0.5, edgecolors='none', label='real')
p = P[~dat.TYPE]
plt.scatter(p[:, 0], p[:, 1], c='r', alpha=0.5, edgecolors='none', label='fake')
plt.scatter(p0[:, 0], p0[:, 1], facecolor='none', cmap='Spectral', alpha=0.5, marker='d', edgecolor='black', label='Predicted')

#print('             ', 'MAD', 'σ', 'AAD', sep='\t')
print('black_diamond: ', round(np.mad(x0),4), round(np.var(x0),4), round(np.aad(x0),4), sep='\t')
plt.scatter(p1[:, 0], p1[:, 1], facecolor='none', cmap='Spectral', alpha=0.5, marker='D', edgecolor='b')
print(' blue_rectang: ', round(np.mad(x1),4), round(np.var(x1),4), round(np.aad(x1),4), sep='\t')
plt.scatter(p2[:, 0], p2[:, 1], facecolor='none', cmap='Spectral', alpha=0.5, marker='s', edgecolor='g')
print('green_square:', round(np.mad(x2),4), round(np.var(x2),4), round(np.aad(x2),4), sep='\t')

plt.legend()
#fig.savefig('fig/se_km/' + tag + '_se_km.png', dpi=150)

plt.show() #block=False)
#time.sleep(5)
#plt.close("all") 
exit()







